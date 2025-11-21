#ifndef SERIALTOMULTICONVERSION_H
#define SERIALTOMULTICONVERSION_H

#include <numeric>
#include <boost/mpi.hpp>


#include <palabos/atomicBlock/dataProcessingFunctional2D.h>
#include <palabos/atomicBlock/dataProcessingFunctional3D.h>
#include <palabos/atomicBlock/dataField2D.h>

#include "palabosmpiutils.h"
#include "plbtypededuction.h"
#include "fieldvalueaccess.h"
#include "dimsizecalc.h"

namespace IPP
{

namespace SerialToMultiConversion
{

template<typename T>
void copyData(const std::vector<T>& arr, const plb::Box2D& domain, plb::ScalarField2D<T>& outputField)
{
    size_t srcIndex = 0;

    for(plb::plint x = domain.x0; x <= domain.x1; ++x)
    {
        for(plb::plint y = domain.y0; y <= domain.y1; ++y)
        {
            const T& value = arr[srcIndex];
            outputField.get(x, y) = value;
            ++srcIndex;
        }
    }

}

template<typename T>
void copyData(const std::vector<T>& arr, const plb::Box3D& domain, plb::ScalarField3D<T>& outputField)
{
    size_t srcIndex = 0;

    for(plb::plint x = domain.x0; x <= domain.x1; ++x)
    {
        for(plb::plint y = domain.y0; y <= domain.y1; ++y)
        {
            for(plb::plint z = domain.z0; z <= domain.z1; ++z)
            {
                const T& value = arr[srcIndex];
                outputField.get(x, y, z) = value;
                ++srcIndex;
            }
        }
    }
}

struct FallbackDataAccess
{
    template<typename T>
    static void copyToVector(const std::vector<T>& input,
                             const std::vector<plb::plint>&,
                             const size_t,
                             std::vector<T>& result)
    {
        result = input;
    }
};


template<size_t dim>
struct DataAccess : public FallbackDataAccess
{

};

template<>
class DataAccess<2>
{
public:

    template<typename T, typename FieldType>
    static void copyToVector(const FieldType& input,
                             const std::vector<plb::plint>& allGlobalIndices,
                             const size_t nProcesses,
                             std::vector<T>& result)
    {
        assert(result.empty());

        constexpr size_t nComp = 4;

        for(size_t iProc = 0; iProc < nProcesses; ++iProc)
        {
            const size_t startIndex = iProc * nComp;
            const plb::plint x0 = allGlobalIndices[startIndex];
            const plb::plint x1 = allGlobalIndices[startIndex + 1];
            const plb::plint y0 = allGlobalIndices[startIndex + 2];
            const plb::plint y1 = allGlobalIndices[startIndex + 3];

            for(plb::plint x = x0; x <= x1; ++x)
            {
                for(plb::plint y = y0; y <= y1; ++y)
                {
                    const T& val = FieldValueAccess<2>::get(input, x, y);
                    result.push_back(val);
                }
            }
        }
    }

};


template<>
class DataAccess<3>
{
public:

    template<typename T, typename FieldType>
    static void copyToVector(const FieldType& input,
                             const std::vector<plb::plint>& allGlobalIndices,
                             const size_t nProcesses,
                             std::vector<T>& result)
    {
        assert(result.empty());

        constexpr size_t nComp = 6;

        for(size_t iProc = 0; iProc < nProcesses; ++iProc)
        {
            const size_t startIndex = iProc * nComp;
            const plb::plint x0 = allGlobalIndices[startIndex];
            const plb::plint x1 = allGlobalIndices[startIndex + 1];
            const plb::plint y0 = allGlobalIndices[startIndex + 2];
            const plb::plint y1 = allGlobalIndices[startIndex + 3];
            const plb::plint z0 = allGlobalIndices[startIndex + 4];
            const plb::plint z1 = allGlobalIndices[startIndex + 5];


            for(plb::plint x = x0; x <= x1; ++x)
            {
                for(plb::plint y = y0; y <= y1; ++y)
                {
                    for(plb::plint z = z0; z <= z1; ++z)
                    {
                        const T& val = FieldValueAccess<3>::get(input, x, y, z);
                        result.push_back(val);
                    }
                }
            }
        }

    }
};


inline void collectIndices(const plb::Dot2D& offset,
                           const plb::Box2D& domain,
                           std::vector<plb::plint>& allGlobalIndices)
{
    allGlobalIndices.push_back(offset.x + domain.x0);
    allGlobalIndices.push_back(offset.x + domain.x1);
    allGlobalIndices.push_back(offset.y + domain.y0);
    allGlobalIndices.push_back(offset.y + domain.y1);
}

inline void collectIndices(const plb::Dot3D& offset,
                           const plb::Box3D& domain,
                           std::vector<plb::plint>& allGlobalIndices)
{
    allGlobalIndices.push_back(offset.x + domain.x0);
    allGlobalIndices.push_back(offset.x + domain.x1);
    allGlobalIndices.push_back(offset.y + domain.y0);
    allGlobalIndices.push_back(offset.y + domain.y1);
    allGlobalIndices.push_back(offset.z + domain.z0);
    allGlobalIndices.push_back(offset.z + domain.z1);
}

template<typename SizeCalc, typename T, size_t dim,
         typename ScalarField, typename Box, typename SerialDataField>
void scatterData(const SerialDataField& serialData, const Box& domain,
                 ScalarField& decomposedMultiData)
{
    const bool isMainProc = plb::global::mpi().isMainProcessor();
    const size_t nProcesses = plb::global::mpi().getSize();

    const PalabosMPIUtils::ProcInfo procInfo(isMainProc, nProcesses);


    std::vector<plb::plint> allGlobalIndices;
    PalabosMPIUtils::gatherAllIndices(procInfo, decomposedMultiData.getLocation(), domain, allGlobalIndices);


    std::vector<int> displs;
    std::vector<int> sndCounts;
    std::vector<T> globalArr;
    if(isMainProc)
    {
        PalabosMPIUtils::transferStrides<SizeCalc>(allGlobalIndices, nProcesses, sndCounts, displs);
        DataAccess<SizeCalc::dim>::copyToVector(serialData, allGlobalIndices, nProcesses, globalArr);
    }

    const size_t nElems = domain.nCells();
    std::vector<T> localArr(nElems);

    const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<T>();

    MPI_Scatterv(globalArr.data(), sndCounts.data(), displs.data(), mpiType,
                 localArr.data(), nElems, mpiType,
                 0, plb::global::mpi().getGlobalCommunicator());


    // copy to output
    copyData(localArr, domain, decomposedMultiData);


}

template<typename T, typename SerialField, size_t dim>
class ScalarFieldConvertT :
        public PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T, dim>::value
{
public:
    ScalarFieldConvertT(const SerialField& serialData)
        : m_serialData(serialData)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::staticVariables;
    }

    virtual ScalarFieldConvertT<T, SerialField, dim>* clone() const
    {
        return new ScalarFieldConvertT<T, SerialField, dim>(*this);
    }

    using BaseClass = typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T, dim>::value;
    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

    virtual void process(Box domain, InputField& decomposedMultiData)
    {
        using SizeCalc = DimSizeCalc<dim>;
        scatterData<SizeCalc, T, dim>(m_serialData, domain, decomposedMultiData);
    }

private:
    const SerialField& m_serialData;
};

}

} // end of namespace IPP

#endif // SERIALTOMULTICONVERSION_H
