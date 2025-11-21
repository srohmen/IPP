#ifndef MULTITOSERIALCONVERSION_H
#define MULTITOSERIALCONVERSION_H

#include <numeric>
#include <boost/mpi.hpp>

#include "palabosmpiutils.h"
#include "plbtypededuction.h"
#include "fieldvalueaccess.h"
#include "dimsizecalc.h"

namespace IPP
{

namespace MultiToSerialConversion
{

struct FallbackDataAccess
{
    template<typename T>
    static void copyData(const size_t,
                         const std::vector<plb::plint>&,
                         const std::vector<T>& globalArr,
                         std::vector<T>& outputField)
    {
        outputField = globalArr;
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
    static void copyData(const size_t nProcs, const std::vector<plb::plint>& allGlobalIndices,
                         const std::vector<T>& globalArr, FieldType& outputField)
    {
        constexpr size_t nComp = 4;

        size_t srcIndex = 0;
        for (size_t iProc = 0; iProc < nProcs; ++iProc)
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
                    const T& value = globalArr[srcIndex];
                    FieldValueAccess<2>::get(outputField, x, y) = value;
                    ++srcIndex;
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
    static void copyData(const size_t nProcs, const std::vector<plb::plint>& allGlobalIndices,
                         const std::vector<T>& globalArr, FieldType& outputField)
    {
        constexpr size_t nComp = 6;

        size_t srcIndex = 0;
        for (size_t iProc = 0; iProc < nProcs; ++iProc)
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
                        const T& value = globalArr[srcIndex];
                        FieldValueAccess<3>::get(outputField, x, y, z) = value;
                        ++srcIndex;
                    }
                }
            }
        }
    }


};

template<typename T>
void copyToVector(const plb::ScalarField2D<T>& input,
                  const plb::Box2D& domain,
                  std::vector<T>& result)
{
    assert(result.empty());

    for(plb::plint x = domain.x0; x <= domain.x1; ++x)
    {
        for(plb::plint y = domain.y0; y <= domain.y1; ++y)
        {
            const T& val = input.get(x, y);
            result.push_back(val);
        }
    }
}

template<typename T>
void copyToVector(const plb::ScalarField3D<T>& input,
                  const plb::Box3D& domain,
                  std::vector<T>& result)
{
    assert(result.empty());

    for(plb::plint x = domain.x0; x <= domain.x1; ++x)
    {
        for(plb::plint y = domain.y0; y <= domain.y1; ++y)
        {
            for(plb::plint z = domain.z0; z <= domain.z1; ++z)
            {
                result.push_back(input.get(x, y, z));
            }
        }
    }
}

template<typename SizeCalc, typename T,
         typename ScalarField, typename Box, typename SerialDataField>
static void gatherData(const ScalarField& decomposedMultiData, const Box& domain,
                       SerialDataField& serialData)
{
    const bool isMainProc = plb::global::mpi().isMainProcessor();
    const size_t nProcesses = plb::global::mpi().getSize();

    const PalabosMPIUtils::ProcInfo procInfo(isMainProc, nProcesses);

    std::vector<plb::plint> allGlobalIndices;
    PalabosMPIUtils::gatherAllIndices(procInfo, decomposedMultiData.getLocation(),
                                      domain, allGlobalIndices);

    std::vector<int> displs;
    std::vector<int> recvcounts;
    if(isMainProc)
    {
        PalabosMPIUtils::transferStrides<SizeCalc>(allGlobalIndices, nProcesses,
                                                   recvcounts, displs);
    }

    const int size = std::accumulate(recvcounts.begin(), recvcounts.end(), 0);
    std::vector<T> globalArr(size);

    const size_t nElems = domain.nCells();


    std::vector<T> inputArr;
    copyToVector(decomposedMultiData, domain, inputArr);

    const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<T>();

    MPI_Gatherv(inputArr.data(), nElems, mpiType,
                globalArr.data(), recvcounts.data(), displs.data(), mpiType,
                0, plb::global::mpi().getGlobalCommunicator());

    if(isMainProc)
    {
        // copy to output
        assert(allGlobalIndices.size() % SizeCalc::nComp == 0);
        DataAccess<SizeCalc::dim>::copyData(nProcesses, allGlobalIndices,
                                            globalArr, serialData);
    }

}

template<typename T, typename SerialField, size_t dim>
class ScalarFieldConvertT
        : public PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T, dim>::value
{
public:
    ScalarFieldConvertT(SerialField& result)
        : m_serialData(result)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
    }

    virtual ScalarFieldConvertT<T,SerialField,dim>* clone() const
    {
        return new ScalarFieldConvertT<T,SerialField,dim>(*this);
    }


    using BaseClass = typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T, dim>::value;
    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

    virtual void process(Box domain, InputField& decomposedMultiData)
    {
        using SizeCalc = DimSizeCalc<dim>;
        gatherData<SizeCalc, T>(decomposedMultiData, domain, m_serialData);
    }

private:
    SerialField& m_serialData;
};



}

} // end of namespace IPP

#endif // MULTITOSERIALCONVERSION_H
