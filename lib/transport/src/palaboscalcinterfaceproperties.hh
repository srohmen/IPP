#ifndef PALABOSCALCINTERFACEPROPERTIES_HH
#define PALABOSCALCINTERFACEPROPERTIES_HH

#include "palaboscalcinterfaceproperties.h"


#include "ippexception.h"

#include "concarrayview.h"
#include "mpitools.h"
#include "fielddecomposition.h"
#include "palaboslatticevalueaccess.h"
#include "palabosconversiontools.h"
#include "dotfindneighbors.h"
#include "permeabiltyflag.h"
#include "setupblockperiodicity.h"

#include "settontensorfunction.h"
#include "setntensorfromarray.h"


namespace IPP
{

template<typename T, typename MaskType, size_t dim>
class AveragePorosWeighted : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD<dim>::value
{
public:
    AveragePorosWeighted(const T& weight)
        : m_weight(weight)
    {

    }

    virtual AveragePorosWeighted<T, MaskType, dim>* clone() const override
    {
        return new AveragePorosWeighted<T, MaskType, dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing; // porosity
        modified[1] = plb::modif::nothing; // perm mask
        modified[2] = plb::modif::staticVariables; // avg poros
    }

private:
    using ScalarField = typename PlbTypeDeduction::GetScalarField_XD<T, dim>::value;
    using MaskField = typename PlbTypeDeduction::GetScalarField_XD<MaskType, dim>::value;

    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD<dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::processGenericBlocks)>;
    using DotList = typename Traits::template argument<1>::type;
    using InputFieldVec = typename Traits::template argument<2>::type;

public:
    virtual void processGenericBlocks(DotList const& dotList, InputFieldVec fields) override
    {
        const ScalarField& porosityField = *dynamic_cast<ScalarField*>(fields[0]);
        const MaskField& permMask = *dynamic_cast<MaskField*>(fields[1]);
        ScalarField& avgPoros = *dynamic_cast<ScalarField*>(fields[2]);

        const auto offsetPermMask = plb::computeRelativeDisplacement(porosityField, permMask);
        const auto offsetAvgPoros = plb::computeRelativeDisplacement(porosityField, avgPoros);

        for(const auto& center : dotList.dots)
        {
            // find neighbors
            // 2D -> n = 4
            // 3D -> n = 6
            const DotList neighList = DotFindNeighbors::find(center);
            assert(neighList.getN() == dim * 2);

            T porosSum = 0.0;
            for(int iNeigh = 0; iNeigh < neighList.getN(); ++iNeigh)
            {
                const auto& neigh = neighList.getDot(iNeigh);

                const MaskType& permVal = DataAccess::get(permMask, neigh + offsetPermMask);
                if(permVal & PF_isPermeable)
                {
                    const T& porosity = DataAccess::get(porosityField, neigh);
                    assert(porosity >= 0.0 && porosity <= 1.0);

                    const T weightPoros = porosity * m_weight;
                    assert(std::isnan(weightPoros) == false);

                    porosSum += weightPoros;
                    assert(std::isnan(porosSum) == false);
                    assert(porosSum >= 0.0 && porosSum <= 1.0);
                }
            }

            DataAccess::get(avgPoros, center + offsetAvgPoros) = porosSum;
        }
    }

private:
    const T m_weight;
};


template<typename T, typename MaskType, typename ArrayViewVec, size_t dim>
class AverageConcWeighted : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD<dim>::value
{
public:
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

    AverageConcWeighted(const T& weight, ArrayViewVec& avgConcField, const Dot& location)
        : m_weight(weight)
        , m_avgConcField(avgConcField)
        , m_location(location)
    {

    }

    virtual AverageConcWeighted<T, MaskType, ArrayViewVec, dim>* clone() const override
    {
        return new AverageConcWeighted<T, MaskType, ArrayViewVec, dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing; // porosField
        modified[1] = plb::modif::nothing; // avgPorosField
        modified[2] = plb::modif::nothing; // concField
        modified[3] = plb::modif::nothing; // permField
    }

private:
    using ScalarField = typename PlbTypeDeduction::GetScalarField_XD<T, dim>::value;
    using MaskField = typename PlbTypeDeduction::GetScalarField_XD<MaskType, dim>::value;

    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD<dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::processGenericBlocks)>;
    using DotList = typename Traits::template argument<1>::type;
    using InputFieldVec = typename Traits::template argument<2>::type;

public:
    virtual void processGenericBlocks(DotList const& dotList, InputFieldVec fields) override
    {
        using NTensorField = typename PlbTypeDeduction::GetNTensorField_XD<T, dim>::value;

        const ScalarField& porosField = *static_cast<ScalarField*>(fields[0]);
        const ScalarField& avgPorosField = *static_cast<ScalarField*>(fields[1]);
        const NTensorField& concField = *static_cast<NTensorField*>(fields[2]);
        const MaskField& permMask = *static_cast<MaskField*>(fields[3]);

        const Dot offsetAvgPoros = plb::computeRelativeDisplacement(porosField, avgPorosField);
        const Dot offsetConc = plb::computeRelativeDisplacement(porosField, concField);
        const Dot offsetPerm = plb::computeRelativeDisplacement(porosField, permMask);

        const Dot offset = porosField.getLocation();

        const size_t nComps = m_avgConcField.bounds()[0];

        for(const auto& center : dotList.dots)
        {
            // find neighbors
            // 2D -> n = 4
            // 3D -> n = 6
            const DotList neighList = DotFindNeighbors::find(center);
            assert(neighList.getN() == dim * 2);

            const T& avgPoros = DataAccess::get(avgPorosField, center + offsetAvgPoros);


            const Dot posOutput = center + offset - m_location;

            using U = typename ArrayViewVec::value_type;
            std::vector<U*> avgConcVec(nComps, nullptr);
            for(size_t iComp = 0; iComp < nComps; ++iComp)
            {
                auto avgConcView = m_avgConcField[iComp];
                auto& avgConc = DataAccess::get(avgConcView, posOutput);
                avgConc = 0.0;
                avgConcVec[iComp] = &avgConc;
            }

            for(int iNeigh = 0; iNeigh < neighList.getN(); ++iNeigh)
            {
                const Dot& neigh = neighList.getDot(iNeigh);

                const MaskType& permVal = DataAccess::get(permMask, neigh + offsetPerm);
                if(permVal & PF_isPermeable)
                {
                    const T& poros = DataAccess::get(porosField, neigh);

                    const T frac = poros * m_weight / avgPoros;
                    assert(std::isnan(frac) == false);
                    assert(std::isinf(frac) == false);

                    const T* concVec = DataAccess::get(concField, neigh + offsetConc);

                    for(size_t iComp = 0; iComp < nComps; ++iComp)
                    {
                        const T& conc = concVec[iComp];
                        const T weightedConc = conc * frac;
                        assert(std::isnan(weightedConc) == false);

                        U& avgConc = *(avgConcVec[iComp]);
                        avgConc += static_cast<U>(weightedConc);
                        assert(std::isnan(avgConc) == false);
                    }
                }
            }
        }
    }


private:
    const T m_weight;
    ArrayViewVec& m_avgConcField;
    const Dot& m_location;

};

template<typename T, typename MaskType, size_t dim>
PalabosCalcInterfaceProperties<T, MaskType, dim>::PalabosCalcInterfaceProperties(const T &weight,
                                                                                 const FieldDecomposition& decomp,
                                                                                 const DotList& interfaceNodes)
    : m_weight(weight)
    , m_decomp(decomp)
    , m_interfaceNodes(interfaceNodes)
    , m_permableMask(nullptr)
    , m_nComps(0)
{

}

template<typename T, typename MaskType, size_t dim>
void PalabosCalcInterfaceProperties<T, MaskType, dim>::initFields(const MaskField* permableMask,
                                                                  const Factory& factory)
{
    assert(m_permableMask == nullptr);
    assert(permableMask != nullptr);
    m_permableMask = permableMask;

    m_porosityField.reset(factory.template createMultiScalarFieldPtr<T>());
    m_avgPorosField.reset(factory.template createMultiScalarFieldPtr<T>());

    IPPCheck::assertCheck(m_nComps > 0);
    m_concField.reset(factory.template createMultiNTensorFieldPtr<T>(m_nComps));
}

template<typename T, typename MaskType, size_t dim>
void PalabosCalcInterfaceProperties<T, MaskType, dim>::setnComps(const size_t nComps)
{
    IPPCheck::assertCheck(m_nComps == 0);
    IPPCheck::assertCheck(nComps > 0);
    m_nComps = nComps;
}

template<typename T, typename MaskType, size_t dim>
void PalabosCalcInterfaceProperties<T, MaskType, dim>::execute(const std::vector<double>& conc,
                                                               const std::vector<double>& porosity,
                                                               std::vector<double>& avgConc,
                                                               std::vector<double>& avgPoros) const
{
    MPI_CHECK_SYNC;

    const IPPBox3DLong& ownDomain = m_decomp.getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    const IPPVector3DLong ownSize = ownDomain.getSize();
    const Dot dotOrigin = PalabosConversionTools::ToPlb<dim>::convertVector(location);

    this->averagePorosity(location, dotOrigin, ownSize, porosity, *m_avgPorosField, avgPoros);

    this->averageConc(location, dotOrigin, ownSize, *m_avgPorosField, conc, avgConc);

    MPI_CHECK_SYNC;
}

template<typename T, typename MaskType, size_t dim>
void PalabosCalcInterfaceProperties<T, MaskType, dim>::averagePorosity(const IPPVector3DLong& location,
                                                                       const Dot& dotOrigin,
                                                                       const IPPVector3DLong& ownSize,
                                                                       const std::vector<double>& porosity,
                                                                       ScalarField& avgPorosField,
                                                                       std::vector<double>& avgPoros) const
{
    using ArrViewConst = av::array_view<const double, dim>;

    const av::bounds<dim> bounds = av::MakeBounds<dim>::get(ownSize[0], ownSize[1], ownSize[2]);
    const ArrViewConst porosView(porosity, bounds);

    using SetFunc = SetScalarsFromArray<const double, dim>;
    plb::setToFunction(*m_porosityField, m_porosityField->getBoundingBox(),
                       SetFunc(porosView, location));

    std::vector<MultiBlock*> fields;
    fields.push_back(&(*m_porosityField));
    fields.push_back(const_cast<MaskField*>(m_permableMask));
    fields.push_back(&avgPorosField);

    using AveragePoros = AveragePorosWeighted<T, MaskType, dim>;
    plb::applyProcessingFunctional(new AveragePoros(m_weight), m_interfaceNodes, fields);

    assert(avgPoros.size() == porosity.size());
    using ArrView = av::array_view<double, dim>;
    ArrView avgPorosView(avgPoros, bounds);
    using GetFunc = GetScalarsToArray<T, decltype(avgPorosView), dim>;
    plb::applyProcessingFunctional(new GetFunc(avgPorosView, dotOrigin),
                                   avgPorosField.getBoundingBox(), avgPorosField);


    MPI_CHECK_SYNC;
}

template<typename T, typename MaskType, size_t dim>
void PalabosCalcInterfaceProperties<T, MaskType, dim>::averageConc(const IPPVector3DLong& location,
                                                                   const Dot& dotOrigin,
                                                                   const IPPVector3DLong& ownSize,
                                                                   const ScalarField& avgPorosField,
                                                                   const std::vector<double>& conc,
                                                                   std::vector<double>& avgConc) const
{

    assert(avgConc.size() == conc.size());
    auto concView = ConcArrayViewTool::makeConcView<dim>(m_nComps, ownSize, conc);
    using ConstConcArrayView = decltype(concView);


    std::vector<MultiBlock*> fields(4);
    fields[0] = &(*m_porosityField);
    fields[1] = &const_cast<ScalarField&>(avgPorosField);

    // conc field is only temporary object.
    // due to performance reasons memory is allocated in init process only once
    NTensorField& concField = const_cast<NTensorField&>(*m_concField);
    fields[2] = &concField;

    fields[3] = const_cast<MaskField*>(m_permableMask);



    using Setter = SetNTensorFromArray<T, ConstConcArrayView, dim>;

    using Sync = SetToNTensorFunction<T, Setter, dim>;
    plb::applyProcessingFunctional(new Sync(Setter(concView, location, m_nComps)),
                                   concField.getBoundingBox(), concField);

    // needs communication in between due non local operation in AverageConc
    auto avgConcView = ConcArrayViewTool::makeConcView<dim>(m_nComps, ownSize, avgConc);
    using ConcArrayView = decltype(avgConcView);

    using AverageConc = AverageConcWeighted<T, MaskType, ConcArrayView, dim>;
    plb::applyProcessingFunctional(new AverageConc(m_weight, avgConcView, dotOrigin),
                                   m_interfaceNodes, fields);

}


}


#endif // PALABOSCALCINTERFACEPROPERTIES_HH
