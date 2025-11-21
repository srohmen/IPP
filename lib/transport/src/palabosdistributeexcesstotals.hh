#ifndef PALABOSDISTRIBUTEEXCESSTOTALS_HH
#define PALABOSDISTRIBUTEEXCESSTOTALS_HH

#include "palabosdistributeexcesstotals.h"

#include "fielddecomposition.h"
#include "palaboslatticevalueaccess.h"
#include "palabosconversiontools.h"
#include "arraydimensionconvert.h"
#include "concarrayview.h"
#include "mpitools.h"
#include "setdotstoscalar.h"
#include "dotfindneighbors.h"

namespace IPP
{

template<typename T, size_t dim>
class DistributeSourceToNeighbors : public PlbTypeDeduction::
        GetScalarFieldDotProcessingFunctionalXD<T,dim>::value
{
public:
    DistributeSourceToNeighbors()
    {

    }

    virtual DistributeSourceToNeighbors<T,dim>* clone() const override
    {
        return new DistributeSourceToNeighbors<T,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing; // porosity
        modified[1] = plb::modif::nothing; // local source
        modified[2] = plb::modif::staticVariables; // dist source
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetScalarFieldDotProcessingFunctionalXD<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotListInput = typename Traits::template argument<1>::type;
    using InputFieldVec = typename Traits::template argument<2>::type;


public:
    virtual void process(DotListInput const& dotList, InputFieldVec fields) override
    {
        const auto& porosityField = *fields[0];
        const auto& sourceField = *fields[1];
        auto& distSourceField = *fields[2];

        const auto offsetSource = plb::computeRelativeDisplacement(porosityField, sourceField);
        const auto offsetDist = plb::computeRelativeDisplacement(porosityField, distSourceField);

        for(const auto& pos : dotList.dots)
        {
            // find neighbors
            // 2D -> n = 4
            // 3D -> n = 6
            auto neighList = DotFindNeighbors::find(pos);
            assert(neighList.getN() == dim * 2);


            // eliminate low poros cells
            const auto removeFunc = [&](const auto& pos)
            {
                const T& porosity = DataAccess::get(porosityField, pos);
                if(porosity < 0.001)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            };
            neighList.dots.erase(std::remove_if(neighList.dots.begin(),
                                                neighList.dots.end(), removeFunc),
                                 neighList.dots.end());


            T porosSum = 0.0;
            for(int iNeigh = 0; iNeigh < neighList.getN(); ++iNeigh)
            {
                const auto& neigh = neighList.getDot(iNeigh);
                const T& porosity = DataAccess::get(porosityField, neigh);
                porosSum += porosity;
            }

            const T& source = DataAccess::get(sourceField, pos + offsetSource);

            // assign weighting with respect to porosity
            for(int iNeigh = 0; iNeigh < neighList.getN(); ++iNeigh)
            {
                const auto& neigh = neighList.getDot(iNeigh);
                const T& porosity = DataAccess::get(porosityField, neigh);
                const T weight = porosity / porosSum;
                const T neighSource = source * weight;

                DataAccess::get(distSourceField, neigh + offsetDist) = neighSource;
            }


        }
    }

private:

};

template<typename T, typename Field>
PalabosDistributeExcessTotals<T, Field>::PalabosDistributeExcessTotals(const FieldDecomposition& decomp,
                                                                       const Factory& factory)
    : m_decomp(decomp)
    , m_factory(factory)
    , m_nComps(0)
{
    MPI_CHECK_SYNC;
}

template<typename T, typename Field>
void PalabosDistributeExcessTotals<T, Field>::init()
{
    m_porosity.reset(m_factory.template createMultiScalarFieldPtr<T>());
}

template<typename T, typename Field>
void PalabosDistributeExcessTotals<T, Field>::setnComps(const size_t nComps)
{
    IPPCheck::assertCheck(m_nComps == 0);
    IPPCheck::assertCheck(nComps > 0);
    m_nComps = nComps;
}

template<typename T, typename Field>
void PalabosDistributeExcessTotals<T, Field>::updatePorosity(const std::vector<double>& porosity)
{
    MPI_CHECK_SYNC;

    using ArrView = av::array_view<const double, dim>;

    const IPPBox3DLong& ownDomain = m_decomp.getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    const IPPVector3DLong size = ownDomain.getSize();
    const av::bounds<dim> bounds = av::MakeBounds<dim>::get(size[0], size[1], size[2]);

    const ArrView porosView(porosity, bounds);

    using SetFunc = SetScalarsFromArray<const double, dim>;
    plb::setToFunction(*m_porosity, m_porosity->getBoundingBox(),
                       SetFunc(porosView, location));

    MPI_CHECK_SYNC;
}

template<typename T, typename Field>
void PalabosDistributeExcessTotals<T, Field>::execute(const std::vector<CellTotalsDiff>& diffPerCell,
                                                      std::vector<double>& distributedSource)
{

    MPI_CHECK_SYNC;

    const IPPBox3DLong& ownDomain = m_decomp.getOwnDomain();
    const IPPVector3DLong& origin = ownDomain.lower;
    const IPPVector3DLong ownSize = ownDomain.getSize();

    distributedSource.resize(ownDomain.getDiagonalSize() * m_nComps);
    auto sourceView = ConcArrayViewTool::makeConcView<dim>(m_nComps, ownSize, distributedSource);

    const auto dotOrigin = PalabosConversionTools::ToPlb<dim>::convertVector(origin);


    using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;
    DotList dotList;

    const ArrayDimensionConvert localConv(ownSize);
    for(const CellTotalsDiff& data : diffPerCell)
    {
        const IPPVector3DInt localCoordInt = localConv.calcCoordinate(data.iCell);
        const IPPVector3DLong localCoord = convert<long>(localCoordInt);
        const IPPVector3DLong globalCoord = localCoord + origin;

        const auto dot = PalabosConversionTools::ToPlb<dim>::convertVector(globalCoord);
        dotList.addDot(dot);
    }



    for(size_t iComp = 0; iComp < m_nComps; ++iComp)
    {
        Field sourceField = m_factory.createMultiScalarField(0.0);
        Field distSourceField = m_factory.createMultiScalarField(0.0);

        MPI_CHECK_SYNC;

        DotList pos;
        T source = -1.0E+9;

        for(const CellTotalsDiff& data : diffPerCell)
        {
            const IPPVector3DInt localCoordInt = localConv.calcCoordinate(data.iCell);
            const IPPVector3DLong localCoord = convert<long>(localCoordInt);
            const IPPVector3DLong globalCoord = localCoord + origin;

            const auto dot = PalabosConversionTools::ToPlb<dim>::convertVector(globalCoord);

            source = data.totalsDiff[iComp];
            pos.addDot(dot);
        }


        using SetVal = SetDotsToScalar_S<T,dim>;
        plb::applyProcessingFunctional(new SetVal(source), pos, sourceField);

        MPI_CHECK_SYNC;

        std::vector<Field*> fields;
        fields.push_back(&(*m_porosity));
        fields.push_back(&sourceField);
        fields.push_back(&distSourceField);

        using Distribute = DistributeSourceToNeighbors<T,dim>;
        plb::applyProcessingFunctional(new Distribute, dotList, fields);

        MPI_CHECK_SYNC;


        auto sourceViewSlice = sourceView[iComp];

        using Sync = GetScalarsToArray<T, decltype(sourceViewSlice), dim>;
        plb::applyProcessingFunctional(new Sync(sourceViewSlice, dotOrigin),
                                       distSourceField.getBoundingBox(), distSourceField);

    }

    MPI_CHECK_SYNC;
}


}

#endif // PALABOSDISTRIBUTEEXCESSTOTALS_HH
