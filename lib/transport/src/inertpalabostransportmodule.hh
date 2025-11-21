#ifndef INERTPALABOSTRANSPORTMODULE_HH
#define INERTPALABOSTRANSPORTMODULE_HH

#include "inertpalabostransportmodule.h"

#include <boost/range/algorithm/set_algorithm.hpp>

#include "palabosboundaryconditions.h"
#include "simulationexchangedata.h"
#include "palaboslatticevalueaccess.h"
#include "configurebouncebackdynamics.h"
#include "dotsetexternalscalar.h"
#include "setdotstoscalar.h"
#include "fielddecomposition.h"
#include "permflagsutils.h"

#include "mpitools.h"

namespace IPP
{


template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::updateFieldPorosity()
{
    using ArrView = av::array_view<const double, dim>;

    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    const IPPVector3DLong size = ownDomain.getSize();
    const av::bounds<dim> bounds = av::MakeBounds<dim>::get(size[0], size[1], size[2]);

    ScalarField& physicalPorosity = this->m_palabosData.getPorosity();

    const Box bbox = physicalPorosity.getBoundingBox();

    const std::vector<double>& localPoros = this->m_simData.getPorosity();
    const ArrView porosView(localPoros, bounds);


    plb::setToFunction(physicalPorosity, bbox,
                       SetScalarsFromArray<const double, dim>(porosView, location));
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::initLatticeTransportProperties()
{
    using ArrView = av::array_view<const Scalar, dim>;

    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    const IPPVector3DLong size = ownDomain.getSize();
    const av::bounds<dim> bounds = av::MakeBounds<dim>::get(size[0], size[1], size[2]);

    ScalarField transScalar = this->m_plbFactory->template createMultiScalarField<Scalar>();

    const Box bbox = transScalar.getBoundingBox();

    const ArrView transScalarView(this->m_transportScalar, bounds);
    plb::setToFunction(transScalar, bbox,
                       SetScalarsFromArray<const Scalar, dim>(transScalarView, location));



    ScalarField& porosity = this->m_palabosData.getPorosity();

    for(DiffusionLatticePtr& lattice : this->m_palabosData.diffLattices)
    {
        using DiffusionDesc = typename TransportTraits::template DiffusionDescriptorT<Scalar>;

        plb::setExternalScalar(*lattice, bbox,
                               DiffusionDesc::ExternalField::porosityBeginsAt,
                               porosity);

        plb::setExternalScalar(*lattice, bbox,
                               DiffusionDesc::ExternalField::transportScalarBeginsAt,
                               transScalar);
    }
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::initModule()
{
    MPI_CHECK_SYNC;

    this->initInertCellsMarker();


    const AbstractBoundaryConditions& bc = *this->m_simData.getBoundaryConditions();
    const std::vector<int>& periodic = bc.getPeriodicBoundaryDimensions();

    if(this->m_conf->flowFunc->isEnabled())
    {
        const AbstractBoundaryConditions::BCVec& advectiveBC = bc.getAdvectiveBoundaryConditions();
        this->initHydrodynamicLattice(periodic, advectiveBC);
        this->prepareResultArrays(this->m_results.flowResults);
    }

    MPI_CHECK_SYNC;

    this->initDiffusionLattices(bc);
    MPI_CHECK_SYNC;


    this->updateFieldPorosity();
    this->updateFieldTransportScalar();
    this->initLatticeTransportProperties();
    MPI_CHECK_SYNC;

    // this->updateDistanceField();
    MPI_CHECK_SYNC;


    this->initLatticeConcentrations();
    MPI_CHECK_SYNC;

    this->finishLatticeInit();
    MPI_CHECK_SYNC;


    // this->initPostTransportConc();
    MPI_CHECK_SYNC;


    this->initBounceBackNodes();
    MPI_CHECK_SYNC;
}

template<typename T>
static size_t getSize(const plb::MultiScalarField2D<T>& field, const size_t dim)
{
    if(dim == 0)
    {
        return field.getNx();
    }
    else if(dim == 1)
    {
        return field.getNy();
    }
    else
    {
        throw std::runtime_error("impossible dim for 2D field: " + std::to_string(dim));
    }
}

template<typename T>
static size_t getSize(const plb::MultiScalarField3D<T>& field, const size_t dim)
{
    if(dim == 0)
    {
        return field.getNx();
    }
    else if(dim == 1)
    {
        return field.getNy();
    }
    else if(dim == 2)
    {
        return field.getNz();
    }
    else
    {
        throw std::runtime_error("impossible dim for 3D field: " + std::to_string(dim));
    }
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::initLatticeConcentrations()
{
    const SimulationExchangeData::DimToDiffSpecies& dimToDiff = this->m_simData.getDimToDiffSpecies();
    assert(dimToDiff.size() == this->m_palabosData.diffLattices.size());

    for(auto it = dimToDiff.begin(); it != dimToDiff.end(); ++it)
    {
        const size_t diffDim = it->first;
        const size_t iComp = it->second;

        DiffusionLatticePtr diffLattice = this->m_palabosData.diffLattices.at(iComp);
        ScalarField concTmp = this->m_plbFactory->template createMultiScalarField(0.0);

        // set to linear gradient
        plb::setToCoordinate(concTmp, concTmp.getBoundingBox(), diffDim);

        const size_t size = getSize(concTmp, diffDim);
        const Scalar gradVal = (Scalar)1.0 / (Scalar)(size - 1);
        plb::multiplyInPlace(concTmp, gradVal);
        ScalarFieldPtr conc = plb::subtract((Scalar)1.0, concTmp);

        // plb::divideInPlace(*conc, this->m_palabosData.getPorosity());

        using InitFunc = typename TransportTraits::PorosityFunc;
        InitFunc::apply(*diffLattice, *conc);
    }
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::setGeomSync(AbstractGeometrySync *geomSync)
{
    // do nothing
}


template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::updatePostReactionState()
{
    // do nothing
}


template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::addPostStreamOperation(DiffusionLattice &lattice)
{
    // do nothing??
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::collectData()
{
    // do nothing
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::initBounceBackNodes()
{
    MPI_CHECK_SYNC;

    DotList permNodes, nonPermNodes;
    this->collectPermNodes(permNodes, nonPermNodes);

    using MaskType = typename TransportData::MaskType;
    using ScalarFieldMask = typename TransportData::ScalarFieldMask;
    ScalarFieldMask& permMask = *this->m_palabosData.permMask;

    using Enable = BoxSetPermMaskFlag<MaskType, dim, EnableFlag>;
    plb::applyProcessingFunctional(new Enable(PF_isPermeable), permMask.getBoundingBox(), permMask);

    using Disable = DotSetPermMaskFlag<MaskType, dim, DisableFlag>;
    plb::applyProcessingFunctional(new Disable(PF_isPermeable), nonPermNodes, permMask);


    auto& diffLattices = this->m_palabosData.diffLattices;

    const SimulationExchangeData::DimToDiffSpecies& dimToDiff =
            this->m_simData.getDimToDiffSpecies();
    IPPCheck::assertCheck(dimToDiff.size() == diffLattices.size());

    // prepare bounce back mask
    using SetValue = SetDotsToScalar_S<Scalar, dim>;
    ScalarField maskTmp = this->m_plbFactory->createMultiScalarField(0);
    const int maskVal = 1;
    plb::applyProcessingFunctional(new SetValue(maskVal), nonPermNodes, maskTmp);
    const auto bbMask = plb::copyConvert<Scalar, int>(maskTmp);

    for(size_t iComp = 0; iComp < diffLattices.size(); ++iComp)
    {
        DiffusionLattice& lattice = *diffLattices[iComp];
        plb::defineDynamics(lattice, *bbMask,
                            new plb::BounceBack<Scalar, TransportTraits::template DiffusionDescriptorT>,
                            maskVal);
    }
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::collectPermNodes(DotList& permNodes,
                                                                    DotList& nonPermNodes)
{
    MPI_CHECK_SYNC;

    ScalarField& porosity = this->m_palabosData.getPorosity();

    using CollectSolidCoords = CollectThreshCoordinates<Scalar, dim>;
    DotList lowPorosNodes;
    plb::applyProcessingFunctional(new CollectSolidCoords(this->m_conf->porosLow, permNodes, lowPorosNodes),
                                   porosity.getBoundingBox(), porosity);

    std::sort(permNodes.dots.begin(), permNodes.dots.end());
    nonPermNodes.dots.erase(std::unique(permNodes.dots.begin(), permNodes.dots.end()),
                            permNodes.dots.end()
                            );


    // appending static inert cells
    const DotList& inertNodes = this->m_palabosData.getInertSolidCells();
    boost::range::set_union(lowPorosNodes.dots, inertNodes.dots, std::back_inserter(nonPermNodes.dots));

    std::sort(nonPermNodes.dots.begin(), nonPermNodes.dots.end());
    nonPermNodes.dots.erase(std::unique(nonPermNodes.dots.begin(), nonPermNodes.dots.end()),
                            nonPermNodes.dots.end()
                            );

    MPI_CHECK_SYNC;
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::saveCheckpoint() const
{
    this->saveLatticeCheckpoint();
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::loadCheckpoint(const size_t iteration)
{
    this->loadLatticeCheckpoint(iteration);
}

template<typename TransportTraits>
void InertPalabosTransportModule<TransportTraits>::prepareResults()
{
    // we have to update the post RM conc when no chemistry is enabled

    auto& diffLattices = this->m_palabosData.diffLattices;
    for(size_t iComp = 0; iComp < diffLattices.size(); ++iComp)
    {
        DiffusionLattice& diffLattice = *diffLattices[iComp];

        ScalarField& conc = this->m_palabosData.postReactionConcentrations[iComp];
        plb::computeDensity(diffLattice, conc, diffLattice.getBoundingBox());
    }
}

}

#endif // INERTPALABOSTRANSPORTMODULE_HH
