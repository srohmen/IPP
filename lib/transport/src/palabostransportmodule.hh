#ifndef PALABOSTRANSPORTMODULE_HH
#define PALABOSTRANSPORTMODULE_HH

#include "palabostransportmodule.h"

#include <chrono>

#include <boost/optional.hpp>
#include <boost/filesystem/operations.hpp>


#include "abstractboundaryconditions.h"
#include "scopedfloatingpointexception.h"

#include "concsetinit.h"

#include "palabosutils.h"
#include "palabosio.h"

#include "simulationexchangedata.h"

#include "countgeomchanged.h"
#include "ippconstants.h"
#include "ippvector.h"


#include "ippconfig.h"
#include "abstractflowfunctor.h"
#include "abstractboundaryconditions.h"

#include "palabostransportdata.h"
#include "palabosresultwriter.h"

#include "ipprenderer.h"
#include "virtualpalabosfieldwrapper.h"

#include "palaboslatticevalueaccess.h"
#include "transportmoduleconfig.h"
#include "palabosobjectfactory.h"


#include "iszerofunc.h"
#include "distancefieldcalc.h"
#include "multitoserialconversion.h"
#include "serialtomulticonversion.h"
#include "filesystemutils.h"
#include "configurebouncebackdynamics.h"
#include "latticeboltzmannunitconversion.h"
#include "bench_tools.h"

#include "scalartransformations.h"


#include "hydrodynamicboundarycondition.h"
#include "advectiondiffusionboundarycondition.h"

#include "concarrayview.h"
#include "palabosboundaryconditions.h"
#include "setdotstoscalar.h"
#include "palabosio.h"

#include "syncperiodicdomain.h"
#include "getlocalboundingbox.h"
#include "palabosfielddecomposition.h"

#include "dotsetexternalscalar.h"
#include "setpermmaskflag.h"
#include "collectpermmaskflagcoordinates.h"

#include "setupblockperiodicity.h"

#include "createfilename.h"

// debug stuff
#include "printbox.h"
#include "printpopulation.h"
#include "printporosity.h"
#include "printdynamics.h"
#include "plbcheckdata.h"


#define HYDRODYNAMIC_DESCRIPTOR TransportTraits::template HydrodynamicDescriptorT
#define DIFFUSION_DESCRIPTOR TransportTraits::template DiffusionDescriptorT

namespace IPP
{


static void initDomainInfos(const std::vector<size_t>& forcedDims,
                            const std::vector<BoundaryConditionData*>& bcs,
                            DomainInfos& domainInfos, size_t& nx, size_t& ny, size_t& nz)
{
    domainInfos.origin = {{0, 0, 0}};

    namespace Setup = PalabosBoundaryConditions::BoundaryConditionSetup;
    domainInfos.boundaryEnvelops = Setup::calcNeededEnvelops(forcedDims, bcs);


    for(BoundaryConditionDomain::BoundaryPosition pos : domainInfos.boundaryEnvelops)
    {
        switch(pos)
        {
        case BoundaryConditionDomain::BP_left:
            --domainInfos.origin[0];
        case BoundaryConditionDomain::BP_right:
            ++nx;
            break;

        case BoundaryConditionDomain::BP_bottom:
            --domainInfos.origin[1];
        case BoundaryConditionDomain::BP_top:
            ++ny;
            break;

        case BoundaryConditionDomain::BP_back:
            --domainInfos.origin[2];
        case BoundaryConditionDomain::BP_front:
            ++nz;
            break;

        default:
        {
            throw std::runtime_error("unknown bc position");
            break;
        }

        }
    }


    domainInfos.totalSize[0] = nx;
    domainInfos.totalSize[1] = ny;
    domainInfos.totalSize[2] = nz;
}


static const std::string s_flowFile = "flow" + IPPConstants::s_dumpExtension;


static void plbMpiInit(const boost::mpi::communicator& comm)
{
    plb::global::mpi().init(comm);
    plb::global::plbRandom<float>().seed(10);
    plb::global::plbRandom<double>().seed(10);
    plb::global::plbRandom<plb::plint>().seed(10);
}

template<typename TransportTraits>
PalabosTransportModule<TransportTraits>::PalabosTransportModule(TransportModuleConfigPtr& config)
    : m_simData(config->simData)
    , m_conf(config)
    , nx(0), ny(0), nz(0)
    , m_transportSyncEnabledBox(PalabosConversionTools::ToPlb<dim>::convertBox(config->ippConf.transportSyncEnabledBox))
    , m_decomp(nullptr)
    , m_offset(PalabosConversionTools::toPlb(config->offset))
    , m_omegaCalc(*config)
    , m_convergenceCheck(config->convergenceTolerance, config->isConverged)
{
    // the init MUST be done before any multi MPI palabos object is constructed
    // otherwise dead lock / undefined behaviour
    plbMpiInit(config->globalComm);
    plb::global::directories().setOutputDir(config->resultDir.string() + "/");

    this->initSize();
}

template<typename TransportTraits>
PalabosTransportModule<TransportTraits>::~PalabosTransportModule()
{

}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::initSize()
{
    MPI_CHECK_SYNC;

    nx = m_conf->nx;
    ny = m_conf->ny;
    nz = m_conf->nz;

    initDomainInfos(m_conf->effDiffDims, m_conf->bcData,
                    m_domainInfos, nx, ny, nz);



    {
        const int tnx = m_domainInfos.totalSize[0];
        const int tny = m_domainInfos.totalSize[1];
        const int tnz = m_domainInfos.totalSize[2];

        IPPBox3DInt wholeDomain({0, 0, 0}, {tnx, tny, tnz});
        for(const IPPBox3DInt& box : m_conf->ippConf.decomp)
        {
            if(wholeDomain.contains(box) == false)
            {
                std::stringstream ss;
                ss << "decomp box not contained in domain: " << box;
                throw std::runtime_error(ss.str());
            }
        }
    }


    m_plbFactory.reset(new PlbFactory(nx, ny, nz, m_conf->ippConf.decomp, m_periodicity));

    // use arbitary dummy multi data field which has the same size
    // and decomposition as all other scalar fields which will constructed later
    ScalarField dummy = m_plbFactory->template createMultiScalarField<Scalar>();
    m_decomp.reset(new PalabosFieldDecomposition<Scalar, ScalarField>(dummy));


    const Box bbox = dummy.getBoundingBox();
    plb::intersect(m_transportSyncEnabledBox, bbox, m_transportSyncEnabledBox);


    const size_t nCells = m_decomp->getLocalNumberOfCells();
    m_transportScalar.resize(nCells, -1.0);
    m_simData.getDistField().resize(nCells, 0.0);

    MPI_CHECK_SYNC;
}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::init()
{    
    const AbstractBoundaryConditions& bc = *this->m_simData.getBoundaryConditions();
    m_periodicity = bc.getPeriodicBoundaryDimensions();


    // TODO: move construction to palabosData
    // init with non permeable value. envelops for non periodic boundarys must be initialized correctly
    // PF_isPermeable == false -> 0
    using MaskType = typename TransportData::MaskType;
    m_palabosData.permMask.reset(m_plbFactory->template createMultiScalarFieldPtr<MaskType>(0));


    const size_t nComp = m_simData.getCompNames().size();
    m_palabosData.init(*m_plbFactory, nComp);

    this->initModule();

    plb::pcout << plb::getMultiBlockInfo(*m_palabosData.diffLattices.front()) << std::endl;

    m_convergenceCheck.init(m_simData.getDimToDiffSpecies().size());
}

template<typename TransportTraits>
NonLocalOperations& PalabosTransportModule<TransportTraits>::getNonLocalOperations()
{
    return m_nonLocalOperations;
}

template<typename TransportTraits>
const DomainInfos& PalabosTransportModule<TransportTraits>::getDomainInfos() const
{
    return m_domainInfos;
}

template<typename TransportTraits>
template<typename Transformation>
void PalabosTransportModule<TransportTraits>::setupDiffBoundaryConditions(
        const AbstractBoundaryConditions& bc,
        const std::string& speciesName,
        const Transformation &trans,
        DiffusionLattice& lattice)
{
    const std::vector<int>& periodic = bc.getPeriodicBoundaryDimensions();
    SetupBlockPeriodicity::definePeriodicity(periodic, lattice);

    if(bc.contains(speciesName))
    {
        namespace PlbBC = PalabosBoundaryConditions;
        namespace Setup = PlbBC::BoundaryConditionSetup;
        using BCSetup = PlbBC::AdvectionDiffusionBoundaryCondition<dim>;

        const AbstractBoundaryConditions::BCVec& diffBCDef = bc.getDiffusiveBoundaryConditions(speciesName);
        Setup::defineBoundaryConditions<BCSetup>(lattice, periodic, diffBCDef, trans);
    }
}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::initHydrodynamicLattice(
        const std::vector<int>& periodic,
        const AbstractBoundaryConditions::BCVec& advectiveBC)
{
    HydrodynamicDynamic* flowBGDyn = new HydrodynamicDynamic(1.0);
    HydrodynamicLatticePtr flowLattice = m_plbFactory->template
            createMultiLattice<Scalar, HYDRODYNAMIC_DESCRIPTOR, HydrodynamicLattice>(flowBGDyn);
    m_palabosData.flowLattice = flowLattice;
    HydrodynamicLattice& lattice = *flowLattice;


    // 1) setup BC

    SetupBlockPeriodicity::definePeriodicity(periodic, lattice);

    using namespace PalabosBoundaryConditions;
    using SetupBC = HydrodynamicBoundaryCondition<dim>;

    namespace Setup = BoundaryConditionSetup;
    Setup::defineBoundaryConditions<SetupBC>(lattice, periodic, advectiveBC,
                                             ScalarTransformations::NoTransform());

    // 2) define initial conditions
    const AbstractFlowFunctor& flowFunc = *m_conf->flowFunc;
    plb::initializeAtEquilibrium(lattice, lattice.getBoundingBox(),
                                 typename PalabosInitFunctors<TransportTraits>::
                                 VelocityAndDensityWrapper(flowFunc));

    // 3) define transport scalar
    using Descriptor = typename HYDRODYNAMIC_DESCRIPTOR<Scalar>;
    PalabosLatticeValueAccess::setTransportScalarIfPossible<Descriptor>(lattice, (Scalar)0.0);


    flowLattice->toggleInternalStatistics(false);
}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::initDiffusionLattices(const AbstractBoundaryConditions& bc)
{
    const std::vector<std::string> &compNames = m_simData.getCompNames();
    const size_t nSpecies = compNames.size();

    typename ResultValCalc::DensityVelocityResultsVec& diffResultsVec = m_results.diffResultsVec;
    diffResultsVec.resize(nSpecies);

    for(size_t iComp = 0; iComp < nSpecies; ++iComp)
    {
        const std::string& speciesName = compNames[iComp];


        const Scalar initVal_omega_cPhi = m_omegaCalc.getInitDynamicsValue();
        DiffusionDynamic* baseDyn = new DiffusionDynamic(initVal_omega_cPhi);
        DiffusionLatticePtr diffLattice = m_plbFactory->template
                createMultiLattice<Scalar, DIFFUSION_DESCRIPTOR, DiffusionLattice>(baseDyn);
        m_palabosData.diffLattices.push_back(diffLattice);
        m_palabosData.componentNames.push_back(speciesName);


        // setup BCs and using dummy transformations
        //        const typename TransformationFactory::TransFwdBackPtr& trans =
        //                m_palabosData.scalarTransformations[iSpecies];
        ScalarTransformations::NoTransformationForwardBackward trans;
        setupDiffBoundaryConditions(bc, speciesName, trans.forward, *diffLattice);


        // setup the transport scalar (diffusion velocity or value to update TRT tuning value)
        using Descriptor = typename DIFFUSION_DESCRIPTOR<Scalar>;
        const Scalar val = m_omegaCalc.template calcTransportScalar<Descriptor>(initVal_omega_cPhi, 1.0);
        PalabosLatticeValueAccess::setTransportScalarIfPossible<Descriptor>(*diffLattice, val);

        PalabosLatticeValueAccess::setPorosityIfPossible<Descriptor>(*diffLattice, val);


        // prepare result arrays
        typename ResultValCalc::DensityVelocityResults& diffResults = diffResultsVec[iComp];
        prepareResultArrays(diffResults);

        // disable internal statistics due to performance reasons
        diffLattice->toggleInternalStatistics(false);


        this->addPostStreamOperation(*diffLattice);
    }
}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::prepareResultArrays(
        typename ResultValCalc::DensityVelocityResults& results)
{
    results.density.reset(m_plbFactory->template createMultiScalarFieldPtr<Scalar>());
    results.velocity.reset(m_plbFactory->template createMultiTensorFieldPtr<Scalar, dim>());
}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::finishLatticeInit()
{
    if(m_palabosData.flowLattice)
    {
        m_palabosData.flowLattice->initialize();
    }

    for(DiffusionLatticePtr& diffLattice : m_palabosData.diffLattices)
    {
        diffLattice->initialize();
    }
}


template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::initInertCellsMarker()
{
    typedef av::array_view<const char, dim> ArrView;

    const IPPBox3DLong& ownDomain = m_decomp->getOwnDomain();
    const Dot& location = PalabosConversionTools::ToPlb<dim>::convertVector(ownDomain.lower);
    const IPPVector3DLong size = ownDomain.getSize();

    const std::vector<char>& inertCells = m_simData.getInertSolidCells();
    const av::bounds<dim> bounds = av::MakeBounds<dim>::get(size[0], size[1], size[2]);
    const ArrView inertCellsView(inertCells, bounds);


    DotList inertCellsTmp;
    for(auto idx : inertCellsView.bounds())
    {
        if(inertCellsView[idx] == 1)
        {
            const Dot dot = PalabosConversionTools::makeDot(idx);
            const Dot absLocation = dot + location;
            inertCellsTmp.addDot(absLocation);
        }
    }

    using MaskType = typename TransportData::MaskType;
    using SetInert = DotSetPermMaskFlag<MaskType, dim, EnableFlag>;
    plb::applyProcessingFunctional(new SetInert(PF_isInertSolid), inertCellsTmp,
                                   *m_palabosData.permMask);


    DotList& inertCellsDots = m_palabosData.getInertSolidCells();
    using Collect = CollectPermMaskFlagCoordinates<MaskType, dim>;
    plb::applyProcessingFunctional(new Collect(PF_isInertSolid, inertCellsDots),
                                   m_palabosData.permMask->getBoundingBox(),
                                   *m_palabosData.permMask);


    std::sort(inertCellsDots.dots.begin(), inertCellsDots.dots.end());
    inertCellsDots.dots.erase(std::unique(inertCellsDots.dots.begin(),
                                          inertCellsDots.dots.end()),
                              inertCellsDots.dots.end());


}


template<typename TransportTraits>
const FieldDecomposition* PalabosTransportModule<TransportTraits>::getDecomposition() const
{
    assert(m_decomp);
    return m_decomp.get();
}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::updateFieldTransportScalar()
{
    using Descriptor = typename DIFFUSION_DESCRIPTOR<Scalar>;
    const auto transScalarCalc = m_omegaCalc.template makeTransportScalarCalcFunc<Descriptor>();

    const std::vector<double>& diffCoefs = m_simData.getDiffusionCoefs();
    for(size_t iCell = 0; iCell < diffCoefs.size(); ++iCell)
    {
        const Scalar relativeDiffCoef = diffCoefs[iCell] / m_conf->diffusionCoefRef;
        const Scalar transScalar = transScalarCalc(relativeDiffCoef);
        m_transportScalar[iCell] = transScalar;
    }



    // TODO: calculate transport scalars from permeability and update hydrodynamic lattice accordingly
    // update hydrodynamic lattice
    // HydrodynamicLatticePtr flowLattice = m_palabosData.flowLattice;
    // ConfigureBounceBackDynamics<TransportTraits>::apply(m_oldNonPermNodesHydr,
    // nonPermCells, *concField, *flowLattice);
    // PlbVA::setTransportScalarIfPossible<HydrodynamicDescriptor>(*flowLattice, *relativeDiff);
}


template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::run()
{
    const ScopedDisableFloatingPointException fpe;

    const size_t nComps = m_palabosData.diffLattices.size();

    HydrodynamicLatticePtr flowLattice = m_palabosData.flowLattice;
    if(flowLattice)
    {
        // update flow field
        // FIXME: update flow lattice until velocity convergence!
        throw std::runtime_error("FIXME");
        flowLattice->collideAndStream();

        const Box totalDomain = flowLattice->getBoundingBox();
        for(size_t iComp = 0; iComp < nComps; ++iComp)
        {
            const DiffusionLatticePtr& diffLattice = m_palabosData.diffLattices[iComp];
            plb::latticeToPassiveAdvDiff(*flowLattice, *diffLattice, totalDomain);
        }
    }

    // transport all concentrations
    for(size_t iComp = 0; iComp < nComps; ++iComp)
    {
        const DiffusionLatticePtr& diffLattice = m_palabosData.diffLattices[iComp];
        diffLattice->collideAndStream();
    }

}


template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::writeDebugData(const ResultsToWrite& resultsToWrite) const
{
    const Scalar dx = m_conf->spatialResolution;
    const std::array<size_t, 3> dims = {{nx, ny, nz }};

    plb::Array<Scalar,3> bcOffset;
    PalabosConversionTools::convertVector(m_domainInfos.origin, bcOffset);
    bcOffset *= dx;
    const plb::Array<Scalar,3> offset = m_offset + bcOffset;

    using Writer = ResultWriter::PalabosResultWriter<TransportTraits>;
    Writer::writeDebugData(resultsToWrite, dx, offset, dims, m_simData, m_palabosData, *m_plbFactory);
}



template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::writeResults(const ResultsToWrite& resultsToWrite)
{
    this->prepareResults();

    const Scalar dx = m_conf->spatialResolution;
    const std::array<size_t, 3> dims = {{nx, ny, nz }};

    plb::Array<Scalar,3> bcOffset;
    PalabosConversionTools::convertVector(m_domainInfos.origin, bcOffset);
    bcOffset *= dx;
    const plb::Array<Scalar,3> offset = m_offset + bcOffset;

    {
        const ScopedDisableFloatingPointException fpe;

        using Writer = ResultWriter::PalabosResultWriter<TransportTraits>;
        Writer::writeResults(resultsToWrite, dx, offset, dims,
                             m_results, m_simData, m_palabosData,
                             *m_plbFactory,
                             m_conf->resultProcessing,
                             m_convergenceCheck);


    }


}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::writeOutletFlux()
{
    const Scalar dx = m_conf->spatialResolution;
    const std::array<size_t, 3> dims = {{nx, ny, nz }};

    plb::Array<Scalar,3> bcOffset;
    PalabosConversionTools::convertVector(m_domainInfos.origin, bcOffset);
    bcOffset *= dx;
    const plb::Array<Scalar,3> offset = m_offset + bcOffset;

    {
        const ScopedDisableFloatingPointException fpe;

        using Writer = ResultWriter::PalabosResultWriter<TransportTraits>;
        Writer::writeOutletFlux(dx, offset, dims,
                                m_simData, m_palabosData,
                                *m_plbFactory,
                                false,
                                m_convergenceCheck);
    }
}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::saveLatticeCheckpoint() const
{

    const std::string folderName = CreateFileName::create(m_simData.getIteration());
    const boost::filesystem::path relativePath = m_conf->checkpointDir / folderName;

    if(plb::global::mpi().isMainProcessor())
    {
        FileSystemUtils::createIfNotExisting(relativePath);
    }

    HydrodynamicLatticePtr flowLattice = m_palabosData.flowLattice;
    if(flowLattice)
    {
        const boost::filesystem::path flow = relativePath / s_flowFile;
        plb::saveBinaryBlock(*flowLattice, flow.string());
    }

    const std::vector<std::string> compNames = m_simData.getCompNames();
    assert(compNames.size() ==  m_palabosData.diffLattices.size());
    for(size_t iComp = 0; iComp < compNames.size(); ++iComp)
    {
        const DiffusionLatticePtr& diffLattice = m_palabosData.diffLattices[iComp];
        const std::string& compName = compNames[iComp];
        const boost::filesystem::path compPath = relativePath / (compName +
                                                                 IPPConstants::s_dumpExtension);
        plb::saveBinaryBlock(*diffLattice, compPath.string());
    }
}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::loadLatticeCheckpoint(size_t iteration)
{
    const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<size_t>();
    MPI_Bcast(&iteration, 1, mpiType, 0, m_conf->globalComm);
    m_simData.setIteration(iteration);

    const std::string folderName = CreateFileName::create(iteration);
    const boost::filesystem::path relativePath = m_conf->checkpointDir / folderName;
    IPPCheck::assertCheck(boost::filesystem::exists(relativePath),
                          "Checkpoint folder is not existing: "
                          + relativePath.string());

    HydrodynamicLatticePtr flowLattice = m_palabosData.flowLattice;
    const boost::filesystem::path flow = relativePath / s_flowFile;

    if(m_conf->flowFunc->isEnabled())
    {
        if(boost::filesystem::exists(flow))
        {
            plb::pcout << "Loading hydrodynamic flow lattice file: " << flow.string() << std::endl;
            plb::loadBinaryBlock(*flowLattice, flow.string());
            plb::pcout << "Finished loading hydrodynamic flow lattice file" << std::endl;
        }
        else
        {
            plb::pcout << "Warning: Hydrodynamic flow lattice file not found: " << flow.string() << std::endl
                       << "Assuming non-advection simulation although was enabled via input file" << std::endl;
        }
    }

    const std::vector<std::string> compNames = m_simData.getCompNames();
    assert(compNames.size() ==  m_palabosData.diffLattices.size());
    for(size_t iComp = 0; iComp < compNames.size(); ++iComp)
    {
        const DiffusionLatticePtr& diffLattice = m_palabosData.diffLattices[iComp];
        const std::string& compName = compNames[iComp];
        const boost::filesystem::path compPath = relativePath / (compName + IPPConstants::s_dumpExtension);
        IPPCheck::assertCheck(boost::filesystem::exists(compPath),
                              "Component diffusion lattice file not found: "
                              + compPath.string());
        plb::loadBinaryBlock(*diffLattice, compPath.string());
    }
}


template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::updateRenderer(IPPRendererPtr& renderer)
{
    IPPResults::DensityVelocityResults<TransportTraits>& results = m_results.diffResultsVec.at(0);    
    typename TransportTraits::ScalarFieldSharedPtr& toPlot = results.density;

    const VirtualPalabosFieldWrapper<Scalar,
            typename TransportTraits::ScalarField, TransportTraits::dim> wrapper(*toPlot);
    renderer->updateImageData(wrapper);
}


template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::updateTransportScalarDref(const double &newPorosRef)
{
    if(newPorosRef != m_omegaCalc.getPorosRef())
    {
        m_omegaCalc.initDrefViaPorosRef(newPorosRef);
        this->updateFieldTransportScalar();
        this->updateLatticeTransportProperties();
    }

}

template<typename TransportTraits>
void PalabosTransportModule<TransportTraits>::updateLatticeTransportProperties()
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


    for(DiffusionLatticePtr& lattice : this->m_palabosData.diffLattices)
    {
        using DiffusionDesc = typename TransportTraits::template DiffusionDescriptorT<Scalar>;
        plb::setExternalScalar(*lattice, bbox,
                               DiffusionDesc::ExternalField::transportScalarBeginsAt,
                               transScalar);
    }
}


} // end of namespace IPP


#endif // PALABOSTRANSPORTMODULE_HH
