#include "transreacfactory.h"

#include <fstream>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/filesystem.hpp>

// serialization stuf for multi dim lookup
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>


#include "palabosphreeqcstate.h"
#include "ippexception.h"
#include "ippconstants.h"

#include "transportmoduleconfig.h"
#include "transportmodulefactory.h"
#include "transportmodule.h"
#include "cellidpassingdissolveonlyfunc.h"



#include "reactionmodule.h"
#include "reactionmoduleconfig.h"
#include "phreeqcreactionmodulefactory.h"
#include "phreeqcreactionmoduledata.h"


// chemistry optimizations
#include "configurephreeqccellsactivity.h"
#include "reactionmoduleoptimization.h"
#include "isinsidevolume.h"
#include "significantchangeprediction.h"



// dynamic config functions
#include "abstractgeometricsicalc.h"
#include "abstractdissolveonlycalc.h"
#include "cellindexcorrection.h"

// cellID -> coordinate wrapper
#include "cellidconverttosifunc.h"

#include "filesystemutils.h"


#include "latticeboltzmannunitconversion.h"

// mpi tools:
#include "mpitools.h"
#include "mpimanager.h"
#include "fielddecomposition.h"
#include "simulationexchangedata.h"
#include "ippstream.h"

#include "abstractsicalc.h"

namespace IPP
{

template<typename Func, typename WrappedFunc>
static const Func* createWrapperFunc(const WrappedFunc calc,
                                     const size_t nx, const size_t ny, const size_t nz)
{
    assert(calc);

    const ArrayDimensionConvert* indexConv = new ArrayDimensionConvert(nx, ny, nz);
    const Func* indexConvWrappedCalc = new Func(indexConv, *calc);
    return indexConvWrappedCalc;

}

static BoundaryConditionsPtr convertAdvectiveBC(const ConfigBoundaryConditions& bc)
{
    BoundaryConditionsPtr result(new BoundaryConditions());

    result->periodicDims = bc.periodicDims;

    // copy advection BC which is directly compatible
    for(size_t iBC = 0; iBC < bc.advectiveBC.size(); ++iBC)
    {
        const ConfigBoundaryConditions::AdvectiveBoundaryConditionData& advBC = bc.advectiveBC[iBC];
        result->advectiveBC.push_back(BoundaryConditions::BoundaryConditionData());
        BoundaryConditions::BoundaryConditionData& resultBC = result->advectiveBC.back();
        resultBC.domain = advBC.domain;
        resultBC.type = advBC.type;
        resultBC.density = advBC.density;
        resultBC.velocity = advBC.velocity;
    }

    return result;
}

static void addDiffEffInertBCs(const size_t nxInput, const size_t nyInput, const size_t nzInput,
                               const IPPConfig::Results::DiffEffInfos& infos,
                               BoundaryConditions::ElementNameToBCVec& elemToBCVec)
{
    IPPCheck::assertCheck(nxInput >= 1);
    IPPCheck::assertCheck(nyInput >= 1);
    IPPCheck::assertCheck(nzInput >= 1);

    const int nx = static_cast<int>(nxInput);
    const int ny = static_cast<int>(nyInput);
    const int nz = static_cast<int>(nzInput);


    for(const IPPConfig::Results::DiffEffInfo& info : infos)
    {
        const std::string name = info.tracerName;
        const size_t dim = info.dim;

        BoundaryConditions::BCVec& bcVec = elemToBCVec.at(name);


        bcVec.push_back(BoundaryConditions::BoundaryConditionData());
        bcVec.push_back(BoundaryConditions::BoundaryConditionData());

        auto it = bcVec.end() - 2;

        BoundaryConditions::BoundaryConditionData& bcData0 = *it;
        bcData0.type = BCT_DensityDirichlet;
        bcData0.density = 1.0;
        bcData0.domain.range.lower = {{ 0, 0, 0 }};

        ++it;
        BoundaryConditions::BoundaryConditionData& bcData1 = *it;
        bcData1.type = BCT_DensityDirichlet;
        // bcData1.density = 1.0e-12;
        bcData1.density = 0.0;
        bcData1.domain.range.upper = {{ nx-1, ny-1, nz-1 }};

        // dim -> axis
        // 0 -> x
        // 1 -> y
        // 2 -> z
        switch(dim)
        {
        case 0:
        {
            bcData0.domain.position = BoundaryConditionDomain::BP_left;
            IPPBox3DInt& range0 = bcData0.domain.range;
            range0.upper = {{ 0, ny-1, nz-1 }};

            bcData1.domain.position = BoundaryConditionDomain::BP_right;
            IPPBox3DInt& range1 = bcData1.domain.range;
            range1.lower = {{ nx-1, 0, 0 }};

            break;
        }

        case 1:
        {
            bcData0.domain.position = BoundaryConditionDomain::BP_bottom;
            IPPBox3DInt& range0 = bcData0.domain.range;
            range0.upper = {{ nx-1, 0, nz-1 }};

            bcData1.domain.position = BoundaryConditionDomain::BP_top;
            IPPBox3DInt& range1 = bcData1.domain.range;
            range1.lower = {{ 0, ny-1, 0 }};

            break;
        }

        case 2:
        {
            bcData0.domain.position = BoundaryConditionDomain::BP_back;
            IPPBox3DInt& range0 = bcData0.domain.range;
            range0.upper = {{ nx-1, ny-1, 0 }};

            bcData1.domain.position = BoundaryConditionDomain::BP_front;
            IPPBox3DInt& range1 = bcData1.domain.range;
            range1.lower = {{ 0, 0, nz-1 }};

            break;
        }

        default:
        {
            throw std::runtime_error("unknown dimension for diffusion coeff calculation: "
                                     + std::to_string(dim));
        }
        }

    }
}

void fixBCposition(const DomainInfos& domainInfos,
                   BoundaryConditionDomain& bc)
{
    switch(bc.position)
    {
        case BoundaryConditionDomain::BP_left:
        {
            // x values
            bc.range[0] = 0;
            bc.range[1] = 0;

            // correct y and z values
            bc.range.lower[1] -= domainInfos.origin[1];
            bc.range.upper[1] -= domainInfos.origin[1];
            bc.range.lower[2] -= domainInfos.origin[2];
            bc.range.upper[2] -= domainInfos.origin[2];

            break;
        }

        case BoundaryConditionDomain::BP_right:
        {
            // x values
            bc.range[0] = domainInfos.totalSize[0] - 1;
            bc.range[1] = domainInfos.totalSize[0] - 1;

            // correct y and z values
            bc.range.lower[1] -= domainInfos.origin[1];
            bc.range.upper[1] -= domainInfos.origin[1];
            bc.range.lower[2] -= domainInfos.origin[2];
            bc.range.upper[2] -= domainInfos.origin[2];

            break;
        }

        case BoundaryConditionDomain::BP_bottom:
        {
            // y values
            bc.range[2] = 0;
            bc.range[3] = 0;

            // correct x and z values
            bc.range.lower[0] -= domainInfos.origin[0];
            bc.range.upper[0] -= domainInfos.origin[0];
            bc.range.lower[2] -= domainInfos.origin[2];
            bc.range.upper[2] -= domainInfos.origin[2];

            break;
        }

        case BoundaryConditionDomain::BP_top:
        {
            // x values
            bc.range[2] = domainInfos.totalSize[1] - 1;
            bc.range[3] = domainInfos.totalSize[1] - 1;

            // correct x and z values
            bc.range.lower[0] -= domainInfos.origin[0];
            bc.range.upper[0] -= domainInfos.origin[0];
            bc.range.lower[2] -= domainInfos.origin[2];
            bc.range.upper[2] -= domainInfos.origin[2];

            break;
        }

        case BoundaryConditionDomain::BP_back:
        {
            // z values
            bc.range[4] = 0;
            bc.range[5] = 0;

            // correct x and y values
            bc.range.lower[0] -= domainInfos.origin[0];
            bc.range.upper[0] -= domainInfos.origin[0];
            bc.range.lower[1] -= domainInfos.origin[1];
            bc.range.upper[1] -= domainInfos.origin[1];

            break;
        }

        case BoundaryConditionDomain::BP_front:
        {
            // z values
            bc.range[4] = domainInfos.totalSize[2] - 1;
            bc.range[5] = domainInfos.totalSize[2] - 1;

            // correct x and y values
            bc.range.lower[0] -= domainInfos.origin[0];
            bc.range.upper[0] -= domainInfos.origin[0];
            bc.range.lower[1] -= domainInfos.origin[1];
            bc.range.upper[1] -= domainInfos.origin[1];

            break;
        }

        default:
        {
            throw std::runtime_error("unknown boundary position: " + std::to_string(bc.position));
            break;
        }
    }
}

IPPBox3DInt findBoundaryBox(const IPPVector3DInt& size,
                            const BoundaryConditionDomain::BoundaryPosition& pos)
{
    IPPBox3DInt box({{0, 0, 0}}, {{size[0]-1, size[1]-1, size[2]-1}});

    switch(pos)
    {
    case BoundaryConditionDomain::BP_left:
        box[1] = box[0];
        break;

    case BoundaryConditionDomain::BP_right:
        box[0] = box[1];
        break;

    case BoundaryConditionDomain::BP_bottom:
        box[3] = box[2];
        break;

    case BoundaryConditionDomain::BP_top:
        box[2] = box[3];
        break;

    case BoundaryConditionDomain::BP_back:
        box[5] = box[4];
        break;

    case BoundaryConditionDomain::BP_front:
        box[4] = box[5];
        break;

    default:
    {
        throw std::runtime_error("unknown boundary position: " + std::to_string(pos));
        break;
    }
    }

    return box;
}

static void determineDecompIndexCorr(const IPPVector3DInt& totalSize,
                                     const FieldDecomposition& decomp,
                                     std::vector<size_t>& globaldToLocalIndices)
{

    globaldToLocalIndices.resize(decomp.getGlobalNumberCells(), -1);
    const ArrayDimensionConvert globalConv(totalSize);

    size_t localIndex = 0;

    const std::vector<IPPBox3DLong>& domains = decomp.getDomains();
    for(size_t iDomain = 0; iDomain < domains.size(); ++iDomain)
    {
        const IPPBox3DLong& domain = domains[iDomain];
        // const IPPVector3DLong& origin = domain.lower;


        for(long x = domain.lower[0]; x <= domain.upper[0]; ++x)
        {
            for(long y = domain.lower[1]; y <= domain.upper[1]; ++y)
            {
                for(long z = domain.lower[2]; z <= domain.upper[2]; ++z)
                {
                    const IPPVector3DLong globalCoord = {{ x, y, z }};
                    const size_t globalIndex = globalConv.calcIndex(globalCoord);

                    globaldToLocalIndices[globalIndex] = localIndex;

                    ++localIndex;
                }
            }
        }
    }
}

static void addPointsAsIndex(const ArrayDimensionConvert& idConv, const IPPVector3DInt& origin,
                             const IPPBox3DInt& disabledBox, std::vector<size_t>& disabledCells)
{
    for(int x = disabledBox.lower[0]; x <= disabledBox.upper[0]; ++x)
    {
        for(int y = disabledBox.lower[1]; y <= disabledBox.upper[1]; ++y)
        {
            for(int z = disabledBox.lower[2]; z <= disabledBox.upper[2]; ++z)
            {
                const IPPVector3DInt posGlobal = {{ x, y, z }};
                // assert(ownDomain.contains(posGlobal));

                const IPPVector3DInt localPos = posGlobal - origin;
                const size_t id = idConv.calcIndex(localPos);
                disabledCells.push_back(id);
            }
        }
    }
}

static void disableChemistryAtInert(const DomainInfos& domainInfos,
                                    const FieldDecomposition& decomp,
                                    const std::vector<char>& inertCells,
                                    const std::vector<IPPBox3DInt>& disabledDomains,
                                    ConfigurePhreeqcCellsActivity& configActivity)
{
    std::vector<size_t> disabledCells;

    const IPPBox3DLong& ownDomainLong = decomp.getOwnDomain();
    IPPBox3DInt ownDomain;
    for(size_t i = 0; i < 6; ++i)
    {
        ownDomain[i] = ownDomainLong[i];
    }

    const IPPVector3DInt& origin = ownDomain.lower;
    const IPPVector3DInt ownSize = ownDomain.getSize();
    const ArrayDimensionConvert idConv(ownSize);

    std::stringstream ss;
    ss << MPIManager::getInstance().getRank() << ": owndomain "
       << ownDomain << std::endl << " envelops: " << std::endl;
    for(const BoundaryConditionDomain::BoundaryPosition& envelop : domainInfos.boundaryEnvelops)
    {
        const IPPBox3DInt box = findBoundaryBox(domainInfos.totalSize, envelop);
        const IPPBox3DInt inter = intersection(box, ownDomain);

        ss << "\tboundary: " << box << "\t intersect: " << inter << std::endl;

        addPointsAsIndex(idConv, origin, inter, disabledCells);
    }

    for(const IPPBox3DInt& box : disabledDomains)
    {
        IPPBox3DInt disabledBox = box;
        disabledBox.lower -= domainInfos.origin;
        disabledBox.upper -= domainInfos.origin;
        const IPPBox3DInt inter = intersection(disabledBox, ownDomain);
        addPointsAsIndex(idConv, origin, inter, disabledCells);
    }

    // ss << MPIManager::getInstance().getRank() << ": " << "inert Cells: ";


    for(size_t iCell = 0; iCell < inertCells.size(); ++iCell)
    {
        if(inertCells[iCell])
        {
            // ss << iCell << " ";
            disabledCells.push_back(iCell);
        }
    }

    // ss << std::endl;

    // std::cout << ss.str();

    configActivity.setStaticDisabledCells(disabledCells);
}

static SignificantChangePrediction* createChemPredict(const IPPConfig::ChemistryPrediction& config)
{

    std::ifstream file(config.filename);

    if(file.is_open())
    {
        SignificantChangePrediction* predict = new SignificantChangePrediction(config.significantPorosDiff);
        SignificantChangePrediction::Interpolation& data = predict->getLookup();

        boost::archive::binary_iarchive ia(file);
        ia >> data;

        return predict;
    }
    else
    {
        throw std::runtime_error("could not open interpolation file: " + config.filename);
    }

    return nullptr;
}

class SIFuncWrapper : public CellID_SI_Func
{
public:
    explicit SIFuncWrapper(AbstractSICalcPtr& wrapped)
        : m_wrapped(wrapped)
    {
        assert(m_wrapped);
    }

    virtual void evaluate(const size_t iCell, const iterator& begin, const iterator& end) const
    {
        assert(m_wrapped);
        m_wrapped->evaluate(iCell, begin, end);
    }

private:
    // needs to save a copy of the smart pointer
    AbstractSICalcPtr m_wrapped;

};

TransReacState *TransReacFactory::create(const IPPConfig &conf)
{
    const boost::mpi::communicator& comm = conf.globalComm;

    PalabosPhreeqcState* state = new PalabosPhreeqcState;
    state->data.reset(new SimulationExchangeData(comm));
    state->data->setIteration(0);
    state->fluxConvergenceTol = conf.fluxConvergenceTolLow;


    using namespace boost::filesystem;
    const path relativePath(conf.outputDir);
    const path outPath = absolute(relativePath);

    const path resultOutPath = outPath / IPPConstants::s_resultOutput;
    const path chemOutPath = outPath / IPPConstants::s_chemOutput;
    const path checkpointOutPath = outPath / IPPConstants::s_checkpointOutput;



    TransportModuleConfigPtr transConf(new TransportModuleConfig(comm, conf,
                                                                 resultOutPath,
                                                                 checkpointOutPath,
                                                                 *(state->data),
                                                                 state->fluxConvergenceTol,
                                                                 state->isConverged));

    state->transModule = TransportModuleFactory::create(transConf);
    assert(state->transModule);

    const DomainInfos& domainInfos = state->transModule->getDomainInfos();
    const size_t nx = domainInfos.totalSize[0];
    const size_t ny = domainInfos.totalSize[1];
    const size_t nz = domainInfos.totalSize[2];

    pcout << "total system size: " << nx << " x " << ny << " x " << nz << std::endl;



    const FieldDecomposition* decomp = state->transModule->getDecomposition();
    state->data->setDecomp(decomp);


    const double dtTmp = LatticeBoltzmannUnitConversion::calcTimeStep(conf);
    state->data->setDtCurr(dtTmp);
    const double& dt = state->data->getDtCurr();


    const path outFilePath = outPath / "out";
    const std::string outFile = outFilePath.string();

    bool success;
    if(comm.rank() == 0)
    {
        FileSystemUtils::createIfNotExisting(outPath);
        FileSystemUtils::createIfNotExisting(chemOutPath);
        FileSystemUtils::createIfNotExisting(resultOutPath);
        FileSystemUtils::createIfNotExisting(checkpointOutPath);

        std::cout << "writing to: " << outPath.string() << std::endl;

        std::ofstream file(outFile);
        if(file.is_open())
        {
            const double Dref = conf.diffusionCoefReference;
            file << "diff coef ref:\t" << Dref << "\tm2/s" << std::endl;
            file << "tau ref:\t" << conf.tauRef << std::endl;
            file << "porosity ref:\t" << conf.porosRef << std::endl;
            file << "porosity low:\t" << conf.porosLow << std::endl;
            file << "time step:\t" << dt << "\ts" << std::endl;
            success = true;
        }
        else
        {
            success = false;
        }
    }

    MPI_Bcast(&success, 1, MPI_CXX_BOOL, 0, comm );
    if(success == false)
    {
        throw std::runtime_error("could not open file: " + outFile);
    }

    if(comm.rank() == 0)
    {
        const std::vector<IPPBox3DLong>& domains = decomp->getDomains();
        std::cout << "MPI decompositions:\n";
        for(size_t iDecomp = 0; iDecomp < domains.size(); ++iDecomp)
        {
            const IPPBox3DLong& domain = domains[iDecomp];
            std::cout << iDecomp << ":\t" << domain << std::endl;
        }
    }



    PhreeqcReactionModuleData& phreeqcData = *state->phreeqcData;

    AbstractDissolveOnlyCalcPtr dissOnlyFunc = conf.dissolveOnlyFunc;
    const std::vector<double>& capillPorosity = state->data->getCapillaryPorosity();
    const std::vector<CellNeighborInfo>& neighInfos = state->data->neighInfos;
    const std::vector<double>& SI = state->data->getSaturationIndices();
    CellIDPassingDissolveOnlyFunc* dissolveOnlyFunc =
            new CellIDPassingDissolveOnlyFunc(dissOnlyFunc,
                                              capillPorosity,
                                              neighInfos,
                                              SI,
                                              phreeqcData.phaseToNuclPhaseID,
                                              phreeqcData.nuclPhaseToMonomerID,
                                              phreeqcData.monomerConc,
                                              state->data->getPhaseNameToInfos());

    phreeqcData.dissolveOnlyFunc = dissolveOnlyFunc;

    AbstractSICalcPtr siCalc = conf.saturationIndexCalc;
    if(siCalc)
    {
        phreeqcData.siFunc = new SIFuncWrapper(siCalc); // createWrapperFunc<CellIDConvertToSIFunc>(siCalc, nx, ny, nz);
    }



    phreeqcData.optim = new ReactionModuleOptimization;

    if(conf.chemPredict.enabled)
    {
        phreeqcData.optim->changePredict = createChemPredict(conf.chemPredict);
        phreeqcData.optim->forceRecalcFreq = conf.chemPredict.forceRecalcFreq;
    }


    // correct BCs by new offset / total size
    {
        ConfigBoundaryConditions::AdvectiveBCVec& advBC = conf.boundaryConditions->advectiveBC;
        for(BoundaryConditionData& bc : advBC)
        {
            fixBCposition(domainInfos, bc.domain);
        }

        ConfigBoundaryConditions::DiffusiveBCVec& diffBC = conf.boundaryConditions->diffusiveBC;
        for(BoundaryConditionData& bc : diffBC)
        {
            fixBCposition(domainInfos, bc.domain);
        }
    }


    BoundaryConditionsPtr bc = convertAdvectiveBC(*conf.boundaryConditions);
    std::shared_ptr<SimulationExchangeData>& exData = state->data;
    exData->setBoundaryConditions(bc);

    // diffusive BCs are converted in reaction module since concentrations needs to be defined


    // cell ids have changed due to BC extensions
    const IPPVector3DInt oldSize = {{ (int)conf.nx, (int)conf.ny, (int)conf.nz }};
    std::vector<size_t> globaldToLocalIndices;
    determineDecompIndexCorr(domainInfos.totalSize, *decomp, globaldToLocalIndices);
    const CellIndexCorrection indexCorr(oldSize, domainInfos.totalSize,
                                        domainInfos.origin,
                                        decomp->getOwnDomain(),
                                        globaldToLocalIndices);

    NonLocalOperations& nonLocalOperation = state->transModule->getNonLocalOperations();

    const ReactionModuleConfig reacConf(conf, indexCorr,
                                        exData, bc->diffusiveBC,
                                        exData->getInertComposition(),
                                        phreeqcData,
                                        chemOutPath,
                                        checkpointOutPath,
                                        nonLocalOperation);

    const size_t nxyz = nx * ny * nz;
    state->reacModule = PhreeqcReactionModuleFactory::create(comm,
                                                             conf.isChemistryEnabled,
                                                             reacConf,
                                                             nxyz);

    state->data->initSyncEarly();



    const double cellVol = conf.is3DSimulation ?
                               conf.spatialResolution * conf.spatialResolution * conf.spatialResolution
                             :
                               conf.spatialResolution * conf.spatialResolution;

    const double cellArea = conf.is3DSimulation ?
                                conf.spatialResolution * conf.spatialResolution
                              :
                                conf.spatialResolution;

    assert(cellVol > 0.0);
    assert(dt > 0.0);
    dissolveOnlyFunc->init(cellVol, cellArea, dt);



    addDiffEffInertBCs(nx, ny, nz, conf.results.diffEffInfos, bc->diffusiveBC);

    for(const IPPConfig::Results::DiffEffInfo& info : conf.results.diffEffInfos)
    {
        std::cout << "effective diffusion tracer: " << info.dim << "\t" << info.tracerName << std::endl;
        exData->setInertTracer(info.dim, info.tracerName);
    }
    state->data->initSyncLate();



    const std::vector<std::string>& phaseNames = state->data->getPhaseNames();

    // initialize target SI and dissolve only functions after initialization of reaction module
    // phase info data are only available afterwards

    if(siCalc)
    {
        const PhaseNameToInfos& phaseInfos = state->data->getPhaseNameToInfos();
        const std::vector<double>& poros = state->data->getPorosity();
        siCalc->init(phaseNames, phaseInfos, &poros);
    }

    // dissolve only behaviour must be reset to initial state when porosity is higher than threshhold again
    // PhaseNameToDissPrecipBehaviour does contain information for specific phases only
    // AbstractDissolveOnlyCalc needs phase index based informations


    std::vector<DissPrecipBehaviour> behavVec(phaseNames.size());

    for(size_t iPhase = 0; iPhase < phaseNames.size(); ++iPhase)
    {
        const std::string& phaseName = phaseNames[iPhase];
        const DissPrecipBehaviour& behav = conf.dissPrecBehav->getBehaviour(phaseName);
        behavVec[iPhase] = behav;
    }

    dissOnlyFunc->setPrecipDissolveOnlyBehaviour(behavVec);


    AbstractGeometrySync* geomSync = state->reacModule->getGeomSync();
    state->transModule->setGeomSync(geomSync);


    if(conf.isChemistryEnabled)
    {
        ConfigurePhreeqcCellsActivity* configActivity = state->reacModule->getConfigCellsActivity();
        assert(configActivity);

        const std::vector<char>& inertCells = exData->getInertSolidCells();
        disableChemistryAtInert(domainInfos, *decomp, inertCells,
                                conf.chemDisabledBoxes, *configActivity);

    }

    // init transport data fields
    state->transModule->init();

    return state;

}

}
