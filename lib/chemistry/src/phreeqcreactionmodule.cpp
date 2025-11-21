#include "phreeqcreactionmodule.h"

#include <assert.h>
#include <fstream>
#include <iomanip>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <numeric>
#include <regex>


#include <boost/algorithm/string.hpp>
#include <boost/mpi.hpp>
#include <boost/filesystem/operations.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "IPhreeqcPhast.h"
#include "cellindexcorrection.h"
#include "disablephreeqccells.h"
#include "reactionmoduleconfig.h"

#include "serialtomultidatasync.h"
#include "multitoserialdatasync.h"

#include "mpi_types.h"
#include "mpitools.h"

#include "simulationexchangedata.h"
#include "ippexception.h"
#include "phreeqcinputfilegenerator.h"
#include "phreeqcprerunhandler.h"
#include "phreeqcconstants.h"
#include "cellid_si_func.h"
#include "phreeqcreactionmoduledata.h"
#include "phreeqcmodify.h"
#include "bench_tools.h"
#include "retrievephaseamounts.h"
#include "isenabledincontainer.h"

#include "handledisabledreactioncells.h"
#include "enablephreeqccells.h"
#include "disablephreeqccells.h"
#include "reactionmoduleoptimization.h"
#include "configurephreeqccellsactivity.h"
#include "handleinterfaceswitch.h"

#include "porositycalcutils.h"
#include "runphreeqccells.h"

#include "phreeqcdump.h"
#include "createfilename.h"
#include "filesystemutils.h"
#include <abstractmultiscalediffusioncalc.h>
#include "phreeqcinitrun.h"
#include "initporosity.h"

#include "phreeqccalcporosity.h"
#include "phreeqccalcdiffcoefs.h"

#include "phreeqcsetfileprefix.h"
#include "geometrytools.h"
#include "ippstream.h"
#include "isinsidevolume.h"
#include "setupphreeqcverbosity.h"
#include "phaseidperelement.h"
#include "indexhelper.h"
#include "phreeqcrmhelper.h"
#include "fixcellsidsfornewfielddecomposition.h"
#include "porosityhelper.h"
#include "findzeroporosityroot.h"
#include "phreeqcglobalinitdata.h"

#include "mpidelegate.h"

#include "celltotalsdiff.h"
#include "retrieveauxvalues.h"
#include "scopedfloatingpointexception.h"

#include "nonlocaloperations.h"
#include "abstractdistributeexcesstotals.h"
#include "abstractcalcinterfaceproperties.h"
#include "phreeqcgeometrysync.h"
#include "calcporositychangefactor.h"
#include "printphaseinfos.h"
#include "processmemoryusage.h"
#include "cellidpassingdissolveonlyfunc.h"
#include "ippconstants.h"
#include "simpleporosityinfos.h"

namespace IPP
{

static void calcNewDiffCoeffs(const AbstractMultiScaleDiffusionCalc& diffCalc,
                              const std::vector<SimplePorosityInfosPtr>& porosInfos,
                              std::vector<double>& diffCoefs,
                              const bool verbose)
{
    const auto calcDiffTimer = std::chrono::steady_clock::now();

    PhreeqcCalcDiffCoefs::calc(diffCalc, porosInfos, diffCoefs);

    if(verbose)
    {
        const auto durationCalc = CalcDuration::calc(calcDiffTimer);
        pcout << "\t-> calc diff coef took:\t" << durationCalc.count() << " ms" << std::endl;
    }

}


static const std::string s_dumpFilename = "phreeqc.dmp";
static const std::string s_decompFilename = "decomp";
static const std::string s_currPorosFilename = "porosity_curr";
static const std::string s_oldPorosFilename = "porosity_old";
static const std::string s_nodeInfosFilename = "node_infos";
static const std::string s_reacAuxFilename = "rm_aux";


PhreeqcReactionModule::PhreeqcReactionModule(const size_t nCellsGlobal,
                                             const boost::mpi::communicator &comm,
                                             const ReactionModuleConfig& conf)
    : m_comm(comm)
    , m_isPrintChemistryOn(conf.ippConfig.write_chem)
    , m_phreeqc(nCellsGlobal, comm)
    , m_simData(conf.data)
    , m_reacData(*m_simData, *conf.ippConfig.porosCalc, *conf.ippConfig.multiScaleDiffusionCalc)
    , m_chemDir(conf.chemOutDir)
    , m_dumpDir(conf.dumpOutDir)
    , m_nonLocalOperations(conf.nonLocalOperation)
    , m_phreeqcData(conf.phreeqcData)
    , m_configActivity(nullptr)
    , m_indexConvLocal(conf.data->getDecomp()->getOwnDomain().getSize())
    , m_verbose(conf.ippConfig.verbose)
{
    MPI_CHECK_SYNC;
    this->init(nCellsGlobal, conf);
}

PhreeqcReactionModule::~PhreeqcReactionModule()
{
    MPI_CHECK_SYNC;
}


void PhreeqcReactionModule::init(const size_t nCellsGlobal, const ReactionModuleConfig& conf)
{

    this->preInit(conf);


    PhreeqcGlobalInitData initData;
    phreeqcMPIDispatch(m_phreeqc, [&]() { this->initRoot(nCellsGlobal, conf, initData); });

    this->pushInitData(initData);
    this->disableLowPorosCells(m_simData->getPorosity(), conf.ippConfig.porosLow);

    phreeqcMPIDispatch(m_phreeqc, [&]() { this->initRootRun(conf, initData); });


    boost::mpi::broadcast(m_comm, initData.genericCellID, 0);
    boost::mpi::broadcast(m_comm, initData.featureMask, 0);


    this->initCells(nCellsGlobal, conf, initData);
    this->initConcentrations();


    const FieldDecomposition& decomp = *m_simData->getDecomp();
    const size_t nCells = decomp.getLocalNumberOfCells();
    const size_t nPhases = m_simData->getPhaseNames().size();
    // check if only eq phases are relevant for SI
    m_simData->getSaturationIndices().resize(nPhases * nCells);

    this->convergePorosity();


    if(m_phreeqcData.dissolveOnlyFunc->needsNeighInfos())
    {
        m_simData->neighInfos.resize(nCells);
    }

    const IPPVector3DLong& origin = decomp.getOwnDomain().lower;
    AbstractGeometrySync::NeighborInfos neighInfos(m_phreeqcData.dissolveOnlyFunc->getLowerPorosityThresh(),
                                                   m_simData->neighInfos);
    m_geomSync = std::make_unique<PhreeqcGeometrySync>(origin, m_indexConvLocal,
                                                       m_initialData,
                                                       m_nodeInfos,
                                                       neighInfos);

#ifdef IPP_DEBUG_ENABLED
    const std::vector<char> initEnabled(decomp.getLocalNumberOfCells(), false);
    m_simData->debugData.enabledCellsSteps.resize(2, initEnabled);
#endif

}

void PhreeqcReactionModule::preInit(const ReactionModuleConfig &conf)
{
    MPI_CHECK_SYNC;

    const IPPConfig& ippConf = conf.ippConfig;

    const int myRank = ippConf.globalComm.rank();
    pcout << "starting preInit. this is rank: " << myRank << std::endl;

    if(m_isPrintChemistryOn)
    {
        PhreeqcRMHelper::muteConsoleOutput(m_phreeqc, m_chemDir);
        m_phreeqc.SetScreenOn(true);
        m_phreeqc.setPrintChemistryOnLocal(true, true, true);
    }
    else
    {
        PhreeqcRMHelper::muteConsoleOutput(m_phreeqc, m_chemDir);
        m_phreeqc.setPrintChemistryOnLocal(false, false, false);
    }



    const double cellVol = ippConf.is3DSimulation ? ippConf.spatialResolution
                                                    * ippConf.spatialResolution
                                                    * ippConf.spatialResolution
                                                  : ippConf.spatialResolution
                                                    * ippConf.spatialResolution;
    m_reacData.invTotalCellVolume = 1.0 / cellVol;


    std::vector<std::string> compNames;
    std::vector<std::string>& phaseNames = m_simData->getPhaseNames();
    PhaseNameToInfos& phaseInfos = m_simData->getPhaseNameToInfos();
    PhaseIdPerElement phaseIdsPerElement;
    PhreeqcPrerunHandler::getAuxDataFromPhreeqC(ippConf, conf.chemOutDir,
                                                PhreeqcConstants::s_calcH2OSeperately,
                                                compNames, phaseNames,
                                                phaseInfos, phaseIdsPerElement,
                                                conf.diffBCs, conf.inertComposition);

    m_simData->setCompNames(compNames);

    assert(phaseNames.size() == phaseInfos.size());
    if(MPIManager::getInstance().isMainProc())
    {
        PrintPhaseInfos::print(compNames, phaseNames, phaseInfos);
    }

    const double& tolerance = ippConf.chemistryTolerance;
    const SignificantChangePrediction* changePredict = conf.phreeqcData.optim->changePredict;
    m_configActivity.reset(new ConfigurePhreeqcCellsActivity(tolerance,
                                                             conf.phreeqcData.optim->forceRecalcFreq,
                                                             changePredict));

    const FieldDecomposition* decomp = m_simData->getDecomp();
    m_phreeqc.init(decomp);


    const size_t nCellsLocal = decomp->getLocalNumberOfCells();
    const size_t nPhases = phaseNames.size();
    const size_t nComps = compNames.size();
    const size_t nPrimaryComps = nComps - PhreeqcConstants::s_primaryCompsBegin;

    std::vector<double>& precipField = m_simData->getPrecipField();
    precipField.resize(nCellsLocal * nPhases, -1.0);



    m_enabledCells.resize(nCellsLocal, true);

    m_currTargetSImap.resize(nCellsLocal * nPhases, 0.0);
    m_currDissolveOnlyCells.resize(nCellsLocal, false);
    m_simData->solidFlags.resize(nCellsLocal, IPPConstants::s_isNotSteadyFlag);

    m_nodeInfos.nonPermNonInterfaceCells.resize(nCellsLocal, false);
    m_nodeInfos.interfaceCells.resize(nCellsLocal, false);
    m_nodeInfos.nInterfaceNodesGlobal = 0;


    m_configActivity->init(nCellsLocal, nComps, nPrimaryComps,
                           nPhases, m_preReacConc,
                           phaseIdsPerElement,
                           m_nodeInfos.nonPermNonInterfaceCells,
                           m_nodeInfos.interfaceCells);

    if(changePredict)
    {
        throw std::runtime_error("currently unsupported");
    } else {
        m_enableCellsFunc.reset(new EnablePhreeqcCells(m_phreeqc, m_enabledCells));
        m_disableCellsFunc.reset(new DisablePhreeqcCells(m_phreeqc, m_enabledCells));
    }

    const std::vector<IPPConfig::Results::NameCommand> &auxStrVec = conf.ippConfig.results.auxData;
    AuxDataVec& auxDataVec = m_simData->getAuxData();
    auxDataVec.resize(auxStrVec.size());
    for(size_t i = 0; i < auxStrVec.size(); ++i)
    {
        const IPPConfig::Results::NameCommand& nameComm = auxStrVec[i];
        AuxDataName& auxData = auxDataVec[i];
        auxData.name = nameComm.name;
    }



    m_reacData.isPermeablePhase.resize(phaseNames.size());
    for(size_t iPhase = 0; iPhase < phaseNames.size(); ++iPhase)
    {
        const std::string phaseName = phaseNames[iPhase];
        const PhaseInfo& info = phaseInfos.at(phaseName);
        m_reacData.isPermeablePhase[iPhase] = info.isPermeable;
    }

    m_nonLocalOperations.setnComps(nComps);


    {
        std::vector<std::string> nuclPhases;
        std::vector<std::string> monomers;
        conf.ippConfig.dissolveOnlyFunc->findNucleationPhases(nuclPhases, monomers);

        const size_t nNuclPhases = nuclPhases.size();

        m_simData->nucleationPhasesVolumeFractions.resize(nNuclPhases,
                                                          std::vector<double>(nCellsLocal, -1.0));

        m_phreeqcData.phaseToNuclPhaseID.resize(nPhases, -1);
        m_phreeqcData.nuclPhaseToMonomerID.resize(nNuclPhases, -1);

        std::vector<std::string> uniqueMonomers = monomers;
        std::sort(uniqueMonomers.begin(), uniqueMonomers.end());
        uniqueMonomers.erase(std::unique(uniqueMonomers.begin(), uniqueMonomers.end()),
                             uniqueMonomers.end());

        for(size_t iNuclPhase = 0; iNuclPhase < nNuclPhases; ++iNuclPhase)
        {
            const std::string& phaseName = nuclPhases[iNuclPhase];
            const std::vector<std::string>::const_iterator phaseIt =
                    std::find(phaseNames.begin(), phaseNames.end(), phaseName);
            IPPCheck::assertCheck(phaseIt != phaseNames.end());
            const size_t iPhase = phaseIt - phaseNames.cbegin();
            assert(m_phreeqcData.phaseToNuclPhaseID.size() > iPhase);
            m_phreeqcData.phaseToNuclPhaseID[iPhase] = iNuclPhase;

            const std::string& monomerName = monomers[iNuclPhase];
            const std::vector<std::string>::const_iterator monomerIt =
                    std::find(uniqueMonomers.begin(), uniqueMonomers.end(), monomerName);
            IPPCheck::assertCheck(monomerIt != monomers.end());
            const size_t iMonomer = monomerIt - uniqueMonomers.cbegin();
            assert(m_phreeqcData.nuclPhaseToMonomerID.size() > iNuclPhase);
            m_phreeqcData.nuclPhaseToMonomerID[iNuclPhase] = iMonomer;
        }

        m_phreeqcData.monomerConc.resize(uniqueMonomers.size() * nCellsLocal, -1.0);
    }

}

void PhreeqcReactionModule::initMapping(const std::vector<int>& grid2chem,
                                        const std::vector<int>& ic)
{
    assert(MPIManager::getInstance().isMainProc());

    //    for(size_t i = 0; i < grid2chem.size(); ++i)
    //    {
    //        std::cout << i << "\t" << grid2chem[i] << std::endl;
    //    }

    //    for(size_t i = 0; i < ic.size(); ++i)
    //    {
    //        std::cout << i << "\t" << ic[i] << std::endl;
    //    }

    IRM_RESULT status;
    status = m_phreeqc.CreateMapping(const_cast<std::vector<int>&>(grid2chem));
    IPPCheck::assertCheck(status == IRM_OK);

    status = m_phreeqc.InitialPhreeqc2Module(const_cast<std::vector<int>&>(ic));
    if (status != IRM_OK)
    {
        // retrieve error messages if needed
        std::cerr << m_phreeqc.GetErrorString();
        throw PhreeqcRMStop();
    }

}

void PhreeqcReactionModule::runCells()
{
    MPI_CHECK_SYNC;

    this->setTime();


    const PhreeqcSetFilePrefix setPrefix(m_phreeqc, m_chemDir,
                                         m_simData->getIteration());

    const SetupPhreeqcVerbosity setupVerbosity(m_phreeqc, m_isPrintChemistryOn);
    const RunPhreeqcCells runCells(m_phreeqc, setPrefix, setupVerbosity);


    {
        auto start = std::chrono::steady_clock::now();

        runCells.run();

        if(m_verbose)
        {
            auto duration = CalcDuration::calc(start);
            pcout << "\t-> chemistry (runCells) took:\t"
                  << duration.count() << " ms" << std::endl;
        }
    }

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::swapToOld()
{
    std::vector<double>& porosityPhysical = m_simData->getPorosity();
    m_oldPorosPhysical.swap(porosityPhysical);
}

void PhreeqcReactionModule::initPrecipMoles(const ReactionModuleConfig &conf,
                                            const std::vector<size_t>& cellToDomain)
{

    //    std::stringstream ss;
    //    for(size_t i = 0; i < cellToDomain.size(); ++i)
    //    {
    //        ss << MPIManager::getInstance().getRank() << " cellToDomain: "
    //           << i << "\t" << cellToDomain[i] << std::endl;
    //    }
    //    std::cout << ss.str();


    const CellIndexCorrection& indexCorr = conf.indexCorr;
    const std::vector<Domain>& domains = conf.ippConfig.domains;

    const std::vector<std::string> phaseNames = m_simData->getPhaseNames();
    const size_t nPhases = phaseNames.size();

    std::map<std::string, size_t> phaseNameToPhaseIndex;
    for(size_t iPhase = 0; iPhase < nPhases; ++iPhase)
    {
        const std::string& phaseName = phaseNames[iPhase];
        assert(phaseNameToPhaseIndex.find(phaseName) == phaseNameToPhaseIndex.end());
        phaseNameToPhaseIndex[phaseName] = iPhase;
    }

    MPI_CHECK_SYNC;


    const std::vector<double>& porosity = m_simData->getPorosity();
    const size_t nCellsLocal = porosity.size();

    std::vector<double>& precipField = m_simData->getPrecipField();
    for(size_t iCellLocal = 0; iCellLocal < nCellsLocal; ++iCellLocal)
    {
        // init with zero
        for(size_t iPhase = 0; iPhase < nPhases; ++iPhase)
        {
            const size_t precipIndex = IndexHelper::toPhreeqcRawIndex(nCellsLocal, iCellLocal, iPhase);
            assert(precipField.size() > precipIndex);
            precipField[precipIndex] = 0.0;
        }

        const size_t iCellGlobal = indexCorr.convertToGlobal(iCellLocal);
        assert(cellToDomain.size() > iCellGlobal);

        const size_t iDomain = cellToDomain[iCellGlobal];

        // leave not initialized cell to zero precip mole
        if(iDomain != (size_t)-1)
        {
            assert(domains.size() > iDomain);
            const Domain& domain = domains[iDomain];
            const Composition& comp = *domain.composition;

            // check for unsupported solution definition:
            // no concentration in solid cell is allowed
            const double& poros = porosity[iCellLocal];
            if(poros < conf.ippConfig.porosLow)
            {
                for(const Composition::ElementConcentration& elem : comp.elemConcVec)
                {
                    if(elem.conc != 0.0)
                    {
                        const std::string errMsg =
                                "at least one concentration in solid composition ("
                                + comp.name + ") is not zero while porosity is lower than porosLow: "
                                + std::to_string(poros);
                        throw std::runtime_error(errMsg);
                    }
                }
            }

            for(const Composition::PhaseDefinition& phaseDef : comp.phases)
            {
                const std::string& phaseName = phaseDef.name;
                const size_t iPhase = phaseNameToPhaseIndex.at(phaseName);
                const size_t precipIndex = IndexHelper::toPhreeqcRawIndex(nCellsLocal, iCellLocal, iPhase);

                assert(precipField.size() > precipIndex);
                precipField[precipIndex] = phaseDef.amount;
            }
        }

    }

    // m_oldPrecipAmount = precipField;
}

void PhreeqcReactionModule::initCells(const size_t nCellsGlobal,
                                      const ReactionModuleConfig &conf,
                                      const PhreeqcGlobalInitData& initData)
{
    MPI_CHECK_SYNC;

    // could be local but for init it does not matter at lot...
    std::vector<size_t> cellToDomain;

    std::vector<int> grid2chem, ic;
    PhreeqcRMHelper::findMappingAndFeatures(conf, initData, nCellsGlobal,
                                            grid2chem, ic, cellToDomain);
    MPI_CHECK_SYNC;

    phreeqcMPIDispatch(m_phreeqc, [&]() { this->initMapping(grid2chem, ic); });

    MPI_CHECK_SYNC;

    this->initPrecipMoles(conf, cellToDomain);

    MPI_CHECK_SYNC;


    phreeqcMPIDispatch(m_phreeqc, [&]()
    {
        m_phreeqc.SetTimeStep(0.0);
        m_phreeqc.SetTime(0.0);
        std::cout << "init RM run" << std::endl;
    });


    this->runCells();

    MPI_CHECK_SYNC;


    // run cells again? why? due to new timestep?
    // check on kinetic stuff if must be reenabled again
    //    phreeqcMPIDispatch(m_phreeqc, [&]()
    //    {
    //        const double dt = m_simData->getDt();
    //        m_phreeqc.SetTimeStep(dt);
    //        m_phreeqc.RunCells();
    //    });

    MPI_CHECK_SYNC;


}

void PhreeqcReactionModule::pushInitData(PhreeqcGlobalInitData& initData)
{    
    MPI_CHECK_SYNC;

    InitPorosity::syncData(initData.tmpGlobalPoros,
                           initData.tmpGlobalPorosCapillary,
                           initData.tmpGlobalDiff,
                           initData.tmpGlobalInertVolFrac,
                           *m_simData);

    MpiTools::syncFromRoot(m_comm, m_reacData.phaseMolarVolumes);

    m_oldPorosPhysical = m_simData->getPorosity();

    // init children calc as well
    if(MPIManager::getInstance().isMainProc() == false)
    {
        m_reacData.porosCalc.setCompNames(&m_simData->getCompNames());
        m_reacData.porosCalc.setPhaseInfos(&m_simData->getPhaseNameToInfos());
        m_reacData.porosCalc.setPhaseNames(&m_simData->getPhaseNames());
    }

}

ConfigurePhreeqcCellsActivity* PhreeqcReactionModule::getConfigCellsActivity()
{
    return m_configActivity.get();
}

AbstractGeometrySync *PhreeqcReactionModule::getGeomSync()
{
    return &(*m_geomSync);
}

void PhreeqcReactionModule::saveCheckpoint()
{
    MPI_CHECK_SYNC;

    const std::string iterationName = CreateFileName::create(m_simData->getIteration());
    const boost::filesystem::path iterationPath = m_dumpDir / iterationName;
    FileSystemUtils::createIfNotExisting(iterationPath);


    savePorosity(iterationPath);
    this->saveOldPorosity(iterationPath);


    phreeqcMPIDispatch(m_phreeqc, [&]() { this->saveDump(iterationPath); });

    this->saveNodeInfos(iterationPath);

    this->saveAuxData(iterationPath);

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::savePorosity(const boost::filesystem::path iterationPath) const
{
    MPI_CHECK_SYNC;

    const std::vector<double>& localPoros = m_simData->getPorosity();

    std::vector<double> globalPoros;
    const FieldDecomposition* decomp = m_simData->getDecomp();
    MultiToSerialDataSync<double> sync(globalPoros, *decomp);
    sync.pullToRoot(localPoros);

    if(MPIManager::getInstance().isMainProc())
    {
        const boost::filesystem::path porosPath = iterationPath / s_currPorosFilename;
        std::ofstream file(porosPath.string());
        boost::archive::binary_oarchive oa(file);
        oa << globalPoros;
    }

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::loadPorosity(const boost::filesystem::path& iterationPath)
{
    if(MPIManager::getInstance().isMainProc())
    {
        const boost::filesystem::path porosPath = iterationPath / s_currPorosFilename;
        IPPCheck::assertCheck(boost::filesystem::exists(porosPath), "porosity file not found: "
                              + porosPath.string());

        std::ifstream file(porosPath.string());
        boost::archive::binary_iarchive ia(file);

        std::vector<double> globalPoros;
        ia >> globalPoros;

        m_phreeqc.SetPorosity(globalPoros);
        m_phreeqc.MpiWorkerBreak();
    }
    else
    {
        m_phreeqc.MpiWorker();
    }

    std::vector<double>& poros = m_simData->getPorosity();
    poros = m_phreeqc.getPorosityLocal();
}

void PhreeqcReactionModule::saveOldPorosity(const boost::filesystem::path iterationPath) const
{
    MPI_CHECK_SYNC;

    std::vector<double> globalPoros;
    const FieldDecomposition* decomp = m_simData->getDecomp();
    MultiToSerialDataSync<double> sync(globalPoros, *decomp);
    sync.pullToRoot(m_oldPorosPhysical);

    if(MPIManager::getInstance().isMainProc())
    {
        const boost::filesystem::path porosPath = iterationPath / s_oldPorosFilename;
        std::ofstream file(porosPath.string());
        boost::archive::binary_oarchive oa(file);
        oa << globalPoros;
    }

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::loadOldPorosity(const boost::filesystem::path& iterationPath)
{
    MPI_CHECK_SYNC;

    std::vector<double> globalPoros;

    if(MPIManager::getInstance().isMainProc())
    {
        const boost::filesystem::path porosPath = iterationPath / s_oldPorosFilename;
        IPPCheck::assertCheck(boost::filesystem::exists(porosPath), "porosity file not found: "
                              + porosPath.string());

        std::ifstream file(porosPath.string());
        boost::archive::binary_iarchive ia(file);

        ia >> globalPoros;

        std::cout << "read old porosity file. array size is: " << globalPoros.size() << std::endl;
    }

    MPI_CHECK_SYNC;

    const FieldDecomposition* decomp = m_simData->getDecomp();
    SerialToMultiDataSync<double> sync(m_oldPorosPhysical, *decomp);
    sync.pushFromRoot(globalPoros);

    MPI_CHECK_SYNC;

}


void PhreeqcReactionModule::saveDump(const boost::filesystem::path& iterationPath)
{
    assert(MPIManager::getInstance().isMainProc());

    const boost::filesystem::path outPath = iterationPath / s_dumpFilename;

    const std::string outFile = outPath.string();
    PhreeqcDump::dump(m_phreeqc, outFile);


    // write decomposition
    const FieldDecomposition* decomp = m_simData->getDecomp();

    const boost::filesystem::path decompPath = iterationPath / s_decompFilename;
    std::ofstream file(decompPath.string());
    boost::archive::text_oarchive oa(file);

    FieldDecomposition copy = *decomp;
    oa << copy;

}

void PhreeqcReactionModule::loadDump(const boost::filesystem::path& iterationPath)
{
    MPI_CHECK_SYNC;

    double vm, rss;

    ProcessMemoryUsage::getMemUsage(vm, rss);
    std::cout << MPIManager::getInstance().getRank()
              << ": mem usage (init): " << vm << " kb" << std::endl;

    {
        std::istream* dump = nullptr;


        const boost::filesystem::path outPath = iterationPath / s_dumpFilename;
        IPPCheck::assertCheck(boost::filesystem::exists(outPath), "Dump file is not existing: "
                              + outPath.string());
        std::ifstream file(outPath.string(), std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_istream inStrm;
        inStrm.push(boost::iostreams::gzip_decompressor());
        inStrm.push(file);


        std::stringstream correctedString;
        boost::iostreams::filtering_istream correctedInStream;
        correctedInStream.push(boost::iostreams::gzip_decompressor());
        correctedInStream.push(correctedString);

        ProcessMemoryUsage::getMemUsage(vm, rss);
        std::cout << MPIManager::getInstance().getRank()
                  << ": mem usage (opening file): " << vm << " kb" << std::endl;

        {
            // MPI decompositions could have changed
            // load old decomposition to which checkpoint was created
            // and adapt the data to the current decomposition
            const boost::filesystem::path decompPath = iterationPath / s_decompFilename;
            const std::string& decompFilename = decompPath.string();
            pcout << "loading decomp: " << decompFilename << std::endl;

            std::ifstream ifs(decompFilename);
            if(ifs.is_open() == false)
            {
                throw std::runtime_error("could not open decomp file: " + decompPath.string());
            }
            boost::archive::text_iarchive ia(ifs);
            FieldDecomposition oldDecomp;
            ia >> oldDecomp;


            const FieldDecomposition* decomp = m_simData->getDecomp();

            bool isIdentical = *decomp == oldDecomp;
            boost::mpi::broadcast(m_comm, isIdentical, 0);

            if(isIdentical)
            {
                pcout << "decomposition (luckily) has not changed: " << std::endl
                      << oldDecomp.getDomains().size() << " -> "
                      << decomp->getDomains().size() << std::endl;
                dump = &inStrm;
            }
            else
            {
                pcout << "decomposition changed: " << std::endl
                      << oldDecomp.getDomains().size() << " -> "
                      << decomp->getDomains().size() << std::endl
                      << "\t-> fixing cell IDs" << std::endl;

                boost::iostreams::filtering_ostream correctedOutStream;
                correctedOutStream.push(boost::iostreams::gzip_compressor());
                correctedOutStream.push(correctedString);

                FixCellsIDsForNewFieldDecomposition::fix(oldDecomp, *decomp,
                                                         inStrm, correctedOutStream);

                dump = &correctedInStream;

                // shall we write the updated dump with the corrected decomp to a file?
            }

        }


        ProcessMemoryUsage::getMemUsage(vm, rss);
        std::cout << MPIManager::getInstance().getRank()
                  << ": mem usage (fixed decomp): " << vm << " kb" << std::endl;


        // maybe compress dump again and decompress on the other side
        // boost::mpi::broadcast(m_comm, corrCompressed, 0);

        pcout << "running dump in PhreeqC" << std::endl;
        const IRM_RESULT result = m_phreeqc.runStringLocal(true, false, false, *dump);
        IPPCheck::assertCheck(result == IRM_OK, "runString error on loading dump");
        pcout << "finished to run the dump in PhreeqC" << std::endl;

        ProcessMemoryUsage::getMemUsage(vm, rss);
        std::cout << MPIManager::getInstance().getRank()
                  << ": mem usage (finished run dump): " << vm << " kb" << std::endl;

    }

    ProcessMemoryUsage::getMemUsage(vm, rss);
    std::cout << MPIManager::getInstance().getRank()
              << ": mem usage (finished loading dump): " << vm << " kb" << std::endl;

    MPI_CHECK_SYNC;

}

void PhreeqcReactionModule::saveNodeInfos(const boost::filesystem::path& iterationPath) const
{    
    MPI_CHECK_SYNC;

    // collect all node infos from all processes
    const FieldDecomposition* decomp = m_simData->getDecomp();
    NodeInfos globalInfos;

    globalInfos.nInterfaceNodesGlobal = m_nodeInfos.nInterfaceNodesGlobal;

    {
        MultiToSerialDataSync<char> sync(globalInfos.interfaceCells, *decomp);
        sync.pullToRoot(m_nodeInfos.interfaceCells);
    }

    MPI_CHECK_SYNC;

    {
        MultiToSerialDataSync<char> sync(globalInfos.nonPermNonInterfaceCells, *decomp);
        sync.pullToRoot(m_nodeInfos.nonPermNonInterfaceCells);
    }

    MPI_CHECK_SYNC;

    if(MPIManager::getInstance().isMainProc())
    {
        const boost::filesystem::path path = iterationPath / s_nodeInfosFilename;
        std::ofstream file(path.string());

        if(file.is_open() == false)
        {
            throw std::runtime_error("could not open node file: " + path.string());
        }

        boost::archive::text_oarchive oa(file);
        oa << globalInfos;
    }

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::loadNodeInfos(const boost::filesystem::path& iterationPath)
{
    MPI_CHECK_SYNC;

    NodeInfos globalInfos;

    if(MPIManager::getInstance().isMainProc())
    {
        const boost::filesystem::path path = iterationPath / s_nodeInfosFilename;
        IPPCheck::assertCheck(boost::filesystem::exists(path), "node file not found: "
                              + path.string());

        std::ifstream file(path.string());

        if(file.is_open() == false)
        {
            throw std::runtime_error("could not open node file: " + path.string());
        }

        boost::archive::text_iarchive ia(file);
        ia >> globalInfos;
    }


    MPI_CHECK_SYNC;

    const FieldDecomposition* decomp = m_simData->getDecomp();

    {
        SerialToMultiDataSync<char> sync(m_nodeInfos.interfaceCells, *decomp);
        sync.pushFromRoot(globalInfos.interfaceCells);
    }

    MPI_CHECK_SYNC;

    {
        SerialToMultiDataSync<char> sync(m_nodeInfos.nonPermNonInterfaceCells, *decomp);
        sync.pushFromRoot(globalInfos.nonPermNonInterfaceCells);
    }

    m_nodeInfos.nInterfaceNodesGlobal = globalInfos.nInterfaceNodesGlobal;
    boost::mpi::broadcast(m_comm, m_nodeInfos.nInterfaceNodesGlobal, 0);

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::saveAuxData(const boost::filesystem::path &iterationPath) const
{
    MPI_CHECK_SYNC;

    // collect all node infos from all processes

    std::vector<double> globalEqConc;
    {
        FieldDecomposition decomp = *m_simData->getDecomp();
        const size_t nComps = m_simData->getCompNames().size();
        decomp.correctByTensorDim(nComps);

        const std::vector<double>& lastEqConc = m_configActivity->getLastPhreeqcEqConc();
        MultiToSerialDataSync<double> sync(globalEqConc, decomp, nComps);
        sync.pullToRoot(lastEqConc);
    }

    MPI_CHECK_SYNC;

    std::vector<double> globalPrecip;
    {
        FieldDecomposition decomp = *m_simData->getDecomp();
        const size_t nPhases = m_simData->getPhaseNames().size();
        decomp.correctByTensorDim(nPhases);

        const std::vector<double>& precip = m_simData->getPrecipField();
        MultiToSerialDataSync<double> sync(globalPrecip, decomp, nPhases);
        sync.pullToRoot(precip);
    }

    MPI_CHECK_SYNC;

    if(MPIManager::getInstance().isMainProc())
    {
        const boost::filesystem::path path = iterationPath / s_reacAuxFilename;
        std::ofstream file(path.string());

        if(file.is_open() == false)
        {
            throw std::runtime_error("could not open last equi file: " + path.string());
        }

        boost::archive::binary_oarchive oa(file);
        // oa << globalPreReacConc;
        oa << globalEqConc;
        oa << globalPrecip;
    }

    MPI_CHECK_SYNC;
}
void PhreeqcReactionModule::loadAuxData(const boost::filesystem::path& iterationPath)
{
    MPI_CHECK_SYNC;

    // std::vector<double> globalPreReacConc;
    std::vector<double> globalEqConc;
    std::vector<double> globalPrecip;

    if(MPIManager::getInstance().isMainProc())
    {
        const boost::filesystem::path path = iterationPath / s_reacAuxFilename;
        IPPCheck::assertCheck(boost::filesystem::exists(path), "last equi file not found: "
                              + path.string());

        std::ifstream file(path.string());

        if(file.is_open() == false)
        {
            throw std::runtime_error("could not open last equi file: " + path.string());
        }

        boost::archive::binary_iarchive ia(file);
        // ia >> globalPreReacConc;
        ia >> globalEqConc;
        ia >> globalPrecip;
    }

    MPI_CHECK_SYNC;


    {
        pcout << "sync lastEqConc: " << globalEqConc.size() << std::endl;
        FieldDecomposition decomp = *m_simData->getDecomp();
        const size_t nComps = m_simData->getCompNames().size();
        decomp.correctByTensorDim(nComps);

        std::vector<double>& lastEqConc = m_configActivity->getLastPhreeqcEqConc();
        // sync module does resize internally
        // assert(lastEqConc.size() == decomp.getLocalNumberOfCells());
        SerialToMultiDataSync<double> sync(lastEqConc, decomp, nComps);
        sync.pushFromRoot(globalEqConc);
    }

    MPI_CHECK_SYNC;

    {
        pcout << "sync precipField: " << globalPrecip.size() << std::endl;

        FieldDecomposition decomp = *m_simData->getDecomp();
        const size_t nPhases = m_simData->getPhaseNames().size();
        decomp.correctByTensorDim(nPhases);

        std::vector<double>& precip = m_simData->getPrecipField();
        SerialToMultiDataSync<double> sync(precip, decomp, nPhases);
        sync.pushFromRoot(globalPrecip);
    }

    MPI_CHECK_SYNC;

}

void PhreeqcReactionModule::loadCheckpoint(const size_t iteration)
{
    MPI_CHECK_SYNC;

    pcout << "loading reaction module checkpoint: " << iteration << std::endl;


    const std::string iterationName = CreateFileName::create(iteration);

    const boost::filesystem::path iterationPath = m_dumpDir / iterationName;
    IPPCheck::assertCheck(boost::filesystem::exists(iterationPath), "checkpoint folder not found: "
                          + iterationPath.string());


    // porosity must be set properly before dump is loaded
    this->loadPorosity(iterationPath);

    this->loadOldPorosity(iterationPath);


    this->loadDump(iterationPath);


    std::vector<double>& conc = m_simData->getPostReacConc();
    // this->getUpdatedConcentrations(m_simData->getPostReacConc());
    std::fill(conc.begin(), conc.end(), -1.0);


    this->loadNodeInfos(iterationPath);

    this->loadAuxData(iterationPath);


    pcout << "loading reaction module checkpoint fully completed" << std::endl;

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::errorDump()
{
    MPI_CHECK_SYNC;

    phreeqcMPIDispatch(m_phreeqc, [&]() {

        assert(MPIManager::getInstance().isMainProc());
        const size_t iteration = m_simData->getIteration();
        const boost::filesystem::path outPath = m_chemDir / ("mainrun_err_" + std::to_string(iteration));
        const std::string outFile = outPath.string();
        PhreeqcDump::dump(m_phreeqc, outFile);


    });


    MPI_CHECK_SYNC;
}


void PhreeqcReactionModule::initConcentrations()
{
    MPI_CHECK_SYNC;

    std::vector<double>& conc = m_simData->getPostReacConc();
    assert(conc.empty());
    m_phreeqc.getConcentrationsLocal(conc);


    // copy the very init concentrations since the lowPoros cells will maybe
    // enabled later. we need to set a proper porosity corrected concentrations there
    // see also updateComponentConcentrations();
    m_initialData.conc = conc;


    std::vector<double>& lastPhreeqcEqConc = m_configActivity->getLastPhreeqcEqConc();
    lastPhreeqcEqConc = conc;


    // set disabled cells conc to zero
    const size_t nCells = m_simData->getDecomp()->getLocalNumberOfCells();
    const size_t nComps = m_simData->getCompNames().size();
    assert(conc.size() == nCells * nComps);
    assert(m_enabledCells.size() == nCells);

    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {
        if(m_enabledCells[iCell] == false)
        {
            for(size_t iComp = 0; iComp < nComps; ++iComp)
            {
                const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iComp);
                conc[index] = 0.0;
            }
        }
    }


    m_simData->getPostTransportConc().resize(conc.size(), -1.0);

    // init with zero because when retrieving conc disabled cells should be set to zero
    m_preReacConc.resize(conc.size(), 0.0);

    m_transportDiff.resize(conc.size(), 0.0);

}

template<typename Iterator>
static void checkPhaseVol(const std::string& phaseName, const double& molVol,
                          const Iterator& begin, const Iterator& end)
{
    if(molVol <= 0.0)
    {
        std::vector<Iterator> bla;
        for(Iterator it = begin; it != end; ++it)
        {
            if(*it > 0.0)
            {
                bla.push_back(it);
            }
        }

        if(bla.empty() == false)
        {

            std::string errMsg = "WARNING: non-zero content for phase. missing molar volume for: "
                    + phaseName + " at cellIDs: ";

            for(size_t i = 0; i < bla.size(); ++i)
            {
                const std::size_t cellID = std::distance(begin, bla[i]);
                errMsg += std::to_string(cellID) + ", ";
            }
            throw std::runtime_error(errMsg);
        }

    }

}

void PhreeqcReactionModule::getUpdatedConcentrations(std::vector<double>& postReacConc)
{    
    std::vector<double> newConc;
    m_phreeqc.getConcentrationsLocal(newConc);


    const size_t nCells = m_enabledCells.size();

    // converting charge equivalents to per liter basis?
    // TODO: use which porosity? physical or the one used for calc chemistry
    // they are different in case of hetero chemistry calc
    //    const std::vector<double>& poros = m_simData->getPorosityPhysical();
    //    for(size_t iCell = 0; iCell < nCells; ++iCell)
    //    {
    //        const size_t iComp = 3;
    //        const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iComp);

    //        assert(newConc.size() > index);
    //        double& conc = newConc[index];
    //        conc /= poros[iCell];
    //    }


    const size_t nComps = m_phreeqc.GetComponentCount();
    std::vector<double>& lastPhreeqcEqConc = m_configActivity->getLastPhreeqcEqConc();
    // std::vector<double>& concDiffCummulative = m_configActivity->getConcDiffCummulative();


    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {
        for(size_t iComp = 0; iComp < nComps; ++iComp)
        {
            const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iComp);

            assert(postReacConc.size() > index);
            double& conc = postReacConc[index];

            assert(lastPhreeqcEqConc.size() > index);
            double& eqConc = lastPhreeqcEqConc[index];

            // m_currDissolveOnlyCells seems to be unused. was formerly used to fix issue with below zero porosity
            // if(m_enabledCells[iCell] && m_currDissolveOnlyCells[iCell] == false)
            if(m_enabledCells[iCell])
            {
                assert(newConc.size() > index);
                const double& val = newConc[index];
                assert(iComp == 3 || val >= 0.0);

                conc = val;
                eqConc = val;
            }
            else
            {                
                // postReacConc is currently at the status of preTransConc of previous iteration
                // we need to update disabled cells properly with post transport conc of current iteration
                // postReacConc was set the post transport conc previously
                conc = m_preReacConc[index];

                // double& diffCumm = concDiffCummulative[index];
                // m_disableCellsFunc->accumulateDiff(eqConc, conc, diffCumm);
            }
        }
    }


}

void PhreeqcReactionModule::setTime()
{
    MPI_CHECK_SYNC;

    const double time = 0.0; //m_simData.getIteration();
    m_phreeqc.setTimeLocal(time);

}

void PhreeqcReactionModule::treatNegativePorosity()
{
    MPI_CHECK_SYNC;

    std::vector<double>& porosityPhysical = m_simData->getPorosity();
    const size_t nCells = porosityPhysical.size();

    std::vector<CellTotalsDiff>& diffPerCell = m_simData->additionalSources;
    diffPerCell.clear();

    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {

        // if(m_currDissolveOnlyCells[iCell] == false)
        {

            const double& porosCurr = porosityPhysical[iCell];

            // const double minPoros = 1.0E-6;

            if(porosCurr < 0.0)
            {
                const double& porosOld = m_oldPorosPhysical[iCell];

                std::cout << "found negative porosity in local cell: " << iCell << "\n\t-> "
                          << "changed from " << porosOld << " to " << porosCurr << std::endl;


                const std::vector<double>& porosityChemVec = m_phreeqc.getPorosityLocal();

                assert(porosityChemVec.size() > iCell);
                const double& porosChem = porosityChemVec[iCell];


                const size_t nComps = m_phreeqc.GetComponentCount();

                std::stringstream ss;
                ss << "-------------------------------------\n";
                ss << MPIManager::getInstance().getRank() << ":\n";

                std::vector<double> totalDiff(nComps);
                for(size_t iComp = 0; iComp < nComps; ++iComp)
                {
                    const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iComp);

                    assert(m_transportDiff.size() > index);
                    totalDiff[iComp] = m_transportDiff[index];

                    ss << iCell << "  " << iComp << " diff: " << totalDiff[iComp] << std::endl;
                }

                ss << "porosChem: " << porosChem << std::endl;

                ss << "-------------------------------------\n";
                std::cout << ss.str();

                FindZeroPorosityRoot::SetTotalsAndRun setAndRun(m_phreeqc, totalDiff,
                                                                porosChem, iCell);

                const RetrievePhaseAmounts retrieveAmounts(m_phreeqc,
                                                           m_enabledCells);

                const PorosityCalc porosCalc(nCells, m_reacData.phaseMolarVolumes,
                                             m_reacData.isPermeablePhase,
                                             m_simData->getInertVolFrac());

                const FindZeroPorosityRoot::Result result =
                        FindZeroPorosityRoot::findRoot(setAndRun, retrieveAmounts, porosCalc,
                                                       porosOld, porosCurr, iCell);
                const double& t = result.t;

                // saving amount supposed to be distributed to neighbors
                diffPerCell.push_back(CellTotalsDiff());
                CellTotalsDiff& data = diffPerCell.back();
                data.iCell = iCell;
                std::vector<double>& diffTotals = data.totalsDiff;
                diffTotals.resize(nComps);

                for(size_t iComp = 0; iComp < nComps; ++iComp)
                {
                    const double& diff = totalDiff[iComp];
                    diffTotals[iComp] = diff * (1.0 - t);
                }

                std::cout << porosCurr << " -> " << result.porosity << std::endl;
            }
        }
    }




    MPI_CHECK_SYNC;

    if(diffPerCell.empty() == false)
    {
        this->collectDataAll();
    }

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::runHomo()
{
    m_configActivity->getInterfaceSwitch().setIsInterfaceRun(false);


    this->configCellsActivity();
    this->runCells();

#ifdef IPP_DEBUG_ENABLED
    m_simData->debugData.enabledCellsSteps[0] = m_enabledCells;
#endif
}

void PhreeqcReactionModule::runHetero()
{        
    auto start = std::chrono::steady_clock::now();

    m_configActivity->getInterfaceSwitch().setIsInterfaceRun(true);


    const std::vector<double>& postReacConc = m_simData->getPostReacConc();
    const std::vector<double>& poros = m_simData->getPorosity();

    m_preReacConc = postReacConc;
    std::vector<double> avgPoros(poros.size(), -1.0);
    m_nonLocalOperations.interfaceProperties->execute(postReacConc, poros, m_preReacConc, avgPoros);

    m_phreeqc.setPorosityLocal(avgPoros);
    m_phreeqc.setConcentrationsLocal(m_preReacConc);


    const std::vector<double>& oldEqConc = m_configActivity->getLastPhreeqcEqConc();
    const size_t nCells = avgPoros.size();

    for(size_t i = 0; i < m_preReacConc.size(); ++i)
    {
        const size_t iCell = i % nCells;
        const double& porosity = avgPoros[iCell];
        const double& preReacConc = m_preReacConc[i];
        const double& oldEq = oldEqConc[i];
        const double diff = preReacConc - oldEq;
        m_transportDiff[i] = diff * porosity;
    }


    if(m_verbose)
    {
        const auto duration = CalcDuration::calc(start);
        pcout << "\t-> calc interface properties took:\t"
              << duration.count() << " ms" << std::endl;
    }

    this->configCellsActivity();
    this->runCells();


#ifdef IPP_DEBUG_ENABLED
    m_simData->debugData.enabledCellsSteps[1] = m_enabledCells;
#endif

}

void PhreeqcReactionModule::run()
{
    MPI_CHECK_SYNC;

    this->swapToOld();

    std::vector<SimplePorosityInfosPtr> homoPorosInfos;

    {
        auto start = std::chrono::steady_clock::now();
        this->runHomo();

        PrecipData homoReacData;
        this->collectData(homoReacData);

        this->calcCapillaryPorosity(homoReacData.localPrecipMol,
                                    homoReacData.solidFractions,
                                    homoPorosInfos);

        this->treatNegativePorosity();

        std::vector<double>& conc = m_simData->getPostReacConc();
        this->correctByPorosChange(conc);

        if(m_verbose)
        {
            const auto duration = CalcDuration::calc(start);
            pcout << "\t-> runHomo took:\t"
                  << duration.count() << " ms" << std::endl;
        }
    }

    if(m_nodeInfos.nInterfaceNodesGlobal > 0)
    {
        auto start = std::chrono::steady_clock::now();
        this->runHetero();

        PrecipData tmpData;
        this->collectData(tmpData);

        std::vector<SimplePorosityInfosPtr> heteroPorosInfos;
        this->calcCapillaryPorosity(tmpData.localPrecipMol,
                                    tmpData.solidFractions,
                                    heteroPorosInfos);

        calcNewDiffCoeffs(m_reacData.diffCalc,
                          heteroPorosInfos,
                          m_simData->getDiffusionCoefs(),
                          m_verbose);

        this->treatNegativePorosity();

        if(m_verbose)
        {
            const auto duration = CalcDuration::calc(start);
            pcout << "\t-> runHetero took:\t"
                  << duration.count() << " ms" << std::endl;
        }
    }
    else
    {
        // if hetero reactions are disabled calc new diff coefs from homo poros infos
        calcNewDiffCoeffs(m_reacData.diffCalc, homoPorosInfos, m_simData->getDiffusionCoefs(), m_verbose);
    }


    if(m_phreeqcData.dissolveOnlyFunc->needsNeighInfos())
    {
        this->collectSteadyState();
    }

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::retrieveConcentrations()
{
    auto start = std::chrono::steady_clock::now();

    std::vector<double>& postReacConc = m_simData->getPostReacConc();
    this->getUpdatedConcentrations(postReacConc);

    if(m_verbose)
    {
        const auto duration = CalcDuration::calc(start);
        pcout << "\t-> retrieve concentrations from PhreeqC took:\t"
              << duration.count() << " ms" << std::endl;
    }
}



void PhreeqcReactionModule::calcCapillaryPorosity(const std::vector<double>& precipMol,
                                                  const std::vector<double>& solidFractions,
                                                  std::vector<SimplePorosityInfosPtr>& porosInfos)
{
    const PhaseNameToInfos& phaseInfos = m_simData->getPhaseNameToInfos();
    const std::vector<double>& inertVolFrac = m_simData->getInertVolFrac();

    PhreeqcCalcPorosity::calc(phaseInfos, solidFractions, precipMol,
                              inertVolFrac, m_reacData.porosCalc, porosInfos);

    const size_t nPhases = m_phreeqcData.phaseToNuclPhaseID.size();
    const size_t nCells = m_simData->getDecomp()->getLocalNumberOfCells();

    for(size_t iPhase = 0; iPhase < nPhases; ++iPhase)
    {
        const size_t iNuclPhaseId = m_phreeqcData.phaseToNuclPhaseID[iPhase];
        if(iNuclPhaseId != (size_t)-1)
        {
            std::vector<double>& volumePerCell = m_simData->nucleationPhasesVolumeFractions[iNuclPhaseId];
            for(size_t iCell = 0; iCell < nCells; ++iCell)
            {
                const size_t srcIndex = IndexHelper::toPhreeqcTransposedIndex(nPhases, iCell, iPhase);
                volumePerCell[iCell] = solidFractions[srcIndex];
            }
        }
    }


    std::vector<double>& capillaryPorosity = m_simData->getCapillaryPorosity();
    assert(capillaryPorosity.size() == porosInfos.size());
    for(size_t iCell = 0; iCell < capillaryPorosity.size(); ++iCell)
    {
        const SimplePorosityInfosPtr& info = porosInfos[iCell];
        capillaryPorosity[iCell] = info->porosityCapillary;
    }



}

void PhreeqcReactionModule::collectDataAll()
{
    PrecipData tmpData;
    this->collectData(tmpData);

    std::vector<SimplePorosityInfosPtr> porosInfos;
    this->calcCapillaryPorosity(tmpData.localPrecipMol,
                                tmpData.solidFractions,
                                porosInfos);
    calcNewDiffCoeffs(m_reacData.diffCalc,
                      porosInfos,
                      m_simData->getDiffusionCoefs(),
                      m_verbose);

    if(m_phreeqcData.dissolveOnlyFunc->needsNeighInfos())
    {
        this->collectSteadyState();
    }
}

void PhreeqcReactionModule::collectSteadyState()
{
    const FieldDecomposition& decomp = *m_simData->getDecomp();
    const size_t nCells = decomp.getLocalNumberOfCells();
    const size_t nPhases = m_simData->getPhaseNames().size();

    const std::vector<double>& precipMol = m_simData->getPrecipField();

    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {
        for(size_t iPhase = 0; iPhase < nPhases; ++iPhase)
        {
            const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iPhase);

            assert(m_oldPrecipAmount.size() > index);
            const double& oldMol = m_oldPrecipAmount[index];

            assert(precipMol.size() > index);
            const double& newMol = precipMol[index];

            const double diff = newMol - oldMol;

            const double currSI = m_currTargetSImap[index];

            if(diff <= 0 && currSI > 0.0)
            {
                m_simData->solidFlags[iCell] = IPPConstants::s_isSteadyFlag;
            }
            else
            {
                m_simData->solidFlags[iCell] = IPPConstants::s_isNotSteadyFlag;
                break;
            }
        }

        // std::cout << iCell << "\t" << m_simData->solidFlags[iCell] << std::endl;
    }
}

void PhreeqcReactionModule::collectData(PrecipData& tmpData)
{
    // IMPORTANT: concentrations must be retrieved before porosity inside of phreeqc gets updated
    this->retrieveConcentrations();

    {
        auto start = std::chrono::steady_clock::now();

        const RetrievePhaseAmounts retrievePhases(m_phreeqc, m_enabledCells);
        std::vector<double>& precipMol = m_simData->getPrecipField();
        m_oldPrecipAmount = precipMol;
        retrievePhases.retrieve(precipMol);


        // TODO: optimize out if last amount is not needed for any optimization
        std::vector<double>& lastAmounts = m_configActivity->getLastPhreeqcPhaseAmount();
        lastAmounts = precipMol;

        const FieldDecomposition* decomp = m_simData->getDecomp();
        const size_t nCellsLocal = decomp->getLocalNumberOfCells();
        const size_t nPhases = m_simData->getPhaseNames().size();
        GeometryTools::transpose(lastAmounts, nCellsLocal, nPhases);

        if(m_verbose)
        {
            const auto duration = CalcDuration::calc(start);
            pcout << "\t-> retrieve phases moles from PhreeqC took:\t"
                  << duration.count() << " ms" << std::endl;
        }

    }

    if(m_phreeqcData.dissolveOnlyFunc->needsSaturationIndices())
    {
        const auto start = std::chrono::steady_clock::now();

        std::vector<double>& satIndices = m_simData->getSaturationIndices();
        this->retrieveSaturationIndices(satIndices);

        if(m_verbose)
        {
            const auto duration = CalcDuration::calc(start);
            pcout << "\t-> retrieve SI data from PhreeqC took:\t"
                  << duration.count() << " ms" << std::endl;
        }
    }

    // TODO: rename the flag in DO func or find a better check
    if(m_phreeqcData.dissolveOnlyFunc->needsNeighInfos())
    {
        const auto start = std::chrono::steady_clock::now();

        this->retrieveMonomersConc();

        if(m_verbose)
        {
            const auto duration = CalcDuration::calc(start);
            pcout << "\t-> retrieve monomer conc from PhreeqC took:\t"
                  << duration.count() << " ms" << std::endl;
        }
    }

    {
        const auto start = std::chrono::steady_clock::now();

        AuxDataVec& auxData = m_simData->getAuxData();

        RetrieveAuxValues retrieveData(m_phreeqc, m_enabledCells);
        retrieveData.retrieve(auxData);

        if(m_verbose)
        {
            const auto duration = CalcDuration::calc(start);
            pcout << "\t-> retrieve aux data from PhreeqC took:\t"
                  << duration.count() << " ms" << std::endl;
        }

    }

    PorosityHelper::CalcPorosityOutput output =
    {
        tmpData.localPrecipMol,
        tmpData.solidFractions,
        m_simData->getPorosity()
    };

    PorosityHelper::calc(m_reacData, output, m_verbose);

}

void PhreeqcReactionModule::collectPorosityConvergenceData()
{
    MPI_CHECK_SYNC;
    this->collectDataAll();
    std::vector<double>& conc = m_simData->getPostReacConc();
    this->correctByPorosChange(conc);
    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::updateComponentConcentrations()
{
    MPI_CHECK_SYNC;

    m_preReacConc = m_simData->getPostTransportConc();


    {
        // sometimes the subtraction of very small values for "diff" leads to INV fpe...
        // ScopedDisableFloatingPointException fpe;

        const std::vector<double>& postReacConcVec = m_simData->getPostReacConc();
        const std::vector<double>& porosityChem = m_phreeqc.getPorosityLocal();

        const size_t nCells = porosityChem.size();

        assert(m_transportDiff.size() == m_preReacConc.size());

        for(size_t i = 0; i < m_preReacConc.size(); ++i)
        {
            // const size_t iComp = i / nCells;

            const size_t iCell = i % nCells;
            const double& porosity = porosityChem[iCell];
            const double& preReacConc = m_preReacConc[i];
            const double& postReacConc = postReacConcVec[i];

            const double diff = preReacConc - postReacConc;
            const double absDiff = diff * porosity;
            m_transportDiff[i] = absDiff;
        }



        // correct for volume effect on charge balance?
        // TODO: check for right poros. see also function where conc is retrieved
        // TODO: check if volume effect must be taken into account for m_transportDiff and negative poros treatment
        // const std::vector<double>& poros = m_simData->getPorosityPhysical();
        // assert(poros.size() == nCells);

        // for(size_t iCell = 0; iCell < nCells; ++iCell)
        // {
        //     const size_t iComp = 3;
        //     const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iComp);

        //     assert(m_preReacConc.size() > index);
        //     double& conc = m_preReacConc[index];
        //     conc *= poros[iCell];
        // }
    }


    IRM_RESULT status = m_phreeqc.setConcentrationsLocal(m_preReacConc);
    IPPCheck::assertCheck(status == IRM_OK);

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::configCellsActivity()
{
    std::vector<size_t> toEnable, toDisable;

    {
        auto start = std::chrono::steady_clock::now();

        const size_t it = m_simData->getIteration();
        m_configActivity->run(it, m_enabledCells, toEnable, toDisable);

        MPI_CHECK_SYNC;

        if(m_verbose)
        {
            const auto duration = CalcDuration::calc(start);
            pcout << "\t-> configure cells activity in PhreeqC took:\t"
                  << duration.count() << " ms" << std::endl;
        }
    }

    {
        const auto start = std::chrono::steady_clock::now();

        //// DEBUG
#ifdef DEBUG_ENABLED_CELLS
        {
            const size_t nCellsLocal = oldEnabledCells.size();
            const size_t nEnabled = std::count(oldEnabledCells.begin(), oldEnabledCells.end(), true);

            const int rank = MPIManager::getInstance().getRank();
            const std::string filename = "enabled_" + std::to_string(rank) + ".csv";
            const boost::filesystem::path enabledData = m_chemDir / filename;

            const size_t it = m_simData->getIteration();

            std::ios_base::openmode flag = std::ios_base::trunc;
            if(it > 1)
            {
                flag = std::ios_base::app;
            }

            std::ofstream file(enabledData.string(), flag);

            //            if(it == 0)
            //            {
            //                file << "it\tenabled\tdisabled\tfrac\n";
            //            }



            file
                    //                    << it
                    //                    << "\t" << nEnabled
                    //                    << "\t" << nCellsLocal - nEnabled
                    //                    << "\t"
                    << (double)nEnabled / (double)nCellsLocal
                    << std::endl;

            file.close();


            //            std::cout << "new enabled: ";
            //            for(size_t i = 0; i < toEnable.size(); ++i)
            //            {
            //                std::cout << toEnable[i] << " ";
            //            }
            //            std::cout << std::endl;

            //            std::cout << "new disabled: ";
            //            for(size_t i = 0; i < toDisable.size(); ++i)
            //            {
            //                std::cout << toDisable[i] << " ";
            //            }
            //            std::cout << std::endl;


            //            const size_t nComp = m_simData->getCompNames().size();

            //            const std::vector<double>& concDiffCummulative =
            // m_configActivity->getConcDiffCummulative();
            //            assert(postTransConc.size() == concDiffCummulative.size());

            //            const size_t chargeIndex = PhreeqcConstants::s_primaryCompsBegin - 1;


            //// DEBUG
            //            const size_t nCells = oldEnabledCells.size();
            //            for(size_t iCell = 0; iCell < nCells; ++iCell)
            //            {
            //                const size_t iSrc = IndexHelper::toPhreeqcRawIndex(nCells, iCell, 4);
            //                std::cout << iCell
            //                          << "\t" << (double)m_simData->getEnabledCells()[iCell]
            //                             << "\t" << postTransConc[iSrc]
            //                                << "\t+\t" << concDiffCummulative[iSrc]
            //                                   << "\t=\t" << postTransConc[iSrc]
            //                                              + concDiffCummulative[iSrc]
            //                                      << std::endl;

            //                for(size_t iComp = 0; iComp < nComp; ++iComp)
            //                {
            //                    const size_t iSrc = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iComp);
            //                    const double& oldConc = postTransConc[iSrc];
            //                    const double& diff = concDiffCummulative[iSrc];

            //                    assert(iComp == chargeIndex ?
            //                               true : (oldConc + diff) >= 0.0);
            //                }
            //            }
        }
        /////////////
#endif

        m_enableCellsFunc->enable(toEnable);
        m_disableCellsFunc->disable(toDisable);


        const size_t nEnabled = std::count(m_enabledCells.begin(), m_enabledCells.end(), 1);

        if(nEnabled == 0)
        {
            // disable phreeqc instance completly
            m_phreeqc.setGlobalEnabled(false);
        }
        else
        {
            m_phreeqc.setGlobalEnabled(true);
        }


        if(m_verbose)
        {
            const auto duration = CalcDuration::calc(start);
            pcout << "\t-> update conc due to cells activity change took:\t"
                  << duration.count() << " ms" << std::endl;
        }

    }

}



void PhreeqcReactionModule::disableLowPorosCells(const std::vector<double>& porosVec,
                                                 const double& lowPorosThresh)
{
    MPI_CHECK_SYNC;

    assert(porosVec.size() == m_enabledCells.size());


    // phreeqc does not like zero porosity upon init, even for disabled cells
    // this would give convergence errors when enabling the cell later
    std::vector<double> porosity = porosVec;
    for(double& poros : porosity)
    {
        poros = std::max(poros, 1.0E-6);
    }
    m_phreeqc.setPorosityLocal(porosity);


    // copy the very init porosity since the lowPoros cells will maybe
    // enabled later. we need to set a proper porosity corrected concentrations there
    // see also updateComponentConcentrations();
    m_initialData.porosity = porosity;


    // eventually disable those cells
    for(size_t iCell = 0; iCell < porosVec.size(); ++iCell)
    {
        const double& poros = porosVec[iCell];
        if(poros < lowPorosThresh)
        {
            m_phreeqc.disableCell(iCell);
            m_enabledCells[iCell] = false;
        }
    }
}

void PhreeqcReactionModule::initRoot(const size_t nCellsGlobal,
                                     const ReactionModuleConfig &conf,
                                     PhreeqcGlobalInitData& initData)
{
    assert(MPIManager::getInstance().isMainProc());

    std::cout << "initialize PhreeqcRM root" << std::endl;


    const std::set<std::string>& waterLimitedPhases = conf.ippConfig.waterLimitedPhases;
    const std::set<std::string>& solidSolutionPhases = conf.ippConfig.solidSolutionPhases;

    std::vector<std::string>& eqPhaseNames = m_simData->getEqPhaseNames();
    const std::vector<std::string>& phaseNames = m_simData->getPhaseNames();
    for(const std::string& name : phaseNames)
    {
        if(waterLimitedPhases.find(name) == waterLimitedPhases.end() &&
                solidSolutionPhases.find(name) == solidSolutionPhases.end())
        {
            eqPhaseNames.push_back(name);
        }
    }
    std::sort(eqPhaseNames.begin(), eqPhaseNames.end());





    // TODO: seperate function

    // disable rebalance
    m_phreeqc.SetRebalanceFraction(0.0);


    const PhaseNameToInfos& phaseInfos = m_simData->getPhaseNameToInfos();
    InitPorosity::extractPhasesMolarVolumes(phaseInfos, phaseNames,
                                            m_reacData.phaseMolarVolumes);



    IRM_RESULT status;

    // throw exception in case of error
    // m_phreeqc.SetErrorHandlerMode(1);

    // no of throw of exception in case of error
    // return error code
    m_phreeqc.SetErrorHandlerMode(0);


    // TODO: initialize correctly
    m_phreeqc.SetComponentH2O(PhreeqcConstants::s_calcH2OSeperately);
    m_phreeqc.SetSpeciesSaveOn(false);

    const boost::filesystem::path prefixPath = m_chemDir / "mainrun_init";
    const std::string mainRunFilePrefix = prefixPath.string();
    m_phreeqc.SetFilePrefix(mainRunFilePrefix);
    m_phreeqc.OpenFiles();



    m_phreeqc.UseSolutionDensityVolume(false);

    m_phreeqc.SetUnitsSolution(PhreeqcConstants::s_unitSolution);
    m_phreeqc.SetUnitsPPassemblage(PhreeqcConstants::s_unitPPassemblage);
    m_phreeqc.SetUnitsKinetics(PhreeqcConstants::s_unitKinetics);
    m_phreeqc.SetUnitsSSassemblage(PhreeqcConstants::s_unitSSassemblage);
    //    m_phreeqc.SetUnitsExchange(0);
    //    m_phreeqc.SetUnitsSurface(0);
    //    m_phreeqc.SetUnitsGasPhase(0);

    std::vector<double> unitVec(nCellsGlobal, 1.0);

    // const std::vector<double> rv(m_nxyz, 0.1);
    m_phreeqc.SetRepresentativeVolume(unitVec);
    m_phreeqc.SetPorosity(unitVec);
    m_phreeqc.SetSaturation(unitVec);


    const double time_conversion = 1.0;
    status = m_phreeqc.SetTimeConversion(time_conversion);
    IPPCheck::assertCheck(status == IRM_OK);


    SetupPhreeqcVerbosity setupVerbosity(m_phreeqc, m_isPrintChemistryOn);
    setupVerbosity.run();


    const IPPConfig& ippConfig = conf.ippConfig;
    status = m_phreeqc.LoadDatabase(ippConfig.getDataBaseRelativePath());


    if (status != IRM_OK)
    {
        std::cerr << m_phreeqc.GetErrorString(); // retrieve error messages if needed
        throw PhreeqcRMStop();
    }

    const std::vector<Domain>& domains = ippConfig.domains;


    // init main proc diff calc here, children must happen later
    const std::vector<std::string>& compNames = m_simData->getCompNames();
    m_reacData.porosCalc.setCompNames(&compNames);


    m_reacData.porosCalc.setPhaseNames(&phaseNames);
    m_reacData.porosCalc.setPhaseInfos(&phaseInfos);


    InitPorosity::run(domains, conf.indexCorr, nCellsGlobal,
                      phaseInfos,
                      m_reacData.phaseMolarVolumes,
                      m_reacData.isPermeablePhase,
                      phaseNames,
                      m_reacData.porosCalc,
                      m_reacData.diffCalc,
                      initData.tmpGlobalPoros,
                      initData.tmpGlobalPorosCapillary,
                      initData.tmpGlobalDiff,
                      initData.tmpGlobalInertVolFrac);




    // set porosity before initRun in order to adjust amount of available water
    // ensure minumum of water to ensure numerical stability
    const double minPorosPhreeqc = 1.0E-4;
    if(ippConfig.porosLow < minPorosPhreeqc)
    {
        std::cout << "WARNING: Minimum porosity may too low for phreeqc: " << ippConfig.porosLow << std::endl;
    }

}

void PhreeqcReactionModule::initRootRun(const ReactionModuleConfig &conf, PhreeqcGlobalInitData& initData)
{

    assert(MPIManager::getInstance().isMainProc());

    const std::string mainRunFilePrefix = m_phreeqc.GetFilePrefix();
    // Run file to define solutions and reactants for initial conditions + selected output
    PhreeqcInitRun::run(m_phreeqc, *m_simData, m_reacData.stoichFac, mainRunFilePrefix,
                        conf, initData.featureMask, initData.genericCellID);



    m_phreeqc.SetTime(0.0);
    m_phreeqc.SetTimeStep(0.0);


}


void PhreeqcReactionModule::convergePorosity()
{
    MPI_CHECK_SYNC;

    double globalMax;

    do
    {
        {
            std::stringstream ss;
            ss << MPIManager::getInstance().getRank() << ": porosity not converged yet" << std::endl;
            std::cout << ss.str();
        }

        const std::vector<double>& currPorosPhysical = m_simData->getPorosity();
        const std::vector<double> oldPoros = currPorosPhysical; // copy intended

        m_phreeqc.setPorosityLocal(currPorosPhysical);

        const std::vector<double>& conc = m_simData->getPostReacConc();
        IRM_RESULT status = m_phreeqc.setConcentrationsLocal(conc);
        assert(status == IRM_OK);
        (void) status;

        MPI_CHECK_SYNC;

        this->swapToOld();
        this->runCells();
        this->collectPorosityConvergenceData();



        double diffMax = 0.0;
        for(size_t i = 0; i < currPorosPhysical.size(); ++i)
        {
            const double diff = currPorosPhysical[i] - oldPoros[i];
            const double absDiff = std::abs(diff);
            if(absDiff > std::abs(diffMax))
            {
                diffMax = diff;
            }
        }

        double absDiffMax = std::abs(diffMax);

        MPI_Allreduce(&absDiffMax, &globalMax, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        {
            std::ios_base::fmtflags f( std::cout.flags() );

            std::stringstream ss;
            ss << MPIManager::getInstance().getRank() << ": maxDiff: "
               << std::scientific << diffMax << std::endl;
            std::cout << ss.str();

            pcout << "diff max: " << std::scientific << globalMax << std::endl;

            std::cout.flags( f );
        }
    }
    while(globalMax > 1.0E-6);

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::updatePreCondition()
{
    MPI_CHECK_SYNC;

    const std::vector<double>& porosity = m_simData->getPorosity();
    m_phreeqc.setPorosityLocal(porosity);

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::updatePostCondition()
{
    MPI_CHECK_SYNC;

    if(m_phreeqcData.siFunc)
    {
        PhreeqcModify::updateTargetSaturationIndices(*m_phreeqcData.siFunc,
                                                     m_simData->getEqPhaseNames(),
                                                     m_phreeqc.getNumberCells(),
                                                     m_phreeqc,
                                                     m_currTargetSImap);
    }

    assert(m_phreeqcData.dissolveOnlyFunc);
    PhreeqcModify::updateDissolveOnlyFlags(*m_phreeqcData.dissolveOnlyFunc,
                                           m_simData->getEqPhaseNames(),
                                           m_phreeqc.getNumberCells(),
                                           m_phreeqc,
                                           m_currDissolveOnlyCells);

    MPI_CHECK_SYNC;
}

void PhreeqcReactionModule::retrieveSaturationIndices(std::vector<double>& satIndices)
{
    MPI_CHECK_SYNC;

    IRM_RESULT status =
            m_phreeqc.SetCurrentSelectedOutputUserNumber(PhreeqcConstants::
                                                         SO_PhasesSaturationIndices);
    IPPCheck::assertCheck(status == IRM_OK);

    std::vector<double> so;
    status = m_phreeqc.getSelectedOutputLocal(so);
    IPPCheck::assertCheck(status == IRM_OK);

    assert(so.size() == satIndices.size());

    // TODO: check sorting of phases -> alphabetical


    const FieldDecomposition* decomp = m_simData->getDecomp();
    const size_t nCellsLocal = decomp->getLocalNumberOfCells();
    const size_t nPhases = m_simData->getPhaseNames().size();
    assert(satIndices.size() == nCellsLocal * nPhases);
    GeometryTools::transpose(so, nCellsLocal, nPhases);

    for(size_t iCell = 0; iCell < nCellsLocal; ++iCell)
    {
        if((bool)m_enabledCells[iCell] == true)
        {
            const size_t offset = iCell * nPhases;
            const auto begin = so.begin() + offset;
            const auto end = begin + nPhases;
            const auto out = satIndices.begin() + offset;
            std::copy(begin, end, out);
        }
    }

}

void PhreeqcReactionModule::retrieveMonomersConc()
{
    MPI_CHECK_SYNC;

    IRM_RESULT status =
            m_phreeqc.SetCurrentSelectedOutputUserNumber(PhreeqcConstants::
                                                         SO_NucleationMonomers);
    IPPCheck::assertCheck(status == IRM_OK);

    std::vector<double> so;
    status = m_phreeqc.getSelectedOutputLocal(so);
    IPPCheck::assertCheck(status == IRM_OK);

    assert(so.size() == m_phreeqcData.monomerConc.size());

    const FieldDecomposition* decomp = m_simData->getDecomp();
    const size_t nCellsLocal = decomp->getLocalNumberOfCells();
    assert(so.size() % nCellsLocal == 0);
    const size_t nMonomers = so.size() / nCellsLocal;

    GeometryTools::transpose(so, nCellsLocal, nMonomers);


    for(size_t iCell = 0; iCell < nCellsLocal; ++iCell)
    {
        if((bool)m_enabledCells[iCell] == true)
        {
            const size_t offset = iCell * nMonomers;
            const auto begin = so.begin() + offset;
            const auto end = begin + nMonomers;
            const auto out = m_phreeqcData.monomerConc.begin() + offset;
            std::copy(begin, end, out);
        }
    }

    MPI_CHECK_SYNC;
}


void PhreeqcReactionModule::correctByPorosChange(std::vector<double>& concVec)
{
    assert(PhreeqcConstants::s_calcH2OSeperately);

    // concentrations must be adapted to porosity/volume change
    const std::vector<double>& currPorosityVec = m_simData->getPorosity();

    const size_t nCells = currPorosityVec.size();
    const size_t nComps = m_phreeqc.GetComponentCount();

    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {
        const double& currPorosity = currPorosityVec[iCell];
        const double& oldPorosity = m_oldPorosPhysical[iCell];

        const double factor = CalcPorosityChangeFactor::calc(oldPorosity, currPorosity);
        assert(factor >= 0.0);

        // exclude water
        for(size_t iComp = 1; iComp < nComps; ++iComp)
        {
            const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iComp);
            double& conc = concVec[index];
            conc *= factor;

            assert(std::isnan(conc) == false);
        }
    }
}




} // end of namespace IPP



