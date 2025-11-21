#include "ippimpl.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <boost/filesystem.hpp>


#include "mpimanager.h"

#include "ipprenderer.h"

#include "transreacfactory.h"
#include "ippcoupling.h"
#include "transportmodule.h"
#include "reactionmodule.h"

#include "simulationexchangedata.h"


#include "cellidconverttosifunc.h"

#include "bench_tools.h"
#include "ippexception.h"
#include "mpitools.h"


#include "latticeboltzmannunitconversion.h"
#include "ippstream.h"

#include "configurephreeqccellsactivity.h"
#include "mpitools.h"
#include "scopedfloatingpointexception.h"

#include "abstractgeometrysync.h"
#include "createfilename.h"
#include "ippconstants.h"
#include "filesystemutils.h"

namespace IPP
{

typedef std::chrono::milliseconds MilliSeconds;

static void printStatus(const size_t it)
{
    IPP::pcout << "starting iteration: " << it << std::endl;
}


IPPimpl::IPPimpl()
    : m_conf()
    , m_state(nullptr)
    , m_lastPorsRefChangeIt(0)
    , m_oldPorosRef(1.0)
{

}

IPPimpl::~IPPimpl()
{
    delete m_state;
    m_state = nullptr;
}

void IPPimpl::writeResult() const
{
    // ugly special case... because when chemistry is not enabled the
    // data are never collected from transport code. for now we do it here.
    // in future there should be a requesting queue or something handling this
    // additionally, collect data is doing more than
    // just collecting the current concentration
    if(m_conf->isChemistryEnabled == false ||
            m_conf->transportOnlyFlag == IPPConfig::TOF_InitChemOnly)
    {
        m_state->transModule->collectData();
    }


    ResultsToWrite toWrite = m_conf->results.resultsToWrite;

    // writing all results for first initial iteration
    if(m_state->data->getIteration() == 0)
    {
        for(size_t i = 0; i < RT_MaxRT; ++i)
        {
            const ResultType type = static_cast<ResultType>(i);
            toWrite.add(type);
        }

        if(m_conf->flowFunc->isEnabled() == false)
        {
            toWrite.remove(RT_FluidDensityAndFlux);
        }
    }


    m_state->transModule->writeResults(toWrite);

}

static void createTimePath(const std::string& outputDir,
                           const size_t iteration,
                           boost::filesystem::path& checkpointPath,
                           boost::filesystem::path& timePath)
{
    using namespace boost::filesystem;
    const path relativePath(outputDir);
    const path outPath = absolute(relativePath);
    const path checkpointDir = outPath / IPPConstants::s_checkpointOutput;

    const std::string iterationName = CreateFileName::create(iteration);
    checkpointPath = checkpointDir / iterationName;

    static const std::string s_timeFileName = "time";
    timePath = checkpointPath / s_timeFileName;
}


void IPPimpl::saveCheckpoint() const
{
    if(MPIManager::getInstance().isMainProc())
    {
        boost::filesystem::path checkpointPath, timePath;
        createTimePath(m_conf->outputDir, m_state->data->getIteration(),
                       checkpointPath, timePath);

        FileSystemUtils::createIfNotExisting(checkpointPath);

        const double time = m_state->data->getTime();
        std::ofstream file(timePath.string());
        file << time;
    }

    m_state->transModule->saveCheckpoint();
    m_state->reacModule->saveCheckpoint();
}

void IPPimpl::initialStep()
{

    // write initial state
    // if(m_conf->outputFrequency > 0)
    {
        this->writeResult();
    }


    pcout << "starting initial iteration" << std::endl;

    printStatus(m_state->data->getIteration());
    auto start = std::chrono::steady_clock::now();

    if(m_conf->isChemistryEnabled)
    {
        ConfigurePhreeqcCellsActivity* configActivity = m_state->reacModule->getConfigCellsActivity();
        configActivity->setGloballyEnabledState(true);
    }

    // TODO: needed here? updatePrecondition is invoked in IPPCoupling, too
    m_state->reacModule->updatePreCondition();


    MPI_CHECK_SYNC;


    // sync transport to chemistry and back due to boundary conditions
    IPPCoupling::runReactionModule(*(m_state->transModule), *(m_state->reacModule), true);

    if(m_conf->isChemistryEnabled)
    {
        ConfigurePhreeqcCellsActivity* configActivity = m_state->reacModule->getConfigCellsActivity();
        configActivity->setGloballyEnabledState(false);
    }


    const auto duration = CalcDuration::calc(start);
    pcout << "\t-> iteration " << m_state->data->getIteration()
          << " finished, calc took:\t" << duration.count()
          << " ms" << std::endl
          << "\ttime: " << m_state->data->getTime();


    const ResultsToWrite toWrite = m_conf->results.resultsToWrite;
    m_state->transModule->writeDebugData(toWrite);
}

void IPPimpl::init(const IPPConfigPtr conf)
{
    // take ownership of config
    m_conf = conf;

    // TODO: check if mpi is needed otherwise implement fallback
    int argc = conf->argc;
    char** argv = conf->argv;
    MPIManager::getInstance().init(argc, argv, m_conf->globalComm);


    const int myRank = MPIManager::getInstance().getRank();
    const pid_t pid = getpid();
    std::cout << "rank to PID: " << myRank << " -> " << pid << std::endl;


    const int nProcs = MPIManager::getInstance().getNProcs();


    if(m_conf->decomp.empty() == false)
    {
        if((int)m_conf->decomp.size() != nProcs)
        {
            throw std::runtime_error("number of decomp boxes are not equal to number of procs: "
                                     + std::to_string(m_conf->decomp.size())
                                     + " vs. " + std::to_string(nProcs));
        }
    }

    m_oldPorosRef = m_conf->porosRef;

    m_state = TransReacFactory::create(*m_conf);

    initialStep();

    if(m_conf->isChemistryEnabled)
    {
        // printing all the solid -> liquid / new interface changes
        AbstractGeometrySync* geomSync = m_state->reacModule->getGeomSync();
        assert(geomSync);
        geomSync->setPrintDebug(true);
    }

    if(m_conf->transportOnlyFlag == IPPConfig::TOF_FakeChem)
    {
        ConfigurePhreeqcCellsActivity* configActivity = m_state->reacModule->getConfigCellsActivity();
        configActivity->setGloballyDisabledState(true);
    }

}

void IPPimpl::loadCheckpoint(const size_t iteration)
{    
    m_lastPorsRefChangeIt = iteration;
    m_state->data->setIteration(iteration);

    {
        boost::filesystem::path checkpointPath, timePath;
        createTimePath(m_conf->outputDir, iteration, checkpointPath, timePath);

        IPPCheck::assertCheck(exists(checkpointPath),
                              "Checkpoint folder is not existing: "
                              + checkpointPath.string());

        IPPCheck::assertCheck(exists(timePath), "time file not found: "
                              + timePath.string());

        double time;
        if(MPIManager::getInstance().isMainProc())
        {
            std::ifstream file(timePath.string());
            std::string timeStr;
            file >> timeStr;
            time = std::stod(timeStr);
        }

        const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<double>();
        MPI_Bcast(&time, 1, mpiType, 0, m_conf->globalComm);
        m_state->data->initTime(time);
    }


    m_state->reacModule->loadCheckpoint(iteration);
    m_state->transModule->loadCheckpoint(iteration);

    pcout << "loading checkpoint fully completed" << std::endl;
}

bool IPPimpl::mustWriteResults(const size_t it) const
{
    return m_conf->outputFrequency > 0 && it % m_conf->outputFrequency == 0;
}

void IPPimpl::runStep()
{
    {
        auto start = std::chrono::steady_clock::now();
        m_state->transModule->run();

        if(m_conf->verbose)
        {
            const auto duration = CalcDuration::calc(start);
            IPP::pcout << "\t-> LBM transport took:\t" << duration.count() << " ms" << std::endl;
        }
    }

    if(m_conf->isChemistryEnabled && m_conf->transportOnlyFlag != IPPConfig::TOF_InitChemOnly)
    {
        auto start = std::chrono::steady_clock::now();

        const size_t it = m_state->data->getIteration();
        const bool enableAllCells = mustWriteResults(it);

        if(enableAllCells)
        {
            ConfigurePhreeqcCellsActivity* configActivity = m_state->reacModule->getConfigCellsActivity();
            configActivity->setGloballyEnabledState(true);
        }

        IPPCoupling::runReactionModule(*(m_state->transModule), *(m_state->reacModule), m_conf->verbose);



        if(enableAllCells)
        {
            ConfigurePhreeqcCellsActivity* configActivity = m_state->reacModule->getConfigCellsActivity();
            configActivity->setGloballyEnabledState(false);
        }

        if(m_conf->verbose)
        {
            const auto duration = CalcDuration::calc(start);
            IPP::pcout << "\t-> total chemistry took:\t" << duration.count() << " ms" << std::endl;
        }
    }


}

void IPPimpl::setPorosRef(const double &newPorosRef)
{
    IPP::pcout << "update porosRef to: " << newPorosRef << std::endl;
    m_state->transModule->updateTransportScalarDref(newPorosRef);

    const double& dx = m_conf->spatialResolution;
    const double& D = m_conf->diffusionCoefReference;
    const double& tau = m_conf->tauRef;
    const double dt = LatticeBoltzmannUnitConversion::calcTimeStepPTRT(dx, D, tau,
                                                                       m_conf->porosLow,
                                                                       newPorosRef,
                                                                       m_conf->is3DSimulation);
    IPP::pcout << "new dt: " << dt << std::endl;
    m_state->data->setDtCurr(dt);
}

bool IPPimpl::updatePorosRef(const size_t it)
{
    // run at least several steps before change poros ref
    const size_t diff = it - m_lastPorsRefChangeIt;
    if(it > 0 && diff > 10000)
    {
        if(m_conf->fluxConvergenceTolFactor > 1.0 && m_oldPorosRef < 1.0)
        {
            double newPorosRef = m_conf->fluxConvergenceTolFactor * m_oldPorosRef;
            if(newPorosRef >= 1.0)
            {
                IPP::pcout << "last step of flux conversion. setting porosRef to 1.0. resetting flux conv tol: "
                           << m_state->fluxConvergenceTol << " -> " << m_conf->fluxConvergenceTolHigh
                           << std::endl;
                newPorosRef = 1.0;
                m_state->fluxConvergenceTol = m_conf->fluxConvergenceTolHigh;
            }

            setPorosRef(newPorosRef);

            m_lastPorsRefChangeIt = it;
            m_oldPorosRef = newPorosRef;
            return false;
        }
        else
        {
            return true;
        }
    }

    return false;
}

void IPPimpl::execute()
{
    auto lastCheckpoint = std::chrono::steady_clock::now();
    auto lastStatus = std::chrono::steady_clock::now();

    const double maxTime = m_conf->maxTime;
    const size_t startIt = m_state->data->getIteration();
    size_t lastStatIt = startIt;

    for(size_t it = startIt; it < m_conf->maxIter || m_state->data->getTime() < maxTime; )
    {
        // default loop
        ++it;
        m_state->data->setIteration(it);

        if(m_conf->verbose)
        {
            printStatus(it);
        }

        auto start = std::chrono::steady_clock::now();

        this->runStep();

        m_state->data->incrementTime();

        const auto duration = CalcDuration::calc(start);
        const double durationMS = duration.count();

        if(m_conf->verbose)
        {
            IPP::pcout << "\t-> iteration " << it
                       << " finished, total calc took:\t" << durationMS << " ms" << std::endl
                       << "\ttime: " << std::setprecision (15) << m_state->data->getTime() << std::endl;
        }

        // print some status each 10 s
        const auto lastStatduration = CalcDuration::calc(lastStatus);
        const double lastStatDurMS = lastStatduration.count();
        if(lastStatDurMS > 10000)
        {
            IPP::pcout << "iterations per second:\t" << 1000 * (it - lastStatIt) / lastStatDurMS << std::endl;
            lastStatus = std::chrono::steady_clock::now();
            lastStatIt = it;
        }

        ///////////////////////
        // update renderer
        if(m_state->renderer.get() != nullptr)
        {
            std::cout << "\tUpdating renderer" << std::endl;
            assert(false);
            // m_mpiFuncs->getCouplingModuleDispatch().updateRenderer();

            if(MPIManager::getInstance().isMainProc())
            {
                m_state->renderer->render();
            }
        }
        ///////////////////////



        const bool lastIteration = it >= m_conf->maxIter;

        {
            // writing results
            if(mustWriteResults(it) || lastIteration)
            {
                MPI_CHECK_SYNC;
                IPP::pcout << "\tWriting result for iteration:\t" << it << std::endl;
                auto startResult = std::chrono::steady_clock::now();
                this->writeResult();
                auto durationResult = std::chrono::duration_cast<MilliSeconds>(std::chrono::steady_clock::now()
                                                                               - startResult);
                IPP::pcout << "\t-> writing result took:\t" << durationResult.count() << " ms" << std::endl;
            }
            else if(m_conf->outletFluxFrequency > 0
                    && it % m_conf->outletFluxFrequency == 0)
            {
                // only write the outlet flux if not have written all results
                // all results writer does already print the outlet flux

                MPI_CHECK_SYNC;
                IPP::pcout << "\tWriting outlet flux" << std::endl;
                auto startResult = std::chrono::steady_clock::now();
                m_state->transModule->writeOutletFlux();
                auto durationResult = std::chrono::duration_cast<MilliSeconds>(std::chrono::steady_clock::now()
                                                                               - startResult);
                IPP::pcout << "\t-> writing outlet flux took:\t" << durationResult.count() << " ms" << std::endl;
            }

            // write debug data anyhow
            const ResultsToWrite& toWrite = m_conf->results.resultsToWrite;
            m_state->transModule->writeDebugData(toWrite);
        }
        ///////////////////////



        // writing checkpoint ?
        bool checkPointFreqHit = false;
        if(m_conf->checkpointingFreq > 0)
        {
            checkPointFreqHit = it % m_conf->checkpointingFreq == 0;
        }

        auto durationSinceLastCheckpoint = std::chrono::duration_cast<MilliSeconds>(std::chrono::steady_clock::now()
                                                                                    - lastCheckpoint);
        const double totalDurationMS = durationSinceLastCheckpoint.count();
        const double totalDurationMinutes = totalDurationMS / (1000 * 60);
        const bool checkPointTimeOver = totalDurationMinutes > m_conf->checkpointingTime;

        bool writeCheckPoint = checkPointFreqHit || checkPointTimeOver || lastIteration || m_state->isConverged;

        // sync from root process since time criteria could be different for difference processes
        // can maybe optimized with using MPI timer check
        const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<bool>();
        MPI_Bcast(&writeCheckPoint, 1, mpiType, 0, m_conf->globalComm);

        if(writeCheckPoint)
        {
            IPP::pcout << "\tWriting checkpoint" << std::endl;
            auto startCheckpoint = std::chrono::steady_clock::now();
            saveCheckpoint();
            auto durationCheckpoint = std::chrono::duration_cast<MilliSeconds>(std::chrono::steady_clock::now()
                                                                               - startCheckpoint);
            IPP::pcout << "\t-> saving checkpoint took:\t" << durationCheckpoint.count() << " ms" << std::endl;

            lastCheckpoint = std::chrono::steady_clock::now();
        }
        ///////////////////////


        if(m_state->isConverged)
        {            
            // hardcoded tauref optim
            const bool allowFinishing = updatePorosRef(it);

            if(allowFinishing)
            {
                pcout << "convergence reached!" << std::endl;
                break;
            }
            else
            {
                m_state->isConverged = false;
            }
        }

    }
}

void IPPimpl::setRenderingModule(IPPRendererPtr renderer)
{
    if(MPIManager::getInstance().isMainProc())
    {
        m_state->renderer = renderer;
    }
}


}
