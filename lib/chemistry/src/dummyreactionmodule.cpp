#include "dummyreactionmodule.h"

#include "phreeqcprerunhandler.h"
#include "initporosity.h"
#include "mpitools.h"
#include "mpimanager.h"

#include "reactionmoduleconfig.h"
#include "ippconfig.h"
#include "simulationexchangedata.h"
#include "phaseidperelement.h"
#include "printphaseinfos.h"

#include "abstractporositycalc.h"
#include "abstractmultiscalediffusioncalc.h"
#include "fielddecomposition.h"

namespace IPP
{

DummyReactionModule::DummyReactionModule(const size_t nCells, const ReactionModuleConfig& conf)
    : m_nxyz(nCells),
      m_nTracer()
{
    this->init(conf);
}

void DummyReactionModule::init(const ReactionModuleConfig &conf)
{
    const IPPConfig& ippConfig = conf.ippConfig;
    m_simData = conf.data;
    m_nTracer = ippConfig.results.diffEffInfos.size();

    BoundaryConditions::ElementNameToBCVec& diffBCs = conf.diffBCs;
    std::vector<std::string> compNames;
    for(const IPPConfig::Results::DiffEffInfo& info : ippConfig.results.diffEffInfos)
    {
        const std::string& compName = info.tracerName;
        compNames.push_back(compName);

        // just add BC
        diffBCs.insert({compName, BoundaryConditions::BCVec()});
    }

    m_simData->setCompNames(compNames);



    // dummys
    BoundaryConditions::ElementNameToBCVec diffBCDummy;
    PhaseIdPerElement phaseIDsPerElementDummy;
    PhreeqcPrerunHandler::getAuxDataFromPhreeqC(conf.ippConfig, conf.chemOutDir,
                                                false, compNamesDummy, phaseNamesDummy,
                                                phaseNameToInfosDummy, phaseIDsPerElementDummy,
                                                diffBCDummy, conf.inertComposition);


    assert(phaseNamesDummy.size() == phaseNameToInfosDummy.size());
    if(MPIManager::getInstance().isMainProc())
    {
        PrintPhaseInfos::print(compNamesDummy, phaseNamesDummy, phaseNameToInfosDummy);
    }


    std::vector<double> phaseMolarVolumes;
    InitPorosity::extractPhasesMolarVolumes(phaseNameToInfosDummy, phaseNamesDummy, phaseMolarVolumes);


    MpiTools::syncFromRoot(conf.ippConfig.globalComm, phaseNamesDummy);

    AbstractPorosityCalcPtr porosCalc = ippConfig.porosCalc;
    AbstractMultiScaleDiffusionCalcPtr diffCalc = ippConfig.multiScaleDiffusionCalc;
    std::vector<double> globalPorosities;
    std::vector<double> globalCapillaryPorosity;
    std::vector<double> globalDiffCoefs;
    std::vector<double> inertVolFrac;

    if(MPIManager::getInstance().isMainProc())
    {
        // init main proc here, children must happen later
        porosCalc->setCompNames(&compNamesDummy);

        const std::vector<Domain>& domains = ippConfig.domains;

        std::vector<char> isPermeablePhase(phaseNamesDummy.size());
        for(size_t iPhase = 0; iPhase < phaseNamesDummy.size(); ++iPhase)
        {
            const std::string phaseName = phaseNamesDummy[iPhase];
            const PhaseInfo& info = phaseNameToInfosDummy.at(phaseName);
            isPermeablePhase[iPhase] = info.isPermeable;
        }


        // init main proc calc here, children must happen later
        porosCalc->setPhaseNames(&phaseNamesDummy);
        porosCalc->setPhaseInfos(&phaseNameToInfosDummy);


        InitPorosity::run(domains, conf.indexCorr, m_nxyz,
                          phaseNameToInfosDummy,
                          phaseMolarVolumes, isPermeablePhase,
                          phaseNamesDummy,
                          *porosCalc, *diffCalc,
                          globalPorosities,
                          globalCapillaryPorosity,
                          globalDiffCoefs, inertVolFrac);
    }

    InitPorosity::syncData(globalPorosities,
                           globalCapillaryPorosity,
                           globalDiffCoefs,
                           inertVolFrac, *m_simData);


    // init children diff calc as well
    if(MPIManager::getInstance().isMainProc() == false)
    {
        porosCalc->setCompNames(&compNamesDummy);
        porosCalc->setPhaseInfos(&phaseNameToInfosDummy);
        porosCalc->setPhaseNames(&phaseNamesDummy);
    }


    // TODO: for what purpose enabled cells is needed?
    const size_t nCells = m_simData->getDecomp()->getLocalNumberOfCells();


    std::vector<double>& conc = m_simData->getPostReacConc();
    conc.resize(nCells * m_nTracer, 0.0);
}

ConfigurePhreeqcCellsActivity* DummyReactionModule::getConfigCellsActivity()
{
    return nullptr;
}

AbstractGeometrySync *DummyReactionModule::getGeomSync()
{
    return nullptr;
}

void DummyReactionModule::saveCheckpoint()
{

}

void DummyReactionModule::loadCheckpoint(const size_t iteration)
{

}

void DummyReactionModule::errorDump()
{

}

void DummyReactionModule::run()
{

}

void DummyReactionModule::updateComponentConcentrations()
{

}

void DummyReactionModule::updatePreCondition()
{

}

void DummyReactionModule::updatePostCondition()
{

}



}
