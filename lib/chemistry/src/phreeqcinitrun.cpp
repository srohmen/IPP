#include "phreeqcinitrun.h"

#include <fstream>
#include <unordered_map>

#include <PhreeqcRM.h>

#include "reactionmoduleconfig.h"
#include "phreeqcinputfilegenerator.h"
#include "ippexception.h"
#include "simulationexchangedata.h"
#include "reactionmoduleoptimization.h"
#include "cellid_dissolveonly_func.h"
#include "phreeqcreactionmoduledata.h"

namespace IPP
{



void PhreeqcInitRun::run(PhreeqcRM& phreeqc,
                         SimulationExchangeData& simData,
                         StoichFactorsPerCompID& stoichFactors,
                         const std::string& mainRunFilePrefix,
                         const ReactionModuleConfig& conf,
                         CellIdToPhreeqcFeatureMask &featureMask,
                         size_t &genericCellID)
{
    // enable totals only when needed since this hurts performance
    const bool needTotals = false;
    const bool needSaturationIndices = conf.phreeqcData.dissolveOnlyFunc->needsSaturationIndices();


    std::stringstream compStrStrm;
    std::stringstream outputStrStrm;
    PhreeqcInputFileGenerator::generateInitString(conf.ippConfig,
                                                  simData.getPhaseNames(),
                                                  needTotals,
                                                  needSaturationIndices,
                                                  featureMask,
                                                  genericCellID,
                                                  compStrStrm,
                                                  outputStrStrm);

    const std::string compInputStr = compStrStrm.str();
    const std::string outputStr = outputStrStrm.str();

    // TODO: write in debug build only?
    std::ofstream f(mainRunFilePrefix + ".pqi");
    IPPCheck::assertCheck(f.is_open());
    f << "DATABASE " << conf.ippConfig.dataBaseName << std::endl;
    f << compInputStr;
    f << outputStr;
    f.close();


    phreeqc.RunString(false, true, false, compInputStr);
    phreeqc.RunString(true, false, false, outputStr);

    // initialize global internal phreeqc values
    const size_t nComps = phreeqc.FindComponents();

    // Print some of the reaction module information
    {
        std::ostringstream oss;
        oss << "Database:                                         " << phreeqc.GetDatabaseFileName() << "\n";
        oss << "Number of threads:                                " << phreeqc.GetThreadCount() << "\n";
        oss << "Number of MPI processes:                          " << phreeqc.GetMpiTasks() << "\n";
        oss << "MPI task number:                                  " << phreeqc.GetMpiMyself() << "\n";
        oss << "File prefix:                                      " << phreeqc.GetFilePrefix() << "\n";
        oss << "Number of grid cells in the user's model:         " << phreeqc.GetGridCellCount() << "\n";
        oss << "Number of chemistry cells in the reaction module: " << phreeqc.GetChemistryCellCount() << "\n";
        oss << "Number of components for transport:               " << phreeqc.GetComponentCount() << "\n";
        oss << "Partioning of UZ solids:                          " << phreeqc.GetPartitionUZSolids() << "\n";
        oss << "Error handler mode:                               " << phreeqc.GetErrorHandlerMode() << "\n";
        phreeqc.OutputMessage(oss.str());

        std::cout << oss.str() << std::endl;
    }


    const std::vector < double > & gfw = phreeqc.GetGfw();
    const std::vector<std::string>& comps = phreeqc.GetComponents();
    const std::vector<std::string>& compNames = simData.getCompNames();


    IPPCheck::assertCheck(comps.size() == compNames.size(),
                          "Component of phreeqc and simulation data are not equal");

    IPPCheck::assertCheck(comps == compNames,
                          "Component of phreeqc and simulation data are not equal");

    std::unordered_map<std::string, size_t> compNameToID;
    for (size_t iComp = 0; iComp < nComps; iComp++)
    {
        std::cout << iComp << "\tM( " << compNames[iComp] << " ) = " << gfw[iComp] << std::endl;

        const std::string& compName = compNames[iComp];
        compNameToID.insert(std::make_pair(compName, iComp));
    }

    const std::vector<std::string> &speciesNames = phreeqc.GetSpeciesNames();
    for (size_t iSpecies = 0; iSpecies < speciesNames.size(); iSpecies++)
    {
        std::cout << iSpecies << "\t" << speciesNames[iSpecies] << std::endl;
    }


    const std::vector<cxxNameDouble>& speciesStoich = phreeqc.GetSpeciesStoichiometry();
    const size_t nSpecies = speciesStoich.size();
    stoichFactors.resize(nComps);
    for (size_t iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    {
        const cxxNameDouble& speciesNameInfo = speciesStoich[iSpecies];

        for (cxxNameDouble::const_iterator it = speciesNameInfo.begin();
             it != speciesNameInfo.end(); ++it)
        {
            const std::string& compName = it->first;
            const double& stoichFac = it->second;

            const size_t compID = compNameToID.at(compName);

            SpeciesIDtoStoichFactorVec& compInfo = stoichFactors[compID];
            compInfo.push_back(SpeciesIDtoStoichFactor(iSpecies, stoichFac));
        }
    }

}


}
