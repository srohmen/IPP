#ifndef PHREEQCINITRUN_H
#define PHREEQCINITRUN_H

#include <string>
#include "reactionmodule.h"
#include "stoichfactors.h"
#include "phreeqcfeatures.h"

class PhreeqcRM;

namespace IPP
{

class SimulationExchangeData;

namespace PhreeqcInitRun
{

void run(PhreeqcRM& phreeqc, SimulationExchangeData& simData,
         StoichFactorsPerCompID& stoichFac,
         const std::string& mainRunFilePrefix,
         const ReactionModuleConfig& conf,
         CellIdToPhreeqcFeatureMask& featureMask,
         size_t &genericCellID);
}

}

#endif // PHREEQCINITRUN_H
