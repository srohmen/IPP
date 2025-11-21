#ifndef REACTIONEXCHANGEDATA_H
#define REACTIONEXCHANGEDATA_H

#include "simulationexchangedata.h"
#include "stoichfactors.h"
#include "abstractporositycalc.h"
#include "abstractmultiscalediffusioncalcfwd.h"

namespace IPP
{

struct ReactionExchangeData
{
    ReactionExchangeData(SimulationExchangeData& simData,
                         AbstractPorosityCalc &porosCalc,
                         AbstractMultiScaleDiffusionCalc& diffCalc);
    ~ReactionExchangeData();

    SimulationExchangeData& simData;

    AbstractPorosityCalc& porosCalc;
    AbstractMultiScaleDiffusionCalc& diffCalc;

    double invTotalCellVolume;
    StoichFactorsPerCompID stoichFac;
    std::vector<double> phaseMolarVolumes; // mol / m^3 * cell volume
    std::vector<char> isPermeablePhase;

};


}

#endif // REACTIONEXCHANGEDATA_H
