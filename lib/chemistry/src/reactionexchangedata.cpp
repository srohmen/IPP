#include "reactionexchangedata.h"

// #include "abstractmultiscalediffusioncalc.h"

namespace IPP
{

ReactionExchangeData::ReactionExchangeData(SimulationExchangeData& simData,
                                           AbstractPorosityCalc& porosCalc,
                                           AbstractMultiScaleDiffusionCalc &diffCalc)
    : simData(simData)
    , porosCalc(porosCalc)
    , diffCalc(diffCalc)
    , invTotalCellVolume(-1.0)
{

}

ReactionExchangeData::~ReactionExchangeData()
{

}

}
