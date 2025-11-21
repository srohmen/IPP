#include "reactionmoduleconfig.h"

namespace IPP
{

ReactionModuleConfig::ReactionModuleConfig(const IPPConfig &geoChemConfig,
                                           const CellIndexCorrection &indexCorr,
                                           std::shared_ptr<SimulationExchangeData> &data,
                                           BoundaryConditions::ElementNameToBCVec &diffBCs,
                                           std::vector<double> &inertComposition,
                                           PhreeqcReactionModuleData& phreeqcData,
                                           const boost::filesystem::path &chemOutDir,
                                           const boost::filesystem::path &dumpOutDir,
                                           NonLocalOperations& nonLocalOperation)

    : ippConfig(geoChemConfig)
    , indexCorr(indexCorr)
    , data(data)
    , diffBCs(diffBCs)
    , inertComposition(inertComposition)
    , phreeqcData(phreeqcData)
    , chemOutDir(chemOutDir)
    , dumpOutDir(dumpOutDir)
    , nonLocalOperation(nonLocalOperation)
{

}


}
