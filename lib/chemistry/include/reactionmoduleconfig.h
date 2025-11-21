#ifndef REACTIONMODULECONFIG_H
#define REACTIONMODULECONFIG_H

#include <memory>

#include <boost/filesystem/path.hpp>

#include "boundaryconditions.h"

namespace IPP
{


class IPPConfig;
class PhreeqcReactionModuleData;
class CellIndexCorrection;
class SimulationExchangeData;
class NonLocalOperations;

class ReactionModuleConfig
{
public:
    ReactionModuleConfig(const IPPConfig& geoChemConfig,
                         const CellIndexCorrection& indexCorr,
                         std::shared_ptr<SimulationExchangeData>& data,
                         BoundaryConditions::ElementNameToBCVec& diffBCs,
                         std::vector<double>& inertComposition,
                         PhreeqcReactionModuleData& phreeqcData,
                         const boost::filesystem::path& chemOutDir,
                         const boost::filesystem::path& dumpOutDir,
                         NonLocalOperations& nonLocalOperation);

    const IPPConfig& ippConfig;
    const CellIndexCorrection& indexCorr;

    std::shared_ptr<SimulationExchangeData> data;
    BoundaryConditions::ElementNameToBCVec& diffBCs;
    std::vector<double>& inertComposition;

    PhreeqcReactionModuleData& phreeqcData;

    const boost::filesystem::path& chemOutDir;
    const boost::filesystem::path& dumpOutDir;

    NonLocalOperations& nonLocalOperation;
};

}

#endif // REACTIONMODULECONFIG_H
