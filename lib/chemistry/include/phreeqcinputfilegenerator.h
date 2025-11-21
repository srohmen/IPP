#ifndef PHREEQCINPUTFILEGENERATOR_H
#define PHREEQCINPUTFILEGENERATOR_H

#include <sstream>

#include "ippconfig.h"
#include "phreeqcfeatures.h"

namespace IPP
{

namespace PhreeqcInputFileGenerator
{

void generateInitString(const IPPConfig& config,
                        const std::vector<std::string>& allPhases,
                        const bool enableTotals, const bool enableSaturationIndex,
                        CellIdToPhreeqcFeatureMask& enabledFeatures,
                        size_t &genericCellID,
                        std::stringstream& compositionSS,
                        std::stringstream& outputSS);

void generateSolutionBody(const Composition& comp,
                          const std::map<ElementSpecies, double>& elemNameToAmount,
                          std::stringstream& ss);

void generateElementConcLine(const std::string& elemName,
                             const double& conc,
                             const std::string& specName,
                             std::stringstream& ss);

void generateEqPhaseBody(const Composition::PhaseDefVec& phases,
                         const AbstractDissPrecipOnlyInfo& phaseBehav,
                         std::stringstream& ss);

void generateInertSoluteDatabaseExt(const std::vector<IPPConfig::Results::DiffEffInfo>& diffEffInfos,
                                    std::stringstream& ss);

void generateConvergenceHacksDefault(std::stringstream& ss);

void generateConvergenceHacksInitRun(std::stringstream& ss);

void generateConvergenceHacks(const size_t iterations, const double& tolerance,
                              const double& convergence_tolerance,
                              const bool diagonalScaling,
                              std::stringstream& ss);

}

} // end of namespace IPP

#endif // PHREEQCINPUTFILEGENERATOR_H
