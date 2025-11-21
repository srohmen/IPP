#ifndef PHREEQCRMHELPER_H
#define PHREEQCRMHELPER_H

#include <boost/filesystem/path.hpp>

#include "phreeqcfeatures.h"

class PhreeqcRM;

namespace IPP
{

class ReactionModuleConfig;
struct PhreeqcGlobalInitData;

namespace PhreeqcRMHelper
{

void muteConsoleOutput(PhreeqcRM& phreeqcRM, const boost::filesystem::path& outDir);


void findMappingAndFeatures(const ReactionModuleConfig &conf,
                            const PhreeqcGlobalInitData& initData,
                            const size_t nxyz,
                            std::vector<int>& grid2chem,
                            std::vector<int>& ic,
                            std::vector<size_t>& cellToDomain);


void findFeatures(const size_t nCells,
                  const CellIdToPhreeqcFeatureMask& featureMaskPerDomain,
                  const std::vector<int>& cells,
                  const size_t iDomain,
                  std::vector<int>& ic);

void printNegativePorosWarnings(const std::vector<double>& porosities);



}


}


#endif // PHREEQCRMHELPER_H
