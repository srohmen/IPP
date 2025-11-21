#ifndef PHREEQCPRERUNHANDLER_H
#define PHREEQCPRERUNHANDLER_H

#include <boost/filesystem/path.hpp>

#include "boundaryconditions.h"
#include "phasenametoinfos.h"

namespace IPP
{

class IPPConfig;
class PhaseIdPerElement;

namespace PhreeqcPrerunHandler
{


void getAuxDataFromPhreeqC(const IPPConfig& conf,
                           const boost::filesystem::path& outDir,
                           const bool calcH2OSeperately,
                           std::vector<std::string> &compNames,
                           std::vector<std::string>& phaseNames,
                           PhaseNameToInfos& phaseToInfos,
                           PhaseIdPerElement& phaseIDsPerElement,
                           BoundaryConditions::ElementNameToBCVec& resultBCs,
                           std::vector<double>& inertComposition);


}

}

#endif // PHREEQCPRERUNHANDLER_H
