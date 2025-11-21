#ifndef PRINTPHASEINFOS_H
#define PRINTPHASEINFOS_H

#include <vector>
#include "phasenametoinfos.h"

namespace IPP
{

namespace PrintPhaseInfos
{

void print(const std::vector<std::string>& compNames,
           const std::vector<std::string>& phaseNames,
           const PhaseNameToInfos& phaseInfos);

}


}


#endif // PRINTPHASEINFOS_H
