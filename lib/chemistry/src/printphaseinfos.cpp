#include "printphaseinfos.h"

#include <iostream>
#include <cassert>

namespace IPP
{

void PrintPhaseInfos::print(const std::vector<std::string>& compNames,
                            const std::vector<std::string>& phaseNames,
                            const PhaseNameToInfos& phaseInfos)
{

    std::cout << "Phase infos:" << std::endl
              << "iPhase\tPhaseName\tgfw\tVm";
    for(size_t iComp = 0; iComp < compNames.size(); ++iComp)
    {
        std::cout << "\t" << compNames[iComp];
    }
    std::cout << std::endl;

    for(size_t iPhase = 0; iPhase < phaseNames.size(); ++iPhase)
    {
        const std::string& phaseName = phaseNames[iPhase];
        const PhaseInfo& info = phaseInfos.at(phaseName);

        std::cout << iPhase
                  << "\t" << phaseName
                  << "\t" << info.gfw
                  << "\t" << info.molarVolume;

        assert(info.stoich.size() == compNames.size());
        for(size_t iComp = 0; iComp < info.stoich.size(); ++iComp)
        {
            std::cout << "\t" << info.stoich[iComp];
        }

        std::cout << std::endl;
    }

}



}
