#include "dissolveonlybase.h"

#include <stdexcept>

namespace IPP
{


DissolveOnlyBase::DissolveOnlyBase()
{

}

void DissolveOnlyBase::init(const PhaseNameToInfos& /*phaseInfos*/,
                            const double& /*cellVol*/,
                            const double& /*cellArea*/,
                            const double& /*dt*/)
{
    // do nothing
}

bool DissolveOnlyBase::evaluate(const DissolveOnlyInfos &infos,
                                const iterator &begin,
                                const iterator &end) const
{
    bool allDissolveOnly = end - begin > 0 ? true : false;

    for(iterator it = begin; it != end; ++it)
    {
        const size_t index = it - begin;
        const PreventPrecipResult result = this->preventPrecip(infos, index);

        switch(result)
        {
        case PPR_AllowPrecip:
            *it = false;
            allDissolveOnly = false;
            break;

        case PPR_PreventPrecip:
            *it = true;
            break;

        case PPR_AllAllowPrecip:
            std::fill(begin, end, false);
            allDissolveOnly = false;
            return allDissolveOnly;
            break;

        case PPR_AllPreventPrecip:
            std::fill(begin, end, true);
            allDissolveOnly = true;
            return allDissolveOnly;
            break;

        default:
            throw std::runtime_error("unknown case");
        }
    }

    return allDissolveOnly;
}

}
