#include "dissolveonlyonbelowthreshold.h"

#include <cassert>

#include "dissolveonlyutils.h"

namespace IPP
{

BelowThreshDissolveOnly::BelowThreshDissolveOnly(const double& porosLowerThresh)
    : m_porosLowerThresh(porosLowerThresh)
{

}

double BelowThreshDissolveOnly::getLowerPorosityThresh() const
{
    return m_porosLowerThresh;
}

bool BelowThreshDissolveOnly::needsNeighInfos() const
{
    return false;
}

bool BelowThreshDissolveOnly::needsSaturationIndices() const
{
    return false;
}

bool BelowThreshDissolveOnly::isNucleation() const
{
    return false;
}

void BelowThreshDissolveOnly::setPrecipDissolveOnlyBehaviour(const std::vector<DissPrecipBehaviour> &behav)
{
    assert(m_behav.empty());
    m_behav = behav;
}

DissolveOnlyBase::PreventPrecipResult BelowThreshDissolveOnly::preventPrecip(const DissolveOnlyInfos &infos,
                                                                             const size_t index) const
{
    const bool isDissolveOnlyBehav = DissolveOnlyUtils::isDissolveOnlyBehav(m_behav, index);

    if(isDissolveOnlyBehav)
    {
        return PPR_PreventPrecip;
    }
    else
    {
        const bool isBelowThresh = DissolveOnlyUtils::isBelowThresh(m_porosLowerThresh, infos.porosity);        

        if(isBelowThresh)
        {
            return PPR_PreventPrecip;
        }
    }

    return PPR_AllowPrecip;
}

BelowThreshNeighborAwareDissolveOnly::BelowThreshNeighborAwareDissolveOnly(const PorosityThresholds& thresholds)
    : m_thresholds(thresholds)
{

}

double BelowThreshNeighborAwareDissolveOnly::getLowerPorosityThresh() const
{
    return m_thresholds.porosLower;
}

bool BelowThreshNeighborAwareDissolveOnly::needsNeighInfos() const
{
    return true;
}

bool BelowThreshNeighborAwareDissolveOnly::needsSaturationIndices() const
{
    return false;
}

bool BelowThreshNeighborAwareDissolveOnly::isNucleation() const
{
    return false;
}

void BelowThreshNeighborAwareDissolveOnly::setPrecipDissolveOnlyBehaviour(const std::vector<DissPrecipBehaviour> &behav)
{
    assert(m_behav.empty());
    m_behav = behav;
}


DissolveOnlyBase::PreventPrecipResult BelowThreshNeighborAwareDissolveOnly::preventPrecip(const DissolveOnlyInfos &infos,
                                                                                          const size_t index) const
{
    const bool isDissolveOnlyBehav = DissolveOnlyUtils::isDissolveOnlyBehav(m_behav, index);

    if(isDissolveOnlyBehav)
    {
        return PPR_PreventPrecip;
    }
    else
    {
        const bool isBelowThresh = DissolveOnlyUtils::isBelowThresh(m_thresholds.porosLower, infos.porosity);
        if(isBelowThresh)
        {
            return PPR_AllPreventPrecip;
        }
        else
        {
            const bool isNeighFilled = DissolveOnlyUtils::isNeighFilled(m_thresholds, infos);

            if(isNeighFilled)
            {
                return PPR_AllAllowPrecip;
            }
            else
            {
                return PPR_AllPreventPrecip;
            }
        }
    }

}
}

