#ifndef DISSOLVEONLYUTILS_H
#define DISSOLVEONLYUTILS_H

#include <vector>
#include <cassert>
#include <cmath>

#include "dissprecipbehaviour.h"
#include "abstractdissolveonlycalc.h"
#include "porositythresholds.h"
#include "dissolveonlyinfos.h"

namespace IPP
{

struct DissolveOnlyInfos;
struct PorosityThresholds;

namespace DissolveOnlyUtils
{

inline bool isDissolveOnlyBehav(const std::vector<DissPrecipBehaviour> &behaves,
                                const size_t phaseIndex)
{
    assert(behaves.size() > phaseIndex);
    const DissPrecipBehaviour& behav = behaves[phaseIndex];

    if(behav == DPB_DissOnly)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline bool isBelowThresh(const double &porosThresh,
                          const double &porosity)
{
    if(porosity <= porosThresh)
    {
        return true;
    }
    else
    {
        return false;
    }

}


inline bool isNeighFilled(const PorosityThresholds &thresh,
                          const DissolveOnlyInfos &infos)
{
    if(infos.neighInfo.location == CellNeighborInfo::L_Surface)
    {
        if(infos.neighInfo.shape == CellNeighborInfo::S_Concave)
        {
            return true;
        }

        if(infos.neighInfo.noPrecip)
        {
            const double r = (double)std::rand() / (double)RAND_MAX;
            const double& poros = infos.porosity ;

            const double dPhi = 0.2;
            const double k = M_LN2 / (1.0 - dPhi);
            const double prop = std::exp( k * (poros - dPhi) ) - 1.0;

            if(r < prop)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return true;
        }
    }

    if(infos.neighInfo.location == CellNeighborInfo::L_Edge)
    {
        const bool isConvexEffect = infos.neighInfo.convexVolume < thresh.neighEdgeVolFrac;
        if(isConvexEffect)
        {
            if(infos.neighInfo.noPrecip)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }

    if(infos.neighInfo.location == CellNeighborInfo::L_Corner)
    {
        const bool isConvexEffect = infos.neighInfo.convexVolume < thresh.neighCornerVolFrac;

        if(isConvexEffect)
        {
            if(infos.neighInfo.noPrecip)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }

    const bool isOwnPorosLow = infos.porosity < thresh.porosUpper;
    if(isOwnPorosLow)
    {
        if(infos.neighInfo.noPrecip)
        {
            return false;
        }
        else
        {
            return true;
        }
    }


    return false;

}


}

}

#endif // DISSOLVEONLYUTILS_H
