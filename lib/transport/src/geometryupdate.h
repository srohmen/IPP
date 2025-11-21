#ifndef GEOMETRYUPDATE_H
#define GEOMETRYUPDATE_H

#include "permeabiltyflag.h"

namespace IPP
{

namespace GeometryUpdate
{

static const double s_relTreshDiff = 0.001;

enum GoemChange
{
    GS_None,
    GS_LiquidToSolid,
    GS_SolidToLiquid
};

template<typename T, typename MaskType>
GoemChange getChange(const T& threshLower, const T& threshUpper,
                     const T& value, const MaskType& permFlag)
{
    if((value < threshLower) &&
       (permFlag & PF_isPermeable))
    {
        // liquid -> solid
        return GS_LiquidToSolid;
    }
    else if((value >= threshUpper) &&
       (permFlag & PF_isPermeable) == false)
    {
        // solid -> liquid
        return GS_SolidToLiquid;
    }
    else
    {
        return GS_None;
    }
}

}


}


#endif // GEOMETRYUPDATE_H
