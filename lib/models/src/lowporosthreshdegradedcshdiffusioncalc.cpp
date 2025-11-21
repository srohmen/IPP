#include "lowporosthreshdegradedcshdiffusioncalc.hpp"

#include "simpleporosityinfos.h"

namespace IPP
{

LowPorosThreshDegradedCSHDiffusionCalc::LowPorosThreshDegradedCSHDiffusionCalc(const CSHGelDiffusionCoefCalc* gelDiffCalc,
                                                                               const EffectiveMediaDiffusionCalc* inclusionModel,
                                                                               const double& freeWaterDiffCoef,
                                                                               const double& porosThreshClog)
    : m_degrDiff(gelDiffCalc, inclusionModel, freeWaterDiffCoef)
    , m_porosThreshClog(porosThreshClog)
{

}

double LowPorosThreshDegradedCSHDiffusionCalc::calc(const SimplePorosityInfos& input) const
{
    if(input.porosityTotal > m_porosThreshClog)
    {
        const double diffDegr = m_degrDiff.calc(input);
        return diffDegr;
    }
    else
    {
        return 0.0;
    }
}

} // namespace IPP
