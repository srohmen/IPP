#include "lowporosdegradedcshdiffusioncalc.h"

#include "simpleporosityinfos.h"
#include <cmath>

namespace IPP
{

LowPorosDegradedCSHDiffusionCalc::LowPorosDegradedCSHDiffusionCalc(const CSHGelDiffusionCoefCalc *gelDiffCalc,
                                                                   const EffectiveMediaDiffusionCalc *inclusionModel,
                                                                   const double& freeWaterDiffCoef,
                                                                   const double& expo,
                                                                   const double& porosThreshArchies,
                                                                   const double &porosThreshClog)
    : m_degrDiff(gelDiffCalc, inclusionModel, freeWaterDiffCoef)
    , m_diffCoef(freeWaterDiffCoef)
    , m_expo(expo)
    , m_porosThreshArchies(porosThreshArchies)
    , m_porosThreshClog(porosThreshClog)
{

}

LowPorosDegradedCSHDiffusionCalc::~LowPorosDegradedCSHDiffusionCalc()
{

}

double LowPorosDegradedCSHDiffusionCalc::calc(const SimplePorosityInfos &input) const
{
    if(input.porosityTotal < m_porosThreshClog)
    {
        return 0.0;
    }
    else
    {
        const double diffDeg = m_degrDiff.calc(input);
        const double diffSpheres = m_diffCoef * std::pow(input.porosityTotal, m_expo);

        const double t = std::min(1.0, input.porosityTotal / m_porosThreshArchies);
        const double diffInterm = diffDeg * t + (1.0 - t) * diffSpheres;
        return diffInterm;
    }
}

} // end of namespace IPP

