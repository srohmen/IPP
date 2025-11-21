#pragma once

#include "abstractmultiscalediffusioncalc.h"
#include "degradedcshdiffusioncalc.h"

namespace IPP {

class LowPorosThreshDegradedCSHDiffusionCalc : public AbstractMultiScaleDiffusionCalc
{
public:
    LowPorosThreshDegradedCSHDiffusionCalc(const CSHGelDiffusionCoefCalc* gelDiffCalc,
                                           const EffectiveMediaDiffusionCalc* inclusionModel,
                                           const double& freeWaterDiffCoef,
                                           const double& porosThreshClog);

    virtual double calc(const SimplePorosityInfos &input) const override;

private:
    const DegradedCSHDiffusionCalc m_degrDiff;
    const double m_porosThreshClog;

};

} // namespace IPP

