#ifndef LOWPOROSDEGRADEDCSHDIFFUSIONCALC_H
#define LOWPOROSDEGRADEDCSHDIFFUSIONCALC_H

#include "abstractmultiscalediffusioncalc.h"
#include "degradedcshdiffusioncalc.h"

namespace IPP
{

class LowPorosDegradedCSHDiffusionCalc : public AbstractMultiScaleDiffusionCalc
{
public:
    LowPorosDegradedCSHDiffusionCalc(const CSHGelDiffusionCoefCalc* gelDiffCalc,
                                     const EffectiveMediaDiffusionCalc* inclusionModel,
                                     const double& freeWaterDiffCoef,
                                     const double &expo,
                                     const double &porosThreshArchies,
                                     const double &porosThreshClog);

    virtual ~LowPorosDegradedCSHDiffusionCalc() override;

    virtual double calc(const SimplePorosityInfos &input) const override;

private:
    const DegradedCSHDiffusionCalc m_degrDiff;
    const double m_diffCoef;
    const double m_expo;
    const double m_porosThreshArchies;
    const double m_porosThreshClog;
};


} // end of namespace IPP


#endif // LOWPOROSDEGRADEDCSHDIFFUSIONCALC_H
