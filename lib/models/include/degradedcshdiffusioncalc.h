#ifndef DEGRADEDCSHDIFFUSIONCALC_H
#define DEGRADEDCSHDIFFUSIONCALC_H

#include "abstractmultiscalediffusioncalc.h"
#include <memory>

namespace IPP
{

class CSHGelDiffusionCoefCalc;
class EffectiveMediaDiffusionCalc;
class AbstractWeightingFunc;

class DegradedCSHDiffusionCalc : public AbstractMultiScaleDiffusionCalc
{
public:  
    DegradedCSHDiffusionCalc(const CSHGelDiffusionCoefCalc* gelDiffCalc,
                             const EffectiveMediaDiffusionCalc* inclusionModel,
                             const double& freeWaterDiffCoef);

    virtual ~DegradedCSHDiffusionCalc();

    virtual double calc(const SimplePorosityInfos &input) const override;

private:
    double calcGelInclusionD(const double& gelPorosNitro, const double& volCSH,
                             const double& volNonPerm, const double& volPoros) const;


    static const double intrGelPorNitroLD;
    static const double intrGelPorNitroHD;
    static const double degradPoreSlope;
    static const double degradPoreIntercept;
    static const double referenceVol; // cm3

    const std::unique_ptr<const CSHGelDiffusionCoefCalc> m_gelDiffCalc;
    const std::unique_ptr<const EffectiveMediaDiffusionCalc> m_inclusionModel;

    const double m_freeWaterD;

};

} // end of namespace IPP


#endif // DEGRADEDCSHDIFFUSIONCALC_H
