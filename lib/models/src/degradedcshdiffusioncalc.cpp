#include "degradedcshdiffusioncalc.h"

#include <cassert>
#include <algorithm>

#include "cshporosityinfos.h"
#include "cshgeldiffusioncoefcalc.h"
#include "effectivemediadiffusioncalc.h"
#include "abstractweightingfunc.h"
#include "ippexception.h"


namespace IPP
{

// nitrogen accessable gel porosity of sound CSH
const double DegradedCSHDiffusionCalc::intrGelPorNitroLD = 0.1782;
const double DegradedCSHDiffusionCalc::intrGelPorNitroHD = 0.0;

// gel porosity increases upon degradataion
const double DegradedCSHDiffusionCalc::degradPoreSlope = 0.065 / (0.8 - 1.667);
const double DegradedCSHDiffusionCalc::degradPoreIntercept = -1.667 * DegradedCSHDiffusionCalc::degradPoreSlope;


const double DegradedCSHDiffusionCalc::referenceVol = 1000.0; // cm3



DegradedCSHDiffusionCalc::DegradedCSHDiffusionCalc(const CSHGelDiffusionCoefCalc *gelDiffCalc,
                                                   const EffectiveMediaDiffusionCalc *inclusionModel,
                                                   const double &freeWaterDiffCoef)
    : m_gelDiffCalc(gelDiffCalc)
    , m_inclusionModel(inclusionModel)
    , m_freeWaterD(freeWaterDiffCoef)
{

}

DegradedCSHDiffusionCalc::~DegradedCSHDiffusionCalc()
{

}


double DegradedCSHDiffusionCalc::calcGelInclusionD(const double& gelPorosNitro,
                                                   const double& volCSH,
                                                   const double& volNonPerm,
                                                   const double& volPoros) const
{
    assert(gelPorosNitro > 0.0);

    const double freeWaterToPoreDiffRatio = 0.1;
    const double Dg = m_freeWaterD * freeWaterToPoreDiffRatio;
    const double gelDiff = m_gelDiffCalc->calc(Dg, gelPorosNitro);


    const double vol = volCSH + volNonPerm + volPoros;

    EffectiveMediaDiffusionCalc::InclusionsData inclusionData(2);
    EffectiveMediaDiffusionCalc::InclusionsVolFracDiffCoef& infoSolid = inclusionData[0];
    infoSolid.diffCoef = 0.0;
    infoSolid.volFrac = volNonPerm / vol;


    EffectiveMediaDiffusionCalc::InclusionsVolFracDiffCoef& infoPoros = inclusionData[1];

    const double porosFrac = volPoros / vol;
    infoPoros.diffCoef = m_freeWaterD;
    infoPoros.volFrac = porosFrac;
    const double matrixDiff = gelDiff;

    const double diffIncl = m_inclusionModel->calcDiffusionCoeffient(matrixDiff, inclusionData);

    assert(diffIncl >= 0.0);

    return diffIncl;
}

double DegradedCSHDiffusionCalc::calc(const SimplePorosityInfos &input) const
{
    const CSHPorosityInfos& concInput = static_cast<const CSHPorosityInfos&>(input);

    double diffCoef;

    if(input.porosityTotal > 1.0e-6)
    {
        const double& volSatCSH = concInput.volSatCSH;

        if(volSatCSH > 0.0)
        {
            const double gelPoreChange = degradPoreSlope * concInput.ratioCaSi
                                         + degradPoreIntercept;


            const double& volFracHD = concInput.volFracCSH_HD;
            const double satVolCSH_HD = volFracHD * volSatCSH;
            const double gelPoreHD = intrGelPorNitroHD + gelPoreChange;
            const double diff_HD = calcGelInclusionD(gelPoreHD, satVolCSH_HD,
                                                     concInput.volNonPerm * volFracHD,
                                                     0.0);



            // all capilary porosity is moved to LD-CSH
            const double volFracLD = 1.0 - volFracHD;
            const double satVolCSH_LD = volFracLD * volSatCSH;
            const double& porosityCapillary = concInput.porosityCapillary;
            const double gelPoreLD = intrGelPorNitroLD + gelPoreChange;
            const double volNonPermLD = concInput.volNonPerm * volFracLD;
            const double volCapillaryPoros = porosityCapillary * referenceVol;
            const double diff_LD = calcGelInclusionD(gelPoreLD, satVolCSH_LD,
                                                     volNonPermLD,
                                                     volCapillaryPoros);



            const double volFracLDAbs = (satVolCSH_LD + volCapillaryPoros + volNonPermLD) / referenceVol;
            const double volFracHDAbs = 1.0 - volFracLDAbs;


            EffectiveMediaDiffusionCalc::InclusionsData inclusionData(1);
            EffectiveMediaDiffusionCalc::InclusionsVolFracDiffCoef& infoMinor = inclusionData.front();

            // LD is the major matrix
            infoMinor.diffCoef = diff_HD;
            infoMinor.volFrac = volFracHDAbs;
            const double diffMatrix = diff_LD;
            diffCoef = m_inclusionModel->calcDiffusionCoeffient(diffMatrix, inclusionData);

        }
        else
        {
            const double volFracNonPerm = concInput.volNonPerm / referenceVol;
            const double porosityCapillary = 1.0 - volFracNonPerm;


            EffectiveMediaDiffusionCalc::InclusionsData inclusionData(1);
            EffectiveMediaDiffusionCalc::InclusionsVolFracDiffCoef& info = inclusionData.back();
            info.diffCoef = 0.0;
            info.volFrac = 1.0 - porosityCapillary;

            diffCoef = m_inclusionModel->calcDiffusionCoeffient(m_freeWaterD, inclusionData);
            IPPCheck::assertCheck(diffCoef >= 0.0, "got negative diff coef");
        }

    }
    else
    {
        // pure solid cells are not diffusive
        diffCoef = 0.0;
    }

    return diffCoef;

}




}
