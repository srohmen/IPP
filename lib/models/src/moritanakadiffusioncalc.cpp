#include "moritanakadiffusioncalc.h"

#include <cassert>

namespace IPP
{

static double calcBeta(const double& matrixDiff, const double& inclusionDiff )
{
    const double beta = (inclusionDiff - matrixDiff) / (inclusionDiff + 2.0 * matrixDiff);
    return beta;
}


MoriTanakaDiffusionCalc::MoriTanakaDiffusionCalc()
{

}

double MoriTanakaDiffusionCalc::calcDiffusionCoeffient(const double& matrixDiff,
                                                       const InclusionsData &inclusionData) const
{
    assert(matrixDiff > 0.0);

    double sum = 0.0;

    for(const InclusionsVolFracDiffCoef& incInfo : inclusionData)
    {
        const double& inclusionDiff = incInfo.diffCoef;
        const double& volFrac = incInfo.volFrac;
        const double beta = calcBeta(matrixDiff, inclusionDiff);
        const double term = beta * volFrac;
        sum += term;
    }

    const double diffResult = matrixDiff * (1.0 + 2.0 * sum) / (1.0 - sum);

    return diffResult;
}

}
