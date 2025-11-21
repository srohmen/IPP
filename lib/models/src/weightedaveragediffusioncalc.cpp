#include "weightedaveragediffusioncalc.h"


namespace IPP
{


WeightedAverageDiffusionCalc::WeightedAverageDiffusionCalc()
{

}

WeightedAverageDiffusionCalc::~WeightedAverageDiffusionCalc()
{

}

double WeightedAverageDiffusionCalc::calcDiffusionCoeffient(const double &a,
                                                            const double &b,
                                                            const double &t) const
{
    const double result = a * t + b * (1.0 - t);
    return result;
}

}
