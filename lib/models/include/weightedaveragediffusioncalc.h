#ifndef WEIGHTEDAVERAGEDIFFUSIONCALC_H
#define WEIGHTEDAVERAGEDIFFUSIONCALC_H

#include "abstractweightingfunc.h"

namespace IPP
{

class WeightedAverageDiffusionCalc : public AbstractWeightingFunc
{
public:
    WeightedAverageDiffusionCalc();
    virtual ~WeightedAverageDiffusionCalc();

    virtual double calcDiffusionCoeffient(const double& a, const double& b, const double& t) const;
};

}

#endif // WEIGHTEDAVERAGEDIFFUSIONCALC_H
