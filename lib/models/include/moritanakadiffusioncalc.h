#ifndef MORITANAKADIFFUSIONCALC_H
#define MORITANAKADIFFUSIONCALC_H

#include "effectivemediadiffusioncalc.h"

namespace IPP
{

class MoriTanakaDiffusionCalc : public EffectiveMediaDiffusionCalc
{
public:
    MoriTanakaDiffusionCalc();


    virtual double calcDiffusionCoeffient(const double &matrixDiff,
                                          const InclusionsData &inclusionData) const;
};

}

#endif // MORITANAKADIFFUSIONCALC_H
