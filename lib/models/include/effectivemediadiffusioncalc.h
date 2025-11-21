#ifndef EFFECTIVEMEDIADIFFUSIONCALC_H
#define EFFECTIVEMEDIADIFFUSIONCALC_H

#include <vector>

namespace IPP
{

class EffectiveMediaDiffusionCalc
{
public:

    struct InclusionsVolFracDiffCoef
    {
        double volFrac;
        double diffCoef;
    };

    typedef std::vector<InclusionsVolFracDiffCoef> InclusionsData;

    EffectiveMediaDiffusionCalc()
    {

    }

    virtual ~EffectiveMediaDiffusionCalc()
    {

    }


    virtual double calcDiffusionCoeffient(const double &matrixDiff,
                                          const InclusionsData &inclusionData) const = 0;
};

}

#endif // EFFECTIVEMEDIADIFFUSIONCALC_H
