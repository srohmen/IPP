#ifndef ARCHIESLAWDIFFUSIONCALC_H
#define ARCHIESLAWDIFFUSIONCALC_H

#include "abstractmultiscalediffusioncalc.h"

namespace IPP
{

class ArchiesLawDiffusionCalc : public AbstractMultiScaleDiffusionCalc
{
public:
    ArchiesLawDiffusionCalc(const double& D0, const double& exponent);
    virtual ~ArchiesLawDiffusionCalc();

    virtual double calc(const SimplePorosityInfos& input) const override;

private:
    const double m_D0;
    const double m_expo;

};

}

#endif // ARCHIESLAWDIFFUSIONCALC_H
