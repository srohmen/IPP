#ifndef THRESHOLDMULTISCALEDIFFUSIONCALC_H
#define THRESHOLDMULTISCALEDIFFUSIONCALC_H

#include "abstractmultiscalediffusioncalc.h"

namespace IPP
{

class ThresholdMultiScaleDiffusionCalc : public AbstractMultiScaleDiffusionCalc
{
public:
    ThresholdMultiScaleDiffusionCalc(const double& thresh, const double reference);


    virtual void setPhaseNames(const std::vector<std::string>* names);

    virtual double calc(const SimplePorosityInfos &input) const override;

private:
    const double m_thresh;
    const double m_reference;
};

}

#endif // THRESHOLDMULTISCALEDIFFUSIONCALC_H
