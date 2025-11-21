#include "thresholdmultiscalediffusioncalc.h"

#include <cassert>

#include "simpleporosityinfos.h"

namespace IPP
{

ThresholdMultiScaleDiffusionCalc::ThresholdMultiScaleDiffusionCalc(const double &thresh, const double reference)
    : m_thresh(thresh),
      m_reference(reference)
{

}

void ThresholdMultiScaleDiffusionCalc::setPhaseNames(const std::vector<std::string>* /*names*/)
{

}

double ThresholdMultiScaleDiffusionCalc::calc(const SimplePorosityInfos &input) const
{    
    const double currPoros = 1.0 - input.porosityTotal;

    if(currPoros < m_thresh)
    {
        return 0.0;
    }
    else
    {
         return m_reference;
    }

}


}
