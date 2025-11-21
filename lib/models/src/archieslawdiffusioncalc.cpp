#include "archieslawdiffusioncalc.h"

#include <cmath>
#include <cassert>
#include "simpleporosityinfos.h"

namespace IPP
{

ArchiesLawDiffusionCalc::ArchiesLawDiffusionCalc(const double& D0,
                                                 const double &exponent)
    : m_D0(D0),
      m_expo(exponent)
{

}

ArchiesLawDiffusionCalc::~ArchiesLawDiffusionCalc()
{

}

double ArchiesLawDiffusionCalc::calc(const SimplePorosityInfos& input) const
{
    const double& porosity = input.porosityTotal;
    double De = 0.0;

    // let De be zero for "negative" porosities
    if(porosity > 0.0)
    {
        const double Dp = m_D0 * std::pow(porosity, m_expo);
        De = Dp * porosity;
    }


    assert(std::isnan(De) == false);
    assert(De >= 0.0);

    return De;

}


}
