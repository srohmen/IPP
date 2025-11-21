#include "cshgeldiffusioncoefcalc.h"

#include <cmath>

namespace IPP
{

CSHGelDiffusionCoefCalc::CSHGelDiffusionCoefCalc()
{

}

double CSHGelDiffusionCoefCalc::calc(const double& Dgel, const double &gelPoros)
{
    // 1.5 corresponds to spherical inclusion
    const double n = 1.5;
    const double tau = std::pow(gelPoros, n);
    const double D = Dgel * tau;
    return D;
}

}
