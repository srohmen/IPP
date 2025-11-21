#ifndef LATTICEBOLTZMANNUNITCONVERSION_H
#define LATTICEBOLTZMANNUNITCONVERSION_H

#include <cstddef>

namespace IPP
{

class IPPConfig;

namespace LatticeBoltzmannUnitConversion
{

double calcTimeStepSRT(const double& dx, const double& D, const double& tau);
double calcTimeStepPTRT(const double& dx, const double& D, const double& tau,
                        const double& porosLow, const double porosRef, const bool is3D);
double calcTimeStep(const IPPConfig& conf);

}

}

#endif // LATTICEBOLTZMANNUNITCONVERSION_H
