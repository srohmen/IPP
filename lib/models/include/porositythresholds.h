#ifndef POROSITYTHRESHOLDS_H
#define POROSITYTHRESHOLDS_H

namespace IPP
{

struct PorosityThresholds
{
    const double porosLower;
    const double porosUpper;

    const double neighEdgeVolFrac;
    const double neighCornerVolFrac;

};

}

#endif // POROSITYTHRESHOLDS_H
