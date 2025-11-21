#ifndef DISTANCETRANSFORM8SSEDT_2D_H
#define DISTANCETRANSFORM8SSEDT_2D_H

#include "distancefieldcalc.h"

namespace IPP
{

class DistanceTransform8SSEDT_2D : public DistanceFieldCalc
{
public:
    DistanceTransform8SSEDT_2D();
    virtual void calc(const AbstractScalarField& input, const double& threshold, AbstractScalarField& output);
};

}

#endif // DISTANCETRANSFORM8SSEDT_2D_H
