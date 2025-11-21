#ifndef DISTANCEFIELDCALC_H
#define DISTANCEFIELDCALC_H


namespace IPP
{
class AbstractScalarField;

class DistanceFieldCalc
{
public:
    DistanceFieldCalc()
    {

    }

    virtual ~DistanceFieldCalc()
    {

    }

    virtual void calc(const AbstractScalarField& input, const double& threshold, AbstractScalarField& output) = 0;

};

}

#endif // DISTANCEFIELDCALC_H
