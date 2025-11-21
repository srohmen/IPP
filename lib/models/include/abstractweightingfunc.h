#ifndef ABSTRACTWEIGHTINGFUNC_H
#define ABSTRACTWEIGHTINGFUNC_H

namespace IPP
{

class AbstractWeightingFunc
{
public:
    virtual ~AbstractWeightingFunc()
    {

    }

    virtual double calcDiffusionCoeffient(const double& a, const double& b, const double& t) const = 0;
};

}


#endif // ABSTRACTWEIGHTINGFUNC_H
