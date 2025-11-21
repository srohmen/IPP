#ifndef ABSTRACTFLOWFUNCTOR_H
#define ABSTRACTFLOWFUNCTOR_H

#include <memory>
#include "ippvector.h"

namespace IPP
{

class AbstractFlowFunctor
{
public:
    AbstractFlowFunctor()
    {

    }

    virtual ~AbstractFlowFunctor()
    {

    }

    virtual bool isEnabled() const = 0;
    virtual void getDefault(double& rho, IPPVector3D& u) const = 0;
    virtual void operator()(int iX, int iY, double& rho, IPPVector2D& u) const = 0;
    virtual void operator()(int iX, int iY, int iZ, double& rho, IPPVector3D& u) const = 0;
};

typedef std::shared_ptr<AbstractFlowFunctor> AbstractFlowFunctorPtr;

}

#endif // ABSTRACTFLOWFUNCTOR_H
