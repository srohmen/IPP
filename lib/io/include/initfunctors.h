#ifndef INITFUNCTORS_H
#define INITFUNCTORS_H

#include "abstractflowfunctor.h"
#include "ippvector.h"

namespace IPP
{


namespace InitFunctors
{

class ConstantVelocityAndDensity : public AbstractFlowFunctor
{
public:
    ConstantVelocityAndDensity(const IPPVector3D& u, const double& rho)
        : m_u(u),
          m_rho(rho)
    {

    }

    virtual bool isEnabled() const override
    {
        return true;
    }

    virtual void getDefault(double& rho, IPPVector3D& u) const override
    {
        rho = m_rho;
        u = m_u;
    }

    virtual void operator()(int iX, int iY, double& rho, IPPVector2D& u) const override
    {
        rho = m_rho;

        u[0] = m_u[0];
        u[1] = m_u[1];
    }

    virtual void operator()(int iX, int iY, int iZ, double& rho, IPPVector3D& u) const override
    {
        rho = m_rho;
        u = m_u;
    }

private:
    const IPPVector3D m_u;
    const double m_rho;

};

} // end of namespace InitFunctors

}

#endif // INITFUNCTORS_H
