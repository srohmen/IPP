#ifndef PALABOSINITFUNCTORS_H
#define PALABOSINITFUNCTORS_H

#include <palabos/dataProcessors/dataInitializerFunctional2D.h>
#include <palabos/core/units.h>

#include "abstractflowfunctor.h"



namespace IPP
{


template<typename TransportTraits>
struct PalabosInitFunctors
{

/// A functional, used to instantiate bounce-back nodes at the locations of the cylinder
template <typename T>
class CylinderShapeDomain2D : public plb::DomainFunctional2D
{
public:
    CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(plb::util::sqr(radius))
    {

    }

    virtual bool operator() (plb::plint iX, plb::plint iY) const
    {
        return plb::util::sqr(iX-cx) + plb::util::sqr(iY-cy) <= radiusSqr;
    }

    virtual void operator() (plb::plint iX, plb::plint iY, T& rho, plb::Array<T, 2>& j) const
    {
        typedef plb::Array<T, 2> Vector2d;
        const Vector2d p(iX, iY);
        const Vector2d center(cx, cy);

        const Vector2d vec = p - center;
        const T length2 = vec[0] * vec[0] + vec[1] * vec[1];

        const T minConc = 0;
        const T maxConc = 1.0E-3;

        if (length2 <= radiusSqr)
        {
            const T interpolation = 1.0 - length2 / radiusSqr;
            PLB_ASSERT(interpolation >= 0.0);
            PLB_ASSERT(interpolation <= 1.0);

            rho = std::max(maxConc * interpolation, minConc);
        }
        else
        {
            rho = minConc;
        }

        j[0] = 0.0;
        j[1] = 0.0;
    }

    virtual CylinderShapeDomain2D<T>* clone() const {
        return new CylinderShapeDomain2D<T>(*this);
    }

private:
    plb::plint cx;
    plb::plint cy;
    plb::plint radiusSqr;
};

template<typename T>
class PoiseuilleVelocity
{
public:
    PoiseuilleVelocity(const plb::IncomprFlowParam<T>& parameters_)
        : parameters(parameters_)
    {
    }

    void operator()(plb::plint /*iX*/, plb::plint iY, plb::Array<T,2>& u) const
    {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }

    void operator()(plb::plint /*iX*/, plb::plint iY, plb::plint /*iZ*/, plb::Array<T,3>& u) const
    {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
        u[2] = T();
    }

private:
    const plb::IncomprFlowParam<T>& parameters;
};

/// A functional, used to initialize a pressure boundary to constant density
class ConstantDensity
{
public:
    using T = typename TransportTraits::Scalar;

    ConstantDensity(T density_)
        : density(density_)
    { }

    T operator()(plb::plint /*iX*/, plb::plint /*iY*/) const
    {
        return density;
    }

    T operator()(plb::plint /*iX*/, plb::plint /*iY*/, plb::plint /*iZ*/) const
    {
        return density;
    }

private:
    T density;
};

/// A functional, used to create an initial condition for the density and velocity
class PoiseuilleVelocityAndDensity
{
public:
    using T = typename TransportTraits::Scalar;

    PoiseuilleVelocityAndDensity(const plb::IncomprFlowParam<T>& parameters_)
        : parameters(parameters_)
    {
    }

    void operator()(plb::plint iX, plb::plint iY, T& rho, plb::Array<T,2>& u) const
    {
        rho = poiseuilleDensity(iX,parameters);

        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }

    void operator()(plb::plint iX, plb::plint iY, plb::plint /*iZ*/, T& rho, plb::Array<T,3>& u) const
    {
        rho = poiseuilleDensity(iX,parameters);

        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
        u[2] = T();
    }

private:
    const plb::IncomprFlowParam<T>& parameters;
};


class ConstantVelocity
{
public:
    using Vector = typename TransportTraits::Vector;

    ConstantVelocity(const Vector& u)
        : m_u(u)
    {

    }

    void operator()(plb::plint, plb::plint, Vector& u) const
    {
        for(size_t i = 0; i < TransportTraits::dim; ++i)
        {
            u[i] = m_u[i];
        }
    }

    void operator()(plb::plint, plb::plint, plb::plint, Vector& u) const
    {
        for(size_t i = 0; i < TransportTraits::dim; ++i)
        {
            u[i] = m_u[i];
        }
    }

private:
    const Vector& m_u;
};


template<typename T>
class ConstantVelocityAndDensity
{
public:
    typedef plb::Array<T,3> VectorType;
    ConstantVelocityAndDensity(const VectorType& u, const T& rho)
        : m_u(u),
          m_rho(rho)
    {

    }

    void operator()(plb::plint /*iX*/, plb::plint /*iY*/, T& rho, plb::Array<T,2>& u) const
    {
        rho = m_rho;

        u[0] = m_u[0];
        u[1] = m_u[1];
    }


    void operator()(plb::plint /*iX*/, plb::plint /*iY*/, plb::plint /*iZ*/, T& rho, plb::Array<T,3>& u) const
    {
        rho = m_rho;
        u = m_u;
    }

private:
    const VectorType m_u;
    const T m_rho;
};


class VelocityAndDensityWrapper
{
public:
    using T = typename TransportTraits::Scalar;

    VelocityAndDensityWrapper(const AbstractFlowFunctor& func)
        : m_func(func)
    {

    }

    void operator()(plb::plint iX, plb::plint iY, T& rho, plb::Array<T,2>& u) const
    {
        double rhoTmp;
        IPPVector2D uTmp;
        m_func(iX, iY, rhoTmp, uTmp);

        rho = rhoTmp;
        u[0] = uTmp[0];
        u[1] = uTmp[1];
    }

    void operator()(plb::plint iX, plb::plint iY, plb::plint iZ, T& rho, plb::Array<T,3>& u) const
    {
        double rhoTmp;
        IPPVector3D uTmp;
        m_func(iX, iY, iZ, rhoTmp, uTmp);

        rho = rhoTmp;
        u[0] = uTmp[0];
        u[1] = uTmp[1];
        u[2] = uTmp[2];
    }

private:
    const AbstractFlowFunctor& m_func;
};


template<typename T, template<class U> class Descriptor>
class SetOmegaDummyFunc : public plb::OneCellIndexedFunctional2D<T, Descriptor>
{
public:
    SetOmegaDummyFunc()
    {

    }

    virtual plb::OneCellIndexedFunctional2D<T,Descriptor>* clone() const
    {
        return new SetOmegaDummyFunc(*this);
    }

    virtual void execute(plb::plint iX, plb::plint iY, plb::Cell<T,Descriptor>& cell) const
    {
        plb::Dynamics<T,Descriptor>& dyn = cell.getDynamics();
        // FIXME
        dyn.setOmega(iX * 1000 + iY);
    }

private:

};

template<typename T, template<class U> class Descriptor>
class SetDiffCoeffFunction : public plb::OneCellIndexedFunctional2D<T, Descriptor>
{
public:
    SetDiffCoeffFunction(const size_t nx, const std::vector<double>& diffCoeffs)
        : nx(nx),
          diffCoeffs(diffCoeffs)
    {

    }

    virtual SetDiffCoeffFunction<T,Descriptor>* clone() const
    {
        return new SetDiffCoeffFunction(*this);
    }

    virtual void execute(plb::plint iX, plb::plint iY, plb::Cell<T,Descriptor>& cell) const
    {
        const size_t index = nx * iY + iX;
        const int diffCoeffOffset = Descriptor<T>::ExternalField::diffCoeffBeginsAt;
        double* diffCoeffFac = cell.getExternal(diffCoeffOffset);
        const double& diffCoeff = diffCoeffs.at(index);

        if(diffCoeff <= std::numeric_limits<double>::epsilon())
        {
            *diffCoeffFac = 0.0;
        }
    }

private:
    const size_t nx;
    const std::vector<double>& diffCoeffs;
};

template<typename T, template<class U> class Descriptor>
class PrintOmegaFunction : public plb::OneCellIndexedFunctional2D<T, Descriptor>
{
public:

    virtual plb::OneCellIndexedFunctional2D<T,Descriptor>* clone() const
    {
        return new PrintOmegaFunction(*this);
    }

    virtual void execute(plb::plint iX, plb::plint iY, plb::Cell<T,Descriptor>& cell) const
    {
        plb::pcout << iX << "\t" << iY << "\t" << cell.getDynamics().getOmega() << std::endl;
    }


private:
};

};

}

#endif // PALABOSINITFUNCTORS_H

