#ifndef DIFFUSIONVELOCITY_H
#define DIFFUSIONVELOCITY_H

#include <palabos/latticeBoltzmann/momentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>

#include "palabosio.h"

namespace IPP
{

namespace DiffusionVelocity
{

template<typename T, typename Descriptor>
T calcRefDiff(const T& omega)
{
    const T refDiff = Descriptor::cs2 * (1.0 / omega - 0.5);
    return refDiff;
}

template<typename T, typename Descriptor>
T calcDiffVelFactor(const T& omega, const T& accelFac)
{
    const T refDiff = calcRefDiff<T, Descriptor>(omega);
    const T diffDifference = refDiff * (accelFac  - 1.0);

    const T fac = omega * Descriptor::invCs2;
    const T tmp = diffDifference * fac;

    const T diffVelFac = tmp / (1.0 + tmp);

    return diffVelFac;
}

template<typename T>
T calcDiffTmp(const T& diffVelFac)
{
    const T tmp = -diffVelFac / (diffVelFac - 1.0);
    return tmp;
}

template<typename T, typename Descriptor>
class CalcDiffVelFactor
{
public:
    CalcDiffVelFactor(const T& omega)
        : omega(omega)
    {

    }

    T operator()(const T& accelFac) const
    {
        return calcDiffVelFactor<T, Descriptor>(omega, accelFac);
    }

private:
    const T omega;
};

template<typename T, template<typename U> class DescriptorT>
void calcJAdvDiff(const plb::Cell<T,DescriptorT>& cell, const T& omega, T& rhoBar,
                  plb::Array<T,DescriptorT<T>::d>& jAdv, plb::Array<T,DescriptorT<T>::d>& jDiff)
{
    typedef DescriptorT<T> Descriptor;

    plb::Array<T,Descriptor::d> j;
    plb::momentTemplates<T, DescriptorT>::get_rhoBar_j(cell, rhoBar, j);

    plb::advectionDiffusionMomentTemplates<T,DescriptorT>::get_jEq(cell, rhoBar, jAdv);

    const T diffVelFac = *cell.getExternal(Descriptor::ExternalField::transportScalarBeginsAt);

    const plb::Array<T,Descriptor::d> jDifference = j - jAdv;
    jDiff = diffVelFac * jDifference;
}

template<typename T, template<typename U> class DescriptorT>
void get_rhoBar_j(plb::Cell<T,DescriptorT> const& cell, const T& omega,
                  T& rhoBar, plb::Array<T,DescriptorT<T>::d>& j)
{
    typedef DescriptorT<T> Descriptor;

    plb::Array<T,Descriptor::d> jCurr;
    plb::momentTemplates<T,DescriptorT>::get_rhoBar_j(cell, rhoBar, jCurr);

    plb::Array<T,Descriptor::d> jAdv;
    plb::advectionDiffusionMomentTemplates<T,DescriptorT>::get_jEq(cell, rhoBar, jAdv);

    const T diffVelFac = *cell.getExternal(Descriptor::ExternalField::transportScalarBeginsAt);
    const plb::Array<T,Descriptor::d> jDifference = jCurr - jAdv;
    const plb::Array<T,Descriptor::d> jDiff = jDifference * diffVelFac;

    j = jCurr + jDiff;
}

template<typename T, template<typename U> class DescriptorT>
void compute_uLb(plb::Cell<T,DescriptorT> const& cell, const T& omega,
                 plb::Array<T,DescriptorT<T>::d>& u)
{
    T rhoBar;
    plb::Array<T,DescriptorT<T>::d> j;
    DiffusionVelocity::get_rhoBar_j(cell, omega, rhoBar, j);
    u = j * DescriptorT<T>::invRho(rhoBar);
}


}
} // end of namespace IPP

#endif // DIFFUSIONVELOCITY_H
