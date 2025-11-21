#ifndef DIFFUSIONVELOCITYDIFFUSIONOMEGACALC_H
#define DIFFUSIONVELOCITYDIFFUSIONOMEGACALC_H

#include "diffusionvelocity.h"

#include "transportmoduleconfig.h"

namespace IPP
{


template<typename T>
class DiffusionVelocityDiffusionOmegaCalc
{
public:
    DiffusionVelocityDiffusionOmegaCalc(TransportModuleConfig &conf)
        : m_omegaRef(1.0 / conf.tauRef)
    {

    }


    T getInitDynamicsValue() const
    {
        return m_omegaRef;
    }

    template<typename Descriptor>
    T calcTransportScalar(const T& omegaRef, const T& accelFac) const
    {
        const T val = DiffusionVelocity::calcDiffVelFactor<T, Descriptor>(omegaRef, accelFac);
        return val;
    }

    template<typename Descriptor>
    DiffusionVelocity::CalcDiffVelFactor<T, Descriptor> makeTransportScalarCalcFunc() const
    {
        return DiffusionVelocity::CalcDiffVelFactor<T, Descriptor>(m_omegaRef);
    }

private:
    const T m_omegaRef;

};

}

#endif // DIFFUSIONVELOCITYDIFFUSIONOMEGACALC_H
