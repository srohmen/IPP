#ifndef POROSITYTRTDIFFUSIONOMEGACALC_H
#define POROSITYTRTDIFFUSIONOMEGACALC_H


#include "transportmoduleconfig.h"

#include "ginzburgtrtconstants.h"
#include "ippstream.h"

namespace IPP
{

namespace PorosityTRT
{

template<typename T>
T calcTransportScalar(const T& cPhi, const T& Dref, const T& relativeDiffCoef)
{
    const T D_LB = relativeDiffCoef * Dref;
    const T tau = (D_LB / cPhi) + 0.5;
    const T omega = 1.0 / tau;
    return omega;
}

template<typename T>
class OmegaMinusCalc
{
public:
    OmegaMinusCalc(const T& cPhi, const T& Dref)
        : m_cPhi(cPhi)
        , m_Dref(Dref)
    {

    }

    T operator()(const T& relativeDiffCoef) const
    {
        const T val = calcTransportScalar(m_cPhi, m_Dref, relativeDiffCoef);
        return val;
    }

private:
    const T m_cPhi;
    const T m_Dref;
};

}

template<typename T>
class PorosityTRTDiffusionOmegaCalc
{
public:
    PorosityTRTDiffusionOmegaCalc(const TransportModuleConfig& conf)
        : m_cPhi(GinzburgTRTConstants::calc_cPhi(conf.porosLow, conf.is3DSimulation))
        , m_tauRef(conf.tauRef)
    {
        IPP::pcout << "cPhi: " << m_cPhi << std::endl;
        initDrefViaPorosRef(conf.porosRef);
    }

    const T& getPorosRef() const
    {
        return m_porosRef;
    }

    void initDrefViaPorosRef(const T& porosRef)
    {
        m_Dref = (m_cPhi / porosRef) * (m_tauRef - 0.5);
        IPP::pcout << "Dref: " << m_Dref << std::endl;
        m_porosRef = porosRef;
    }

    T getInitDynamicsValue() const
    {
        return m_cPhi;
    }

    template<typename Descriptor>
    T calcTransportScalar(const T& /*omegaRef*/, const T& relativeDiffCoef) const
    {
        return PorosityTRT::calcTransportScalar(m_cPhi, m_Dref, relativeDiffCoef);
    }

    template<typename Descriptor>
    PorosityTRT::OmegaMinusCalc<T> makeTransportScalarCalcFunc() const
    {
        return PorosityTRT::OmegaMinusCalc<T>(m_cPhi, m_Dref);
    }

private:
    const T m_cPhi;
    const T m_tauRef;
    T m_porosRef;
    T m_Dref;
};


}

#endif // POROSITYTRTDIFFUSIONOMEGACALC_H
