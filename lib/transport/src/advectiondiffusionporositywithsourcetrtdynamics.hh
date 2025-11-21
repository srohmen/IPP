#ifndef ADVECTIONDIFFUSIONPOROSITYWITHSOURCETRTDYNAMICS_HH
#define ADVECTIONDIFFUSIONPOROSITYWITHSOURCETRTDYNAMICS_HH

#include "advectiondiffusionporositywithsourcetrtdynamics.h"

#include <palabos/core/dynamicsIdentifiers.h>
#include <palabos/core/latticeStatistics.h>

#include <palabos/latticeBoltzmann/dynamicsTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>

#include <palabos/latticeBoltzmann/momentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>

#include <palabos/latticeBoltzmann/indexTemplates.h>

#include "ginzburgtrt.h"


namespace IPP
{

template<typename T, template<typename U> class Descriptor>
const T AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::s_magic = 1.0/4.0;


///////// ordinary palabos housekeeping stuff /////////
template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::id =
        plb::meta::registerGeneralDynamics<T,Descriptor,AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor> >("ADSPorosity_TRT");

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>::AdvectionDiffusionPorosityWithSourceTRTDynamics(const T& cPhi)
    : plb::TRTdynamics<T, Descriptor>(1.0)
    , m_cPhi(cPhi)
{

}

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>::AdvectionDiffusionPorosityWithSourceTRTDynamics(plb::HierarchicUnserializer& unserializer)
    : plb::TRTdynamics<T, Descriptor>(T())
    , m_cPhi(-1.0)
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>* AdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>::clone() const
{
    return new AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>::getId() const
{
    return id;
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::serialize(plb::HierarchicSerializer& serializer) const
{
    plb::TRTdynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(m_cPhi);
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::unserialize(plb::HierarchicUnserializer& unserializer)
{
    plb::TRTdynamics<T, Descriptor>::unserialize(unserializer);
    m_cPhi = unserializer.readValue<T>();
}

//////////////////////////////////////////////////////



////// here comes the algorithmic implementation:

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::collide(plb::Cell<T, Descriptor> &cell,
                                                                            plb::BlockStatistics &statistics)
{
    const T* pTrans = cell.getExternal(Descriptor<T>::ExternalField::transportScalarBeginsAt);
    const T& sMinus = *pTrans;
    const T tauMinus = 1.0 / sMinus;
    const T tauPlus = s_magic / (tauMinus - 0.5) + 0.5;
    const T sPlus = 1.0 / tauPlus;


    const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
    const T& poros = *pPoros;

    const T* pSrc = cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
    const T& absSource = *pSrc;

    const T total = plb::momentTemplates<T,Descriptor>::get_rhoBar(cell);
    const T rhoBar = total / poros;

    const T* pVelocity = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
    plb::Array<T,Descriptor<T>::d> jEq;
    jEq.from_cArray(pVelocity);
    T jSqr = normSqr(jEq);

    plb::Array<T,Descriptor<T>::q>& f = cell.getRawPopulations();

    plb::Array<T,Descriptor<T>::q> eq;
    GinzburgTRT::equilibrium<T, Descriptor<T>>(m_cPhi, rhoBar, jEq, jSqr, poros, eq);
//    plb::dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria(rhoBar, 1.0, j, jSqr, eq);


    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    plb::Array<T,Descriptor<T>::q/2+1> eq_plus, eq_minus, f_plus, f_minus;

    for (plb::plint i=1; i<=Descriptor<T>::q/2; ++i)
    {
        const plb::plint iOpp = i + Descriptor<T>::q / 2;
        eq_plus[i]  = 0.5*(eq[i] + eq[iOpp]);
        eq_minus[i] = 0.5*(eq[i] - eq[iOpp]);
        f_plus[i]   = 0.5*(f[i] + f[iOpp]);
        f_minus[i]  = 0.5*(f[i] - f[iOpp]);
    }


    const T w = 0.5;
    const T srcOther = w * m_cPhi * (absSource / poros);
    const T src0 = absSource - 2 * Descriptor<T>::d * srcOther;


    f[0] += sPlus * (eq[0] - f[0]) + src0;

    for (plb::plint i = 1; i <= Descriptor<T>::q/2; ++i)
    {
        const plb::plint iOpp = i + Descriptor<T>::q/2;

        const T termPlus    = -sPlus   * ( f_plus[i]    - eq_plus[i]    ) + srcOther;
        const T termMinus   = -sMinus  * ( f_minus[i]   - eq_minus[i]   );
        f[i]     += termPlus + termMinus;
        f[iOpp]  += termPlus - termMinus;
    }


    if (cell.takesStatistics())
    {        
        const T invRho = Descriptor<T>::invRho(rhoBar);
        plb::gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho );
    }
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::collideExternal(plb::Cell<T, Descriptor> &cell, T rhoBar,
                                                                                    const plb::Array<T, Descriptor<T>::d> &j,
                                                                                    T thetaBar, plb::BlockStatistics &stat)
{
    PLB_ASSERT(false);
}

template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::computeEquilibrium(plb::plint iPop, T rhoBar,
                                                                                    const plb::Array<T, Descriptor<T>::d> &j,
                                                                                    T jSqr, T /*thetaBar*/) const
{
    // porosity is not available in this function
    // assumuming rhoBar and j is macroscopic concentration and flux (porosity = 1)
    // and NOT in TRT scales with respect to porosity
    //
    // !!!
    // equilibrium value of rest velocity f[0] must be
    // corrected by porosity after function return
    // !!!

    if(iPop != 0)
    {
        const T f = GinzburgTRT::equilibrium<T, Descriptor<T>>(iPop, m_cPhi, rhoBar, j, jSqr, 1.0);
        return f;
    }
    else
    {
        plb::Array<T, Descriptor<T>::q> eq;
        GinzburgTRT::equilibrium<T, Descriptor<T>>(m_cPhi, rhoBar, j, jSqr, 1.0, eq);
        return eq[0];
    }
}

template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::computeDensity(plb::Cell<T,Descriptor> const& cell) const
{
    const T rhoBarTmp = plb::momentTemplates<T,Descriptor>::get_rhoBar(cell);

    const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
    const T& poros = *pPoros;

    assert(std::isnan(rhoBarTmp) == false);
    const T rhoBar = rhoBarTmp / poros;

    // poros and rhoBar could be zero at the same time
    // -> yields NaN
    // assert(std::isnan(rhoBar) == false);


    return rhoBar;
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::computeRhoBarJ(plb::Cell<T,Descriptor> const& cell,
                                                                                   T& rhoBar,
                                                                                   plb::Array<T,Descriptor<T>::d>& j ) const
{
    plb::momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);

    const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
    const T& poros = *pPoros;


    rhoBar /= poros;

    for(plb::plint i = 0; i < Descriptor<T>::d; ++i)
    {
        j[i] /= poros;
    }
}

}


#endif // ADVECTIONDIFFUSIONPOROSITYWITHSOURCETRTDYNAMICS_HH


