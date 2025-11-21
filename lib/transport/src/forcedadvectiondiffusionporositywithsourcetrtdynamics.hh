#include "forcedadvectiondiffusionporositywithsourcetrtdynamics.h"

#include <palabos/core/dynamicsIdentifiers.h>

#include <palabos/latticeBoltzmann/dynamicsTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>

#include <palabos/latticeBoltzmann/momentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>

#include <palabos/latticeBoltzmann/indexTemplates.h>

#include "ginzburgtrt.h"



namespace IPP
{

template<typename T, template<typename U> class Descriptor>
const T ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::s_magic = 1.0/4.0;


///////// ordinary palabos housekeeping stuff /////////
template<typename T, template<typename U> class Descriptor>
int ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::id =
        plb::meta::registerGeneralDynamics<T,Descriptor,
ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor> >("ForcedADSPorosity_TRT");

template<typename T, template<typename U> class Descriptor>
ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>::
ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics(const T &omega, const T& minPoros)
    : plb::TRTdynamics<T, Descriptor>(omega)
    , m_cPhi(minPoros / 3.0) // TODO: make general for 3D lattices as well
{

}

template<typename T, template<typename U> class Descriptor>
ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>::
ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics(plb::HierarchicUnserializer& unserializer)
    : plb::TRTdynamics<T, Descriptor>(T())
    , m_cPhi(-1.0)
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>*
ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>::clone() const
{
    return new ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>::getId() const
{
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::
serialize(plb::HierarchicSerializer& serializer) const
{
    plb::TRTdynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(m_cPhi);
}

template<typename T, template<typename U> class Descriptor>
void ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::
unserialize(plb::HierarchicUnserializer& unserializer)
{
    plb::TRTdynamics<T, Descriptor>::unserialize(unserializer);
    m_cPhi = unserializer.readValue<T>();
}

//////////////////////////////////////////////////////



////// here comes the algorithmic implementation:

template<typename T, template<typename U> class Descriptor>
void ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::
collide(plb::Cell<T, Descriptor> &cell,
        plb::BlockStatistics &statistics)
{
    const T* pTrans = cell.getExternal(Descriptor<T>::ExternalField::transportScalarBeginsAt);
    const T& tauMinus = *pTrans;
    const T sMinus = 1.0 / tauMinus;
    const T tauPlus = s_magic / (tauMinus - 0.5) + 0.5;
    const T sPlus = 1.0 / tauPlus;

    const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
    const T& poros = *pPoros;

    const T* pSrc = cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
    const T& source = *pSrc;

    plb::Array<T,Descriptor<T>::d> jInit;
    T rhoBar;
    plb::momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, jInit);
    rhoBar /= poros;

    const T* pVelocity = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
    plb::Array<T,Descriptor<T>::d> u;
    u.from_cArray(pVelocity);
    const plb::Array<T,Descriptor<T>::d> jEq = u * rhoBar;
    const T jSqr = plb::normSqr(jEq);

    plb::Array<T,Descriptor<T>::q>& f = cell.getRawPopulations();
    plb::Array<T, Descriptor<T>::q> eq;

    for(plb::plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
    {
        const T fEq = plb::advectionDiffusionDynamicsTemplates<T,Descriptor>::
                bgk_ma1_equilibrium(iPop, rhoBar, jEq);
        eq[iPop] = fEq;
    }

    plb::Array<T,Descriptor<T>::q> fSource;
    fSource.resetToZero();

    const T length = plb::norm(jInit);

    if(length != 0.0)
    {
        if(source != 0)
        {
            std::cout << jInit << std::endl;
        }

        plb::Array<T,Descriptor<T>::d> jNorm;
        jNorm = jInit / length;

        const plb::Array<T,Descriptor<T>::d> jSource = std::abs(source) * jNorm;

        for (plb::plint i = 0; i < Descriptor<T>::q; ++i)
        {
            T& sourceTerm = fSource[i];
            for(plb::plint d = 0; d < Descriptor<T>::d; ++d)
            {
                sourceTerm += jSource[d] * Descriptor<T>::c[i][d];
            }
            sourceTerm *= Descriptor<T>::t[i] * Descriptor<T>::invCs2;
        }


        for (plb::plint i = 0; i < Descriptor<T>::q; ++i)
        {
            T& sourceTerm = fSource[i];
            if(source > 0.0 && sourceTerm < 0.0)
            {
                const plb::plint iOpp = plb::indexTemplates::opposite<Descriptor<T>>(i);
                fSource[iOpp] -= sourceTerm;
                sourceTerm = 0.0;
            }
            else if(source < 0.0 && sourceTerm > 0.0)
            {
                const plb::plint iOpp = plb::indexTemplates::opposite<Descriptor<T>>(i);
                fSource[iOpp] -= sourceTerm;
                sourceTerm = 0.0;
            }
        }


        plb::Array<T,Descriptor<T>::d> jZero;
        jZero.resetToZero();

        // relax directional source according to relaxation parameter later
        for (plb::plint iPop=0; iPop < Descriptor<T>::q; ++iPop)
        {
            const T t = sMinus - 1;
            const T fEq = plb::advectionDiffusionDynamicsTemplates<T,Descriptor>:: bgk_ma1_equilibrium(iPop, source, jZero);

            fSource[iPop] *= t;
            fSource[iPop] += (1.0 - t) * fEq;

        }
    }
    else
    {
        plb::Array<T,Descriptor<T>::d> jZero;
        jZero.resetToZero();

        for (plb::plint iPop=0; iPop < Descriptor<T>::q; ++iPop)
        {
            const T fEq = plb::advectionDiffusionDynamicsTemplates<T,Descriptor>:: bgk_ma1_equilibrium(iPop, source, jZero);
            fSource[iPop] = fEq;
        }
    }

    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    plb::Array<T,Descriptor<T>::q> eq_plus, eq_minus, f_plus, f_minus;

    for (plb::plint i = 0; i < Descriptor<T>::q; ++i)
    {
        const plb::plint iOpp = plb::indexTemplates::opposite<Descriptor<T>>(i);
        eq_plus[i]  = 0.5 * (eq[i] + eq[iOpp]);
        eq_minus[i] = 0.5 * (eq[i] - eq[iOpp]);
        f_plus[i]   = 0.5 * (f[i] + f[iOpp]);
        f_minus[i]  = 0.5 * (f[i] - f[iOpp]);
    }

    T sumPlus = 0;
    T sumMinus = 0;

    T sumfPlus = 0;
    T sumfMinus = 0;
    T sumF = 0;

    for (plb::plint i = 0; i < Descriptor<T>::q; ++i)
    {
        const T termPlus    = sPlus   * ( f_plus[i]    - eq_plus[i]    );
        const T termMinus   = sMinus  * ( f_minus[i]   - eq_minus[i]   );

        sumF += f[i];
        sumfPlus += f_plus[i];
        sumfMinus += f_minus[i];
        sumPlus += termPlus;
        sumMinus += termMinus;

        f[i] += -termPlus - termMinus;
    }

    f += fSource;

    if (cell.takesStatistics())
    {
        const T invRho = Descriptor<T>::invRho(rhoBar);
        const T uSqr = jSqr * invRho * invRho;
        plb::gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::
collideExternal(plb::Cell<T, Descriptor> &cell, T rhoBar,
                const plb::Array<T, Descriptor<T>::d> &j,
                T thetaBar, plb::BlockStatistics &stat)
{
    PLB_ASSERT(false);
}

template<typename T, template<typename U> class Descriptor>
T ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::
computeEquilibrium(plb::plint iPop, T rhoBar,
                   const plb::Array<T, Descriptor<T>::d> &j,
                   T jSqr, T /*thetaBar*/) const
{
    // assumuming rhoBar and j is already in TRT scales with respect to porosity

    if(iPop != 0)
    {
        const T f = GinzburgTRT::equilibrium<T, Descriptor<T>>(iPop, m_cPhi, rhoBar,
                                                               j, jSqr, 1.0);
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
T ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::
computeDensity(plb::Cell<T,Descriptor> const& cell) const
{
    T rhoBar = plb::momentTemplates<T,Descriptor>::get_rhoBar(cell);

    const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
    const T& poros = *pPoros;
    rhoBar /= poros;

    return rhoBar;
}

template<typename T, template<typename U> class Descriptor>
void ForcedAdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>::
computeRhoBarJ(plb::Cell<T,Descriptor> const& cell,
               T& rhoBar,
               plb::Array<T,Descriptor<T>::d>& j ) const
{
    plb::momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);

    const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
    const T& poros = *pPoros;
    rhoBar /= poros;

    for(plb::plint i = 0; i <= Descriptor<T>::d; ++i)
    {
        j[i] /= poros;
    }
}

}

