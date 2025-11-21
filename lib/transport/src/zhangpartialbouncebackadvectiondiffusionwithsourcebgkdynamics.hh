#include "zhangpartialbouncebackadvectiondiffusionwithsourcebgkdynamics.h"

#include "partialbouncebackcorrection.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
int ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor> >
("ZhangPartialBounceBackAdvectionDiffusionWithSource_BGK");

template<typename T, template<typename U> class Descriptor>
ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::
ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics(T omega)
    : plb::AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>(omega)
{

}

template<typename T, template<typename U> class Descriptor>
ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::
ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics(plb::HierarchicUnserializer& unserializer)
    : plb::AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>*
ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::clone() const
{
    return new ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::getId() const
{
    return id;
}

// FIXME: this implementation is incomplete

const double ns = 0.0;

template<typename T, typename Descriptor>
T equilibriumZhang(const plb::plint iPop,
              const T& rhoBar,
              const plb::Array<T,Descriptor::d>& jTmp,
              const T& ns)
{
    T jd = 0;
    for(plb::plint d = 0; d < Descriptor::d; ++d)
    {
        jd += Descriptor::c[iPop][d] * jTmp[d];
    }

    const T poroRhoBar = (1.0 - ns) * rhoBar;

    const T eq = (poroRhoBar + 3.5 * jd) / 7.0;

    return eq;

}

template<typename T, typename Descriptor>
void calcJTmp(const T& omega, const T& ns,
              const plb::Array<T,Descriptor::d>& j,
              plb::Array<T,Descriptor::d>& jTmp)
{
    const T tau = 1.0 / omega;
    const T tmp = tau / (tau - 2.0 * ns * ((1 - 0.5 * tau) + ns * ( tau - 1.0)) );
    jTmp = j / tmp;
    // TODO: if ns == 1 -> jTmp = 0
}

template<typename T, typename Descriptor>
void equilibria(const T rhoBar,
                const plb::Array<T,Descriptor::d>& j,
                const T& omega, const T& ns,
                plb::Array<T,Descriptor::q>& eq)
{
    plb::Array<T,Descriptor::d> jTmp;
    calcJTmp<T, Descriptor>(omega, ns, j, jTmp);

    for(plb::plint iPop = 1; iPop < Descriptor::q; ++iPop)
    {
        eq[iPop] = equilibriumZhang<T, Descriptor>(iPop, rhoBar, jTmp, ns);
    }

}


template<typename T, template<typename U> class Descriptor>
void ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::collide(
        plb::Cell<T,Descriptor>& cell,
        plb::BlockStatistics& statistics )
{
    const plb::Array<T,Descriptor<T>::q> fOld = cell.getRawPopulations();

    plb::Array<T,Descriptor<T>::d> jDiff;
    T rhoBar;
    plb::momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, jDiff);

    T rhoBarScaled = rhoBar / (1.0 - ns);

    plb::Array<T,Descriptor<T>::d> jAdv(0, 0);

    const T omega = this->getOmega();

    plb::Array<T,Descriptor<T>::q> eq;
    equilibria<T,Descriptor<T>>(rhoBarScaled, jAdv, omega, ns, eq);


    for (plb::plint iPop=0; iPop < Descriptor<T>::q; ++iPop)
    {
        cell[iPop] *= (T)1 - omega;
        cell[iPop] += omega * eq[iPop];
    }

    plb::Array<T,Descriptor<T>::q>& fNew = cell.getRawPopulations();

    PartialBounceBack::Walsh<T,Descriptor>::correct(ns, fOld, fNew);
}


template<typename T, template<typename U> class Descriptor>
T ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::computeEquilibrium(
        plb::plint iPop, T rhoBar, plb::Array<T,Descriptor<T>::d> const& j,
        T /*jSqr*/, T /*thetaBar*/) const
{
    // TODO: find solution for define ns

    const T omega = this->getOmega();

    T rhoBarScaled = rhoBar / (1.0 - ns);
    plb::Array<T,Descriptor<T>::d> jTmp;
    calcJTmp<T, Descriptor<T>>(omega, ns, j, jTmp);

    const T eq = equilibriumZhang<T,Descriptor<T>>(iPop, rhoBarScaled, jTmp, ns);
    return eq;
}

}
