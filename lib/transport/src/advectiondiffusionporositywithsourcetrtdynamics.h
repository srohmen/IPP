#ifndef ADVECTIONDIFFUSIONPOROSITYWITHSOURCETRTDYNAMICS_H
#define ADVECTIONDIFFUSIONPOROSITYWITHSOURCETRTDYNAMICS_H


#include <palabos/complexDynamics/trtDynamics.h>

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionPorosityWithSourceTRTDynamics : public plb::TRTdynamics<T, Descriptor>
{
public:
    AdvectionDiffusionPorosityWithSourceTRTDynamics(const T& cPhi);

    AdvectionDiffusionPorosityWithSourceTRTDynamics(plb::HierarchicUnserializer& unserializer);

    virtual AdvectionDiffusionPorosityWithSourceTRTDynamics<T,Descriptor>* clone() const;

    virtual int getId() const;


    virtual void serialize(plb::HierarchicSerializer& serializer) const;
    virtual void unserialize(plb::HierarchicUnserializer& unserializer);


    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(plb::Cell<T,Descriptor>& cell, T rhoBar,
                                 plb::Array<T,Descriptor<T>::d> const& j,
                                 T thetaBar, plb::BlockStatistics& stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plb::plint iPop, T rhoBar, plb::Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;


    virtual T computeDensity(plb::Cell<T,Descriptor> const& cell) const;


    virtual void computeRhoBarJ(plb::Cell<T,Descriptor> const& cell,
                                T& rhoBar, plb::Array<T,Descriptor<T>::d>& j ) const;

private:
    static const T s_magic;
    static int id;

    T m_cPhi;

};

} // end of namespace IPP

#ifdef EXPLICIT_INSTANTS
#include "externtemplatedescriptor.h"
EXTERN_TEMPLATE_SCALAR_DIFFPOROSITYDESCRIPTOR(class IPP::AdvectionDiffusionPorosityWithSourceTRTDynamics)
#else
#include "advectiondiffusionporositywithsourcetrtdynamics.hh"
#endif

#endif // ADVECTIONDIFFUSIONPOROSITYWITHSOURCETRTDYNAMICS_H
