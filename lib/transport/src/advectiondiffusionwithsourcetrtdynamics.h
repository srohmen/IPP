#ifndef ADVECTIONDIFFUSIONWITHSOURCETRTDYNAMICS_H
#define ADVECTIONDIFFUSIONWITHSOURCETRTDYNAMICS_H

#include <palabos/complexDynamics/trtDynamics.h>

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionWithSourceTRTDynamics : public plb::TRTdynamics<T, Descriptor>
{
public:
    AdvectionDiffusionWithSourceTRTDynamics(const T& omega);

    virtual AdvectionDiffusionWithSourceTRTDynamics<T,Descriptor>* clone() const;

    virtual int getId() const;

    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(plb::Cell<T,Descriptor>& cell, T rhoBar,
                                 plb::Array<T,Descriptor<T>::d> const& j,
                                 T thetaBar, plb::BlockStatistics& stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plb::plint iPop, T rhoBar, plb::Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;




private:
    static const T s_magic;
    static int id;

};


}

#endif // ADVECTIONDIFFUSIONWITHSOURCETRTDYNAMICS_H
