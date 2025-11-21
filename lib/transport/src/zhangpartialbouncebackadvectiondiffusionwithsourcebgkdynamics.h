#ifndef ZHANGPARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H
#define ZHANGPARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H

#include <palabos/complexDynamics/advectionDiffusionDynamics.hh>
#include <palabos/core/dynamicsIdentifiers.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>

namespace IPP
{


template<typename T, template<typename U> class Descriptor>
class ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics
        : public plb::AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>
{
public:
    ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics(T omega);

    ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics(plb::HierarchicUnserializer& unserializer);

    virtual ZhangPartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>* clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Collision step
    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics );
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plb::plint iPop, T rhoBar, plb::Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;


private:
    static int id;
};

}

#endif // ZHANGPARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H
