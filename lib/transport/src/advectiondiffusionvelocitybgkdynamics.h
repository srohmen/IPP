#ifndef ADVECTIONDIFFUSIONVELOCITYBGKDYNAMICS_H
#define ADVECTIONDIFFUSIONVELOCITYBGKDYNAMICS_H

#include <palabos/complexDynamics/advectionDiffusionDynamics.h>
#include <palabos/core/dynamicsIdentifiers.h>

#include "diffusionvelocity.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionVelocityBGKdynamics : public plb::AdvectionDiffusionBGKdynamics <T,Descriptor>
{
public:
    AdvectionDiffusionVelocityBGKdynamics(const T& omega)
        : plb::AdvectionDiffusionBGKdynamics <T,Descriptor>(omega)
    {

    }

    AdvectionDiffusionVelocityBGKdynamics(plb::HierarchicUnserializer& unserializer)
        : plb::AdvectionDiffusionBGKdynamics <T,Descriptor>(T())
    {
        this->unserialize(unserializer);
    }

    virtual AdvectionDiffusionVelocityBGKdynamics<T,Descriptor>* clone() const
    {
        return new AdvectionDiffusionVelocityBGKdynamics<T,Descriptor>(*this);
    }


    /// Return a unique ID for this class.
    virtual int getId() const
    {
        return id;
    }

    /// Collision step
    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics)
    {
        const T omega = this->getOmega();

        T rhoBar;
        plb::Array<T,Descriptor<T>::d> jAdv;
        plb::Array<T,Descriptor<T>::d> jDiff;
        DiffusionVelocity::calcJAdvDiff(cell, omega, rhoBar, jAdv, jDiff);
        plb::Array<T,Descriptor<T>::d> jTot = jAdv + jDiff;

        T uSqr = plb::advectionDiffusionDynamicsTemplates<T,Descriptor>::
                no_corr_bgk_collision(cell, rhoBar, jTot, omega);

        const T jSqr = plb::VectorTemplate<T,Descriptor>::normSqr(jAdv);
        const T invRho = Descriptor<T>::invRho(rhoBar);
        uSqr = jSqr*invRho*invRho;

        if (cell.takesStatistics())
        {
            gatherStatistics(statistics, rhoBar, uSqr);
        }
    }

    virtual void collideExternal (
            plb::Cell<T,Descriptor>& cell, T rhoBar,
            plb::Array<T,Descriptor<T>::d> const& j, T thetaBar, plb::BlockStatistics& stat )
    {
        PLB_ASSERT(false);
    }

    virtual void computeVelocity( plb::Cell<T,Descriptor> const& cell,
                                  plb::Array<T,Descriptor<T>::d>& u ) const
    {
        DiffusionVelocity::compute_uLb(cell, this->getOmega(), u);
    }

    virtual void computeRhoBarJ(plb::Cell<T,Descriptor> const& cell,
                                T& rhoBar, plb::Array<T,Descriptor<T>::d>& j ) const
    {
        DiffusionVelocity::get_rhoBar_j(cell, this->getOmega(), rhoBar, j);
    }


private:
    static int id;

};


template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionVelocityBGKdynamics<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,AdvectionDiffusionVelocityBGKdynamics<T,Descriptor> >
("AdvectionDiffusionVelocity_BGK");


}


#endif // ADVECTIONDIFFUSIONVELOCITYBGKDYNAMICS_H
