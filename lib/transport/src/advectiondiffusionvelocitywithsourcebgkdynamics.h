#ifndef ADVECTIONDIFFUSIONVELOCITYWITHSOURCEBGKDYNAMICS_H
#define ADVECTIONDIFFUSIONVELOCITYWITHSOURCEBGKDYNAMICS_H

#include <palabos/complexDynamics/advectionDiffusionDynamics.h>
#include <palabos/latticeBoltzmann/momentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>
#include <palabos/core/dynamicsIdentifiers.h>

#include "diffusionvelocity.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionVelocityWithSourceBGKdynamics : public plb::AdvectionDiffusionWithSourceBGKdynamics<T, Descriptor>
{
public:
    AdvectionDiffusionVelocityWithSourceBGKdynamics(const T& omega)
        : plb::AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>(omega)
    {

    }

    AdvectionDiffusionVelocityWithSourceBGKdynamics(plb::HierarchicUnserializer& unserializer)
        : plb::AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>(T())
    {
        this->unserialize(unserializer);
    }


    virtual AdvectionDiffusionVelocityWithSourceBGKdynamics<T,Descriptor>* clone() const
    {
        return new AdvectionDiffusionVelocityWithSourceBGKdynamics<T,Descriptor>(*this);
    }


    /// Return a unique ID for this class.
    virtual int getId() const
    {
        return id;
    }

    /// Collision step
    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics )
    {
        const T omega = this->getOmega();

        T rhoBar;
        plb::Array<T,Descriptor<T>::d> jAdv;
        plb::Array<T,Descriptor<T>::d> jDiff;
        IPP::DiffusionVelocity::calcJAdvDiff(cell, omega, rhoBar, jAdv, jDiff);
        plb::Array<T,Descriptor<T>::d> jTot = jAdv + jDiff;

        const T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

        T uSqr = plb::advectionDiffusionDynamicsTemplates<T,Descriptor>::
                no_corr_bgk_collision(cell, rhoBar, jTot, omega, sourceTerm);

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
        // FIXME: this is not consistent with the other advection diffusion dynamics
        // why special case? check with external and internal (diffusion) velocity
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
int AdvectionDiffusionVelocityWithSourceBGKdynamics<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,AdvectionDiffusionVelocityWithSourceBGKdynamics<T,Descriptor> >
("AdvectionDiffusionVelocityWithSource_BGK");

}

#endif // ADVECTIONDIFFUSIONVELOCITYWITHSOURCEBGKDYNAMICS_H
