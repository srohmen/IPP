#ifndef ADVECTIONDIFFUSIONVELOCITYWITHSOURCERLBDYNAMICS_H
#define ADVECTIONDIFFUSIONVELOCITYWITHSOURCERLBDYNAMICS_H


#include <palabos/complexDynamics/advectionDiffusionDynamics.h>
#include <palabos/latticeBoltzmann/momentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>
#include <palabos/core/dynamicsIdentifiers.h>

#include "diffusionvelocity.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionVelocityWithSourceRLBdynamics : public plb::AdvectionDiffusionWithSourceRLBdynamics<T, Descriptor>
{
public:
    AdvectionDiffusionVelocityWithSourceRLBdynamics(const T& omega)
        : plb::AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>(omega)
    {

    }


    virtual AdvectionDiffusionVelocityWithSourceRLBdynamics<T,Descriptor>* clone() const
    {
        return new AdvectionDiffusionVelocityWithSourceRLBdynamics<T,Descriptor>(*this);
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
        plb::Array<T,Descriptor<T>::d> jTot;
        DiffusionVelocity::calcJAdvDiff(cell, omega, rhoBar, jTot);

        const T sourceTerm = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

        T uSqr = plb::advectionDiffusionDynamicsTemplates<T,Descriptor>::
                no_corr_bgk_collision(cell, rhoBar, jTot, omega, sourceTerm);

        if (cell.takesStatistics())
        {
            gatherStatistics(statistics, rhoBar, uSqr);
        }


    }



private:
    static int id;

};


template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionVelocityWithSourceRLBdynamics<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,AdvectionDiffusionVelocityWithSourceRLBdynamics<T,Descriptor> >
("AdvectionDiffusionVelocityWithSource_BGK");

}

#endif // ADVECTIONDIFFUSIONVELOCITYWITHSOURCERLBDYNAMICS_H
