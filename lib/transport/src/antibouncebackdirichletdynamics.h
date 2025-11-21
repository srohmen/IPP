#ifndef ANTIBOUNCEBACKDIRICHLETDYNAMICS_H
#define ANTIBOUNCEBACKDIRICHLETDYNAMICS_H

#include <palabos/complexDynamics/advectionDiffusionDynamics.h>
#include <palabos/latticeBoltzmann/indexTemplates.h>
#include <palabos/latticeBoltzmann/momentTemplates.h>

#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>

namespace IPP
{


template<typename T, template<typename U> class Descriptor>
class AntiBounceBackDirichletDynamics : public plb::AdvectionDiffusionBGKdynamics<T,Descriptor>
{
public:
    AntiBounceBackDirichletDynamics(T omega_)
        : plb::AdvectionDiffusionBGKdynamics<T,Descriptor>(omega_)
    {

    }

    AntiBounceBackDirichletDynamics(plb::HierarchicUnserializer& unserializer)
        : plb::AdvectionDiffusionBGKdynamics<T,Descriptor>(T())
    {
        this->unserialize(unserializer);
    }

    virtual AntiBounceBackDirichletDynamics<T,Descriptor>* clone() const
    {
        return new AntiBounceBackDirichletDynamics<T,Descriptor>(*this);
    }

    virtual int getId() const
    {
        return id;
    }

    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics )
    {
        const T* pSrc = cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
        const T& source = *pSrc;


        plb::Array<T,Descriptor<T>::d> j;
        T rhoBar;
        plb::momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);

        const T* pU = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
        plb::Array<T,Descriptor<T>::d> jEq;
        jEq.from_cArray(pU);
        plb::advectionDiffusionDynamicsTemplates<T,Descriptor>::no_corr_bgk_collision(cell, rhoBar, jEq, this->getOmega());

        const T rhoBarBoundary = rhoBar + source;

        for (plb::plint i = 0; i <= Descriptor<T>::q/2; ++i)
        {
            const plb::plint iOpp = plb::indexTemplates::opposite<Descriptor<T>>(i);
            const T fi = cell[i];
            const T fOpp = cell[iOpp];

            cell[i] += 2.0 * Descriptor<T>::t[i] * rhoBarBoundary - fOpp;
            cell[iOpp] += 2.0 * Descriptor<T>::t[iOpp] * rhoBarBoundary - fi;
        }
    }

//    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& jEq,
//                                 T jSqr, T thetaBar=T()) const;

private:
    static int id;

};

template<typename T, template<typename U> class Descriptor>
int AntiBounceBackDirichletDynamics<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,AntiBounceBackDirichletDynamics<T,Descriptor> >("AntiBounceBackDirichlet");


}



#endif // ANTIBOUNCEBACKDIRICHLETDYNAMICS_H
