#ifndef CONSTANTRHOADVECTIONDIFFUSIONBGKDYNAMICS_H
#define CONSTANTRHOADVECTIONDIFFUSIONBGKDYNAMICS_H

#include <palabos/boundaryCondition/boundaryDynamics.h>

#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>
#include <palabos/complexDynamics/advectionDiffusionDynamics.h>

#include <palabos/core/dynamicsIdentifiers.h>
#include <palabos/core/latticeStatistics.h>


namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class ConstantRhoAdvectionDiffusionBGKdynamics : public plb::StoreDensityDynamics<T,Descriptor>
{
public:

    typedef plb::StoreDensityDynamics<T,Descriptor> BaseClass;
    typedef plb::AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor> DelegateClass;

    ConstantRhoAdvectionDiffusionBGKdynamics(T omega)
        : BaseClass(new DelegateClass(omega))
    {

    }

    ConstantRhoAdvectionDiffusionBGKdynamics(plb::HierarchicUnserializer& unserializer)
        : BaseClass(new DelegateClass(T()))
    {
          this->unserialize(unserializer);
    }

    virtual ConstantRhoAdvectionDiffusionBGKdynamics<T,Descriptor>* clone() const
    {
        return new ConstantRhoAdvectionDiffusionBGKdynamics<T,Descriptor>(*this);
    }

    /// Return a unique ID for this class.
    virtual int getId() const
    {
        return id;
    }

    virtual void prepareCollision(plb::Cell<T,Descriptor>& cell)
    {
        const T* velArr = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
        plb::Array<T,Descriptor<T>::d> j;
        j.from_cArray(velArr);
        const T jSqr = normSqr(j);

        plb::Array<T,Descriptor<T>::q> fEq;
        this->computeEquilibria(fEq, this->rhoBar, j, jSqr);

        const plb::Array<T,Descriptor<T>::q> fOld = cell.getRawPopulations();

        const plb::Array<T,Descriptor<T>::q> fNeq = fOld - fEq;

        const T omega = this->getOmega();
        const T oneMinusOmega = 1.0 - omega;
        const plb::Array<T,Descriptor<T>::q> fNew = fEq + oneMinusOmega * fNeq;
        cell.setPopulations(fNew);
    }

private:
    static int id;
};


template<typename T, template<typename U> class Descriptor>
int ConstantRhoAdvectionDiffusionBGKdynamics<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,ConstantRhoAdvectionDiffusionBGKdynamics<T,Descriptor> >
("ConstantRhoAdvectionDiffusion_BGK");


} // end of namespace LBGeoChem

#endif // CONSTANTRHOADVECTIONDIFFUSIONBGKDYNAMICS_H
