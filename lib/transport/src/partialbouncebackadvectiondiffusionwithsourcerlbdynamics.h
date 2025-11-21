#ifndef PARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCERLBDYNAMICS_H
#define PARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCERLBDYNAMICS_H


#include <palabos/complexDynamics/advectionDiffusionDynamics.h>
#include <palabos/core/dynamicsIdentifiers.h>
#include "partialbouncebackcorrection.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class PartialBounceBackAdvectionDiffusionWithSourceRLBdynamics
        : public plb::AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>
{
public:
    typedef plb::AdvectionDiffusionWithSourceRLBdynamics<T,Descriptor> BaseClass;

    PartialBounceBackAdvectionDiffusionWithSourceRLBdynamics(T omega)
        : BaseClass(omega)
    {

    }

    PartialBounceBackAdvectionDiffusionWithSourceRLBdynamics(plb::HierarchicUnserializer& unserializer)
        : BaseClass(T())
    {
        this->unserialize(unserializer);
    }

    virtual PartialBounceBackAdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>* clone() const
    {
        return new PartialBounceBackAdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>(*this);
    }

    /// Return a unique ID for this class.
    virtual int getId() const
    {
        return id;
    }

    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics )
    {
        const plb::Array<T, Descriptor<T>::q> fOld = cell.getRawPopulations();
        BaseClass::collide(cell, statistics);
        PartialBounceBack::Walsh<T, Descriptor>::correct(fOld, cell);
    }

private:
    static int id;
};


template<typename T, template<typename U> class Descriptor>
int PartialBounceBackAdvectionDiffusionWithSourceRLBdynamics<T,Descriptor>::id =
        plb::meta::registerGeneralDynamics<T,Descriptor,PartialBounceBackAdvectionDiffusionWithSourceRLBdynamics<T,Descriptor> >
("PartialBounceBackAdvectionDiffusionWithSource_RLB");


}

#endif // PARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCERLBDYNAMICS_H
