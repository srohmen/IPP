#ifndef PARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H
#define PARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H

#include <palabos/complexDynamics/advectionDiffusionDynamics.hh>
#include <palabos/core/dynamicsIdentifiers.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>


#include "partialbouncebackcorrection.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class PartialBounceBackAdvectionDiffusionWithSourceBGKdynamics
        : public plb::AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>
{
public:
    typedef plb::AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor> BaseClass;

    PartialBounceBackAdvectionDiffusionWithSourceBGKdynamics(T omega)
        : BaseClass(omega)
    {

    }    

    PartialBounceBackAdvectionDiffusionWithSourceBGKdynamics(plb::HierarchicUnserializer& unserializer)
        : BaseClass(T())
    {
        this->unserialize(unserializer);
    }

    virtual PartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>* clone() const
    {
        return new PartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>(*this);
    }

    /// Return a unique ID for this class.
    virtual int getId() const
    {
        return id;
    }

    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics )
    {
        const plb::Array<T,Descriptor<T>::q> fOld = cell.getRawPopulations();
        BaseClass::collide(cell, statistics);
        PartialBounceBack::Walsh<T,Descriptor>::correct(fOld, cell);
    }


private:
    static int id;
};


template<typename T, template<typename U> class Descriptor>
int PartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,PartialBounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor> >
("PartialBounceBackAdvectionDiffusionWithSource_BGK");


}

#endif // PARTIALBOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H
