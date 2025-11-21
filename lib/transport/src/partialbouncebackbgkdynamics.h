#ifndef PARTIALBOUNCEBACKBGKDYNAMICS_H
#define PARTIALBOUNCEBACKBGKDYNAMICS_H

#include <palabos/basicDynamics/isoThermalDynamics.h>
#include <palabos/latticeBoltzmann/dynamicsTemplates.h>
#include <palabos/core/latticeStatistics.h>
#include <palabos/core/dynamicsIdentifiers.h>

#include "partialbouncebackcorrection.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class PartialBounceBackBGKDynamics : public plb::BGKdynamics<T,Descriptor>
{
public:
    typedef plb::BGKdynamics<T,Descriptor> BaseClass;

    PartialBounceBackBGKDynamics(T omega)
        : plb::BGKdynamics<T,Descriptor>(omega)
    {

    }

    PartialBounceBackBGKDynamics(plb::HierarchicUnserializer& unserializer)
        : plb::BGKdynamics<T,Descriptor>(T())
    {
        this->unserialize(unserializer);
    }


    /// Clone the object on its dynamic type.
    virtual PartialBounceBackBGKDynamics<T,Descriptor>* clone() const
    {
        return new PartialBounceBackBGKDynamics<T,Descriptor>(*this);
    }

    /// Return a unique ID for this class.
    virtual int getId() const
    {
        return id;
    }

    /// Implementation of the collision step
    virtual void collide(plb::Cell<T,Descriptor>& cell,
                         plb::BlockStatistics& statistics)
    {
        const plb::Array<T,Descriptor<T>::q> fOld = cell.getRawPopulations();
        BaseClass::collide(cell, statistics);
        PartialBounceBack::Walsh<T,Descriptor>::correct(fOld, cell);
    }

private:
    static int id;
};


}


template<typename T, template<typename U> class Descriptor>
int IPP::PartialBounceBackBGKDynamics<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,PartialBounceBackBGKDynamics<T,Descriptor> >("PartialBounceBackBGK");

#endif // PARTIALBOUNCEBACKBGKDYNAMICS_H
