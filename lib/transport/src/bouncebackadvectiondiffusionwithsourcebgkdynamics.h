#ifndef BOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H
#define BOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H

#include <palabos/complexDynamics/advectionDiffusionDynamics.hh>
#include <palabos/latticeBoltzmann/indexTemplates.h>

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class BounceBackAdvectionDiffusionWithSourceBGKdynamics
        : public plb::AdvectionDiffusionBGKdynamics<T,Descriptor>
{
public:
    BounceBackAdvectionDiffusionWithSourceBGKdynamics(T omega)
        : plb::AdvectionDiffusionBGKdynamics<T,Descriptor>(omega)
    {

    }

    BounceBackAdvectionDiffusionWithSourceBGKdynamics(plb::HierarchicUnserializer& unserializer)
        : plb::AdvectionDiffusionBGKdynamics<T,Descriptor>(T())
    {
        this->unserialize(unserializer);
    }

    virtual BounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>* clone() const
    {
        return new BounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>(*this);
    }

    /// Return a unique ID for this class.
    virtual int getId() const
    {
        return id;
    }

    virtual void collide(plb::Cell<T, Descriptor> &cell, plb::BlockStatistics & /*statistics*/)
    {
        plb::Array<T, Descriptor<T>::q>& f = cell.getRawPopulations();
        for (plb::plint iPop=1; iPop <= Descriptor<T>::q/2; ++iPop)
        {
            std::swap(cell[iPop], f[iPop+Descriptor<T>::q/2]);
        }
    }


private:
    static int id;
};


template<typename T, template<typename U> class Descriptor>
int BounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor>::id =
        plb::meta::registerGeneralDynamics<T,Descriptor,BounceBackAdvectionDiffusionWithSourceBGKdynamics<T,Descriptor> >
("BounceBackAdvectionDiffusionWithSource_BGK");


}

#endif // BOUNCEBACKADVECTIONDIFFUSIONWITHSOURCEBGKDYNAMICS_H
