#ifndef POROSITYBOUNCEBACK_H
#define POROSITYBOUNCEBACK_H

#include <palabos/core/dynamics.h>
#include <palabos/latticeBoltzmann/momentTemplates.h>
#include <palabos/core/dynamicsIdentifiers.h>

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class PorosityBounceBack : public plb::BounceBack<T,Descriptor>
{
public:
    PorosityBounceBack(T rho_=T())
        : plb::BounceBack<T,Descriptor>(rho_)
    {

    }

    PorosityBounceBack(plb::HierarchicUnserializer& unserializer)
    {
        this->unserialize(unserializer);
    }

    virtual PorosityBounceBack<T,Descriptor>* clone() const override
    {
        return new PorosityBounceBack<T,Descriptor>(*this);
    }

    void serialize(plb::HierarchicSerializer& serializer) const override
    {
        plb::BounceBack<T,Descriptor>::serialize(serializer);
    }

    void unserialize(plb::HierarchicUnserializer& unserializer) override
    {
        PLB_PRECONDITION( unserializer.getId() == this->getId() );
        plb::BounceBack<T,Descriptor>::unserialize(unserializer);
    }


    virtual int getId() const override
    {
        return id;
    }


    virtual T computeDensity(plb::Cell<T,Descriptor> const& cell) const override
    {
        const T rhoBarTmp = plb::momentTemplates<T,Descriptor>::get_rhoBar(cell);

        const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
        const T& poros = *pPoros;
        T rhoBar;

        if(poros == 0.0)
        {
            rhoBar = 0.0;
        }
        else
        {
            rhoBar = rhoBarTmp / poros;
        }

        return rhoBar;
    }

private:
    static int id;
};

template<typename T, template<typename U> class Descriptor>
int PorosityBounceBack<T,Descriptor>::id =
    plb::meta::registerGeneralDynamics<T,Descriptor,PorosityBounceBack<T,Descriptor> >("PorosityBounceBack");

}

#endif // POROSITYBOUNCEBACK_H
