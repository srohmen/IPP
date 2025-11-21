#ifndef POROSITYACCESS_H
#define POROSITYACCESS_H

#include <palabos/core/cell.h>

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
struct PorosityAccessT
{
    static T& get(plb::Cell<T, Descriptor>& cell)
    {
        T* cellVal = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
        return *cellVal;
    }

    static const T& get(const plb::Cell<T, Descriptor>& cell)
    {
        const T* cellVal = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
        return *cellVal;
    }

};


}

#endif // POROSITYACCESS_H
