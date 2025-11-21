#ifndef PARTIALBOUNCEBACK_H
#define PARTIALBOUNCEBACK_H

#include <palabos/core/cell.h>
#include <palabos/latticeBoltzmann/geometricOperationTemplates.h>
#include <palabos/latticeBoltzmann/indexTemplates.h>

namespace IPP
{

namespace PartialBounceBack
{

template<typename T, template<typename U> class Descriptor>
struct Walsh
{    
    static void correct(const T& ns,
                        const plb::Array<T,Descriptor<T>::q>& preCollision,
                        plb::Array<T,Descriptor<T>::q>& fNew)
    {
        for (plb::plint iPop=1; iPop <= Descriptor<T>::q/2; ++iPop)
        {
            const plb::plint iOpp = plb::indexTemplates::opposite<Descriptor<T>>(iPop);
            fNew[iPop] += ns * (preCollision[iOpp] - fNew[iPop]);
            fNew[iOpp] += ns * (preCollision[iPop] - fNew[iOpp]);
        }

    }

    static void correct(const plb::Array<T,Descriptor<T>::q>& preCollision,
                        plb::Cell<T,Descriptor>& cell)
    {
        const T* pTransportScalar = cell.getExternal(Descriptor<T>::ExternalField::transportScalarBeginsAt);
        const T transportScalar = *pTransportScalar ;

        plb::Array<T,Descriptor<T>::q>& fNew = cell.getRawPopulations();

        correct(transportScalar, preCollision, fNew);
    }

};
}


}

#endif // PARTIALBOUNCEBACK_H
