#ifndef CONFIGUREBOUNCEBACKDYNAMICS_H
#define CONFIGUREBOUNCEBACKDYNAMICS_H

#include <palabos/dataProcessors/dataInitializerWrapper2D.h>
#include <palabos/dataProcessors/dataInitializerWrapper3D.h>

#include "collectthreshcoordinates.h"
#include "collectnonpermeablecoordinates.h"
#include "instantiatebouncebackdynamic.h"


namespace IPP
{

template<typename TransportTraits>
struct ConfigureBounceBackDynamics
{
    template<typename T,
             typename DotList,
             typename ScalarField,
             template<class U> class Descriptor,
             template<class U, template<class V> class Desc> class LatticeType>
    static void apply(const DotList& newPermNodes, const DotList& newNonPermNodes,
                      ScalarField& concField, LatticeType<T, Descriptor>& lattice)
    {
        const plb::Dynamics<T, Descriptor>& backgroundDyn = lattice.getBackgroundDynamics();
        plb::defineDynamics(lattice, newPermNodes, backgroundDyn.clone());

        using InstantFunc = InstantiateBounceBackDynamic_LS<T, Descriptor>;
        plb::applyProcessingFunctional(new InstantFunc, newNonPermNodes, lattice, concField);
    }

};

}


EXTERN_TEMPLATE_TRANSPORT_TRAITS(struct IPP::ConfigureBounceBackDynamics)
#ifndef EXPLICIT_INSTANTS
#include "configurebouncebackdynamics.hh"
#endif


#endif // CONFIGUREBOUNCEBACKDYNAMICS_H
