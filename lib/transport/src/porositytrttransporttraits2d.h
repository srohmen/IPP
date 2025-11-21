#ifndef POROSITYTRTTRANSPORTTRAITS2D_H
#define POROSITYTRTTRANSPORTTRAITS2D_H

#include "transporttraits2d.h"

#include "transporttypes.h"
#include "porositytrttraits.h"
#include "porositytrtdescriptor.h"

namespace IPP
{

struct PorosityTRTTransportTraits2D : public TransportTraits2D<
        TransportFloatingPointType,
        PorosityTRTDescriptor::HydrodynamicDescriptor2D,
        PorosityTRTDescriptor::DiffusionDescriptor2D,
        PorosityTRTTraits::HydrodynamicDynamics,
        PorosityTRTTraits::DiffusionDynamics,
        PorosityTRTTraits::HelperFunctions
        >
{
    // TODO: for some reason the descriptors must be on top level of classes. find solution to put that into derived traits!

    //    template<typename T>
    //    using HydrodynamicDescriptorT = TransportScalarD2Q9Descriptor<T>;

    //    template<typename T>
    //    using DiffusionDescriptorT = TransportScalarPorosityAdvectionDiffusionWithSourceD2Q5Descriptor<T>;
};

}

#endif // POROSITYTRTTRANSPORTTRAITS2D_H
