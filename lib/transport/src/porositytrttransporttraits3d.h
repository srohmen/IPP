#ifndef POROSITYTRTTRANSPORTTRAITS3D_H
#define POROSITYTRTTRANSPORTTRAITS3D_H

#include "transporttraits3d.h"

#include "transporttypes.h"
#include "porositytrttraits.h"
#include "porositytrtdescriptor.h"

namespace IPP
{

struct PorosityTRTTransportTraits3D : public TransportTraits3D<
        TransportFloatingPointType,
        PorosityTRTDescriptor::HydrodynamicDescriptor3D,
        PorosityTRTDescriptor::DiffusionDescriptor3D,
        PorosityTRTTraits::HydrodynamicDynamics,
        PorosityTRTTraits::DiffusionDynamics,
        PorosityTRTTraits::HelperFunctions
        >
{
    //    template<typename T>
    //    using HydrodynamicDescriptorT = TransportScalarD3Q19Descriptor<T>;

    //    template<typename T>
    //    using DiffusionDescriptorT = TransportScalarPorosityAdvectionDiffusionWithSourceD3Q7Descriptor<T>;
};

}

#endif // POROSITYTRTTRANSPORTTRAITS3D_H
