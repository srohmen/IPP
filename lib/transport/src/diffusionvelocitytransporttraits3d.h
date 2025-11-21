#ifndef DIFFUSIONVELOCITYTRANSPORTTRAITS3D_H
#define DIFFUSIONVELOCITYTRANSPORTTRAITS3D_H

#include "transporttraits3d.h"

#include "transporttypes.h"
#include "diffusionvelocitytraits.h"
#include "diffusionvelocitydescriptor.h"

namespace IPP
{

//typedef TransportTraits3D
//<
//double,
//TransportScalarD3Q19Descriptor,
//TransportScalarAdvectionDiffusionWithSourceD3Q7Descriptor
//>
//DiffusionVelocityTransportTraits3D;

struct DiffusionVelocityTransportTraits3D : public TransportTraits3D<
        TransportFloatingPointType,
        DiffusionVelocityDescriptor::HydrodynamicDescriptor3D,
        DiffusionVelocityDescriptor::DiffusionDescriptor3D,
        DiffusionVelocityTraits::HydrodynamicDynamics,
        DiffusionVelocityTraits::DiffusionDynamics,
        DiffusionVelocityTraits::HelperFunctions
        >
{
//    template<typename T>
//    using HydrodynamicDescriptorT = TransportScalarD3Q19Descriptor<T>;

//    template<typename T>
//    using DiffusionDescriptorT = TransportScalarAdvectionDiffusionWithSourceD3Q7Descriptor<T>;

};

}


#endif // DIFFUSIONVELOCITYTRANSPORTTRAITS3D_H
