#ifndef DIFFUSIONVELOCITYTRANSPORTTRAITS2D_H
#define DIFFUSIONVELOCITYTRANSPORTTRAITS2D_H

#include "transporttraits2d.h"

#include "transporttypes.h"
#include "diffusionvelocitytraits.h"
#include "diffusionvelocitydescriptor.h"

namespace IPP
{

struct DiffusionVelocityTransportTraits2D : public TransportTraits2D<
        TransportFloatingPointType,
        DiffusionVelocityDescriptor::HydrodynamicDescriptor2D,
        DiffusionVelocityDescriptor::DiffusionDescriptor2D,
        DiffusionVelocityTraits::HydrodynamicDynamics,
        DiffusionVelocityTraits::DiffusionDynamics,
        DiffusionVelocityTraits::HelperFunctions
        >
{
//    template<typename T>
//    using HydrodynamicDescriptorT = TransportScalarD2Q9Descriptor<T>;

//    template<typename T>
//    using DiffusionDescriptorT = TransportScalarAdvectionDiffusionWithSourceD2Q5Descriptor<T>;

};


}

#endif // DIFFUSIONVELOCITYTRANSPORTTRAITS2D_H
