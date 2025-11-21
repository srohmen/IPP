#ifndef DIFFUSIONVELOCITYDESCRIPTOR_H
#define DIFFUSIONVELOCITYDESCRIPTOR_H


#include "transportscalardescriptor.h"
#include "transportscalaradvectiondiffusiondescriptor.h"

namespace IPP
{

struct DiffusionVelocityDescriptor
{
    template<typename T>
    using HydrodynamicDescriptor2D = TransportScalarD2Q9Descriptor<T>;

    template<typename T>
    using HydrodynamicDescriptor3D = TransportScalarD3Q19Descriptor<T>;

    template<typename T>
    using DiffusionDescriptor2D = TransportScalarAdvectionDiffusionWithSourceD2Q5Descriptor<T>;

    template<typename T>
    using DiffusionDescriptor3D = TransportScalarAdvectionDiffusionWithSourceD3Q7Descriptor<T>;
};

}

#endif // DIFFUSIONVELOCITYDESCRIPTOR_H
