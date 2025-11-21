#ifndef POROSITYTRTDESCRIPTOR_H
#define POROSITYTRTDESCRIPTOR_H

#include "transportscalardescriptor.h"
#include "transportscalaradvectiondiffusiondescriptor.h"

namespace IPP
{

struct PorosityTRTDescriptor
{
    template<typename T>
    using HydrodynamicDescriptor2D = TransportScalarD2Q9Descriptor<T>;

    template<typename T>
    using HydrodynamicDescriptor3D = TransportScalarD3Q19Descriptor<T>;

    template<typename T>
    using DiffusionDescriptor2D = TransportScalarPorosityAdvectionDiffusionWithSourceD2Q5Descriptor<T>;

    template<typename T>
    using DiffusionDescriptor3D = TransportScalarPorosityAdvectionDiffusionWithSourceD3Q7Descriptor<T>;
};

}

#endif // POROSITYTRTDESCRIPTOR_H
