#ifndef TRANSPORTSCALARDESCRIPTOR_H
#define TRANSPORTSCALARDESCRIPTOR_H

#include <palabos/latticeBoltzmann/nearestNeighborLattices2D.h>
#include <palabos/latticeBoltzmann/nearestNeighborLattices3D.h>

struct TransportScalarDescriptor2D
{
    static const int numScalars = 1;
    static const int numSpecies = 1;
    static const int transportScalarBeginsAt = 0;
    static const int sizeOfTransportScalar = 1;
    static const int sizeOfForce      = 0;
};

struct TransportScalarDescriptor2DBase2D
{
    typedef TransportScalarDescriptor2D ExternalField;
};

/// D2Q9 lattice
template <typename T>
struct TransportScalarD2Q9Descriptor
        : public plb::descriptors::D2Q9DescriptorBase<T>, public TransportScalarDescriptor2DBase2D
{
    static const char name[];
};

template<typename T>
const char TransportScalarD2Q9Descriptor<T>::name[] = "TransportScalarD2Q9Descriptor";



///////////////////
// 3D types

struct TransportScalarDescriptor3D
{
    static const int numScalars = 1;
    static const int numSpecies = 1;
    static const int transportScalarBeginsAt = 0;
    static const int sizeOfTransportScalar = 1;
    static const int sizeOfForce      = 0;
};

struct TransportScalarDescriptor3DBase3D
{
    typedef TransportScalarDescriptor3D ExternalField;
};

/// D2Q19 lattice
template <typename T>
struct TransportScalarD3Q19Descriptor
        : public plb::descriptors::D3Q19DescriptorBase<T>, public TransportScalarDescriptor3DBase3D
{
    static const char name[];
};

template<typename T>
const char TransportScalarD3Q19Descriptor<T>::name[] = "TransportScalarD3Q19Descriptor";



#endif // TRANSPORTSCALARDESCRIPTOR_H
