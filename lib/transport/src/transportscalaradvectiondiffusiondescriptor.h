#ifndef TRANSPORTSCALARADVECTIONDIFFUSIONDESCRIPTOR_H
#define TRANSPORTSCALARADVECTIONDIFFUSIONDESCRIPTOR_H

#include <palabos/latticeBoltzmann/advectionDiffusionLattices.hh>
#include <palabos/latticeBoltzmann/nearestNeighborLattices2D.h>
#include "noroundoffd2q5descriptorbase.h"
#include "noroundoffd2q9descriptorbase.h"
#include "noroundoffd3q7descriptorbase.h"

namespace IPP
{

//////////////////////////////////////
// 2d part

struct TransportScalarAdvectionDiffusionDescriptor2D
{
    // V = velocity
    // S = source
    // T = transport
    //         0 1   2
    // layout: V V - T

    static const int numScalars = 3;
    static const int numSpecies = 2;
    static const int velocityBeginsAt = 0;
    static const int sizeOfVelocity   = 2;
    static const int transportScalarBeginsAt   = 2;
    static const int sizeOfTransportScalar = 1;
    static const int sizeOfForce      = 0;
};

struct TransportScalarAdvectionDiffusionWithSourceDescriptor2D
{
    // V = velocity
    // S = source
    // T = transport
    //         0 1   2   3
    // layout: V V - S - T

    static const int numScalars = 4;
    static const int numSpecies = 3;
    static const int velocityBeginsAt = 0;
    static const int sizeOfVelocity   = 2;
    static const int scalarBeginsAt   = 2;
    static const int sizeOfScalar     = 1;
    static const int transportScalarBeginsAt   = 3;
    static const int sizeOfTransportScalar = 1;
    static const int sizeOfForce      = 0;
};

struct TransportScalarPorosityAdvectionDiffusionWithSourceDescriptor2D
{
    // V = velocity
    // S = source
    // T = transport
    // P = porosity
    //         0 1   2   3   4
    // layout: V V - S - T - P

    static const int numScalars = 5;
    static const int numSpecies = 4;
    static const int velocityBeginsAt = 0;
    static const int sizeOfVelocity   = 2;
    static const int scalarBeginsAt   = 2;
    static const int sizeOfScalar     = 1;
    static const int transportScalarBeginsAt   = 3;
    static const int sizeOfTransportScalar = 1;
    static const int porosityBeginsAt = 4;
    static const int sizeOfPorosity = 1;
    static const int sizeOfForce      = 0;
};




struct TransportScalarAdvectionDiffusionDescriptor2DBase2D
{
    typedef TransportScalarAdvectionDiffusionDescriptor2D ExternalField;
};

struct TransportScalarAdvectionDiffusionWithSourceDescriptor2DBase2D
{
    typedef TransportScalarAdvectionDiffusionWithSourceDescriptor2D ExternalField;
};

struct TransportScalarPorosityAdvectionDiffusionWithSourceDescriptor2DBase2D
{
    typedef TransportScalarPorosityAdvectionDiffusionWithSourceDescriptor2D ExternalField;
};



/// AD D2Q5 lattice
template <typename T>
struct TransportScalarAdvectionDiffusionD2Q5Descriptor
        : public NoRoundOffD2Q5DescriptorBase<T>, public TransportScalarAdvectionDiffusionDescriptor2DBase2D
{
    static const char name[];
};

template<typename T>
const char TransportScalarAdvectionDiffusionD2Q5Descriptor<T>::name[] = "TransportScalarAdvectionDiffusionD2Q5";

/// ADS D2Q5 lattice
template <typename T>
struct TransportScalarAdvectionDiffusionWithSourceD2Q5Descriptor
        : public NoRoundOffD2Q5DescriptorBase<T>, public TransportScalarAdvectionDiffusionWithSourceDescriptor2DBase2D
{
    static const char name[];
};

template<typename T>
const char TransportScalarAdvectionDiffusionWithSourceD2Q5Descriptor<T>::name[] = "TransportScalarAdvectionDiffusionWithSourceD2Q5";

/// ADPS D2Q5 lattice
template <typename T>
struct TransportScalarPorosityAdvectionDiffusionWithSourceD2Q5Descriptor
        : public NoRoundOffD2Q5DescriptorBase<T>, public TransportScalarPorosityAdvectionDiffusionWithSourceDescriptor2DBase2D
{
    static const char name[];
};

template<typename T>
const char TransportScalarPorosityAdvectionDiffusionWithSourceD2Q5Descriptor<T>::name[] = "TransportScalarPorosityAdvectionDiffusionWithSourceD2Q5";


/// AD D2Q9 lattice
template <typename T>
struct TransportScalarAdvectionDiffusionD2Q9Descriptor
        : public NoRoundOffD2Q9DescriptorBase<T>, public TransportScalarAdvectionDiffusionDescriptor2DBase2D
{
    static const char name[];
};

template<typename T>
const char TransportScalarAdvectionDiffusionD2Q9Descriptor<T>::name[] = "TransportScalarAdvectionDiffusionD2Q9";


/// ADS D2Q9 lattice
template <typename T>
struct TransportScalarAdvectionDiffusionWithSourceD2Q9Descriptor
        : public NoRoundOffD2Q9DescriptorBase<T>, public TransportScalarAdvectionDiffusionWithSourceDescriptor2DBase2D
{
    static const char name[];
};

template<typename T>
const char TransportScalarAdvectionDiffusionWithSourceD2Q9Descriptor<T>::name[] = "TransportScalarAdvectionDiffusionWithSourceD2Q9";


//////////////////////////////////////
// 3d part

struct TransportScalarAdvectionDiffusionDescriptor3D
{
    // V = velocity
    // S = source
    // T = porosity
    //         0 1 2   3
    // layout: V V V - T

    static const int numScalars = 4;
    static const int numSpecies = 2;
    static const int velocityBeginsAt = 0;
    static const int sizeOfVelocity   = 3;
    static const int transportScalarBeginsAt   = 3;
    static const int sizeOfTransportScalar = 1;
    static const int sizeOfForce      = 0;
};

struct TransportScalarAdvectionDiffusionWithSourceDescriptor3D
{
    // V = velocity
    // S = source
    // T = porosity
    //         0 1 2   3   4
    // layout: V V V - S - T

    static const int numScalars = 5;
    static const int numSpecies = 3;
    static const int velocityBeginsAt = 0;
    static const int sizeOfVelocity   = 3;
    static const int scalarBeginsAt   = 3; // source scalar
    static const int sizeOfScalar     = 1;
    static const int transportScalarBeginsAt   = 4;
    static const int sizeOfTransportScalar = 1;
    static const int sizeOfForce      = 0;
};



struct TransportScalarPorosityAdvectionDiffusionWithSourceDescriptor3D
{
    // V = velocity
    // S = source
    // T = transport
    // P = porosity
    //         0 1 2   3   4   5
    // layout: V V V - S - T - P

    static const int numScalars = 6;
    static const int numSpecies = 4;
    static const int velocityBeginsAt = 0;
    static const int sizeOfVelocity   = 3;
    static const int scalarBeginsAt   = 3;
    static const int sizeOfScalar     = 1;
    static const int transportScalarBeginsAt   = 4;
    static const int sizeOfTransportScalar = 1;
    static const int porosityBeginsAt = 5;
    static const int sizeOfPorosity = 1;
    static const int sizeOfForce      = 0;
};



struct TransportScalarAdvectionDiffusionDescriptor3DBase3D
{
    typedef TransportScalarAdvectionDiffusionDescriptor3D ExternalField;
};


struct TransportScalarAdvectionDiffusionWithSourceDescriptor3DBase3D
{
    typedef TransportScalarAdvectionDiffusionWithSourceDescriptor3D ExternalField;
};

struct TransportScalarPorosityAdvectionDiffusionWithSourceDescriptor3DBase3D
{
    typedef TransportScalarPorosityAdvectionDiffusionWithSourceDescriptor3D ExternalField;
};


/// AD D3Q7 lattice
template <typename T>
struct TransportScalarAdvectionDiffusionD3Q7Descriptor
        : public NoRoundOffD3Q7DescriptorBase<T>, public TransportScalarAdvectionDiffusionDescriptor3DBase3D
{
    static const char name[];
};

template<typename T>
const char TransportScalarAdvectionDiffusionD3Q7Descriptor<T>::name[] = "TransportScalarAdvectionDiffusionD3Q7";

/// AD D3Q7 lattice
template <typename T>
struct TransportScalarAdvectionDiffusionWithSourceD3Q7Descriptor
        : public NoRoundOffD3Q7DescriptorBase<T>, public TransportScalarAdvectionDiffusionWithSourceDescriptor3DBase3D
{
    static const char name[];
};

template<typename T>
const char TransportScalarAdvectionDiffusionWithSourceD3Q7Descriptor<T>::name[] = "TransportScalarAdvectionDiffusionWithSourceD3Q7";




/// ADPS D3Q7 lattice
template <typename T>
struct TransportScalarPorosityAdvectionDiffusionWithSourceD3Q7Descriptor
        : public NoRoundOffD3Q7DescriptorBase<T>, public TransportScalarPorosityAdvectionDiffusionWithSourceDescriptor3DBase3D
{
    static const char name[];
};

template<typename T>
const char TransportScalarPorosityAdvectionDiffusionWithSourceD3Q7Descriptor<T>::name[] = "TransportScalarPorosityAdvectionDiffusionWithSourceD3Q7";




///// AD D3Q19 lattice
//template <typename T>
//struct TransportScalarAdvectionDiffusionD3Q19Descriptor
//        : public NoRoundOffD3Q19DescriptorBase<T>, public TransportScalarAdvectionDiffusionDescriptor3DBase3D
//{
//    static const char name[];
//};

//template<typename T>
//const char TransportScalarAdvectionDiffusionD3Q19Descriptor<T>::name[] = "TransportScalarAdvectionDiffusionD3Q19";


///// AD D3Q19 lattice
//template <typename T>
//struct TransportScalarAdvectionDiffusionWithSourceD3Q19Descriptor
//        : public NoRoundOffD3Q19DescriptorBase<T>, public TransportScalarAdvectionDiffusionWithSourceDescriptor3DBase3D
//{
//    static const char name[];
//};

//template<typename T>
//const char TransportScalarAdvectionDiffusionWithSourceD3Q19Descriptor<T>::name[] = "TransportScalarAdvectionDiffusionWithSourceD3Q19";

}

#endif // TRANSPORTSCALARADVECTIONDIFFUSIONDESCRIPTOR_H
