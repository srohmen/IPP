#ifndef ADVECTIONDIFFUSIONVELOCITYDESCRIPTOR_H
#define ADVECTIONDIFFUSIONVELOCITYDESCRIPTOR_H

#include <palabos/latticeBoltzmann/advectionDiffusionLattices.h>

struct VelocityDiffusionDescriptor2D
{
    static const int numScalars = 3;
    static const int numSpecies = 2;
    static const int velocityBeginsAt = 0;
    static const int sizeOfVelocity   = 2;
    static const int diffCoeffBeginsAt   = 2;
    static const int sizeOfDiffCoeff = 1;
    static const int sizeOfForce      = 0;
};

struct VelocityDiffusionDescriptor2DBase2D {
    typedef VelocityDiffusionDescriptor2D ExternalField;
};


namespace plb
{
namespace descriptors
{
/// AD D2Q5 lattice
template <typename T>
struct AdvectionDiffusionVelocityD2Q5Descriptor
        : public D2Q5DescriptorBase<T>, public VelocityDiffusionDescriptor2DBase2D
{
    static const char name[];
};

template<typename T>
const char AdvectionDiffusionVelocityD2Q5Descriptor<T>::name[] = "AdvectionDiffusionVelocityD2Q5";

}
}

#endif // ADVECTIONDIFFUSIONVELOCITYDESCRIPTOR_H
