#ifndef ADVECTIONDIFFUSIOND2Q9DESCRIPTOR_H
#define ADVECTIONDIFFUSIOND2Q9DESCRIPTOR_H

#include "noroundoffd2q9descriptorbase.h"


namespace IPP
{

// AD D2Q9 lattice
template <typename T>
struct AdvectionDiffusionD2Q9Descriptor
        : public NoRoundOffD2Q9DescriptorBase<T>, public plb::descriptors::Velocity2dDescriptorBase
{
    static const char name[];
};

template<typename T>
const char AdvectionDiffusionD2Q9Descriptor<T>::name[] = "AdvectionDiffusionD2Q9";


// AD D2Q9 lattice
template <typename T>
struct AdvectionDiffusionWithSourceD2Q9Descriptor
        : public NoRoundOffD2Q9DescriptorBase<T>, public plb::descriptors::VelocityAndScalar2dBase
{
    static const char name[];
};

template<typename T>
const char AdvectionDiffusionWithSourceD2Q9Descriptor<T>::name[] = "AdvectionDiffusionWithSourceD2Q9";

}

#endif // ADVECTIONDIFFUSIOND2Q9DESCRIPTOR_H
