#ifndef NOROUNDOFFD3Q7DESCRIPTORBASE_H
#define NOROUNDOFFD3Q7DESCRIPTORBASE_H

#include <palabos/latticeBoltzmann/advectionDiffusionLattices.h>

namespace IPP
{

template <typename T> struct NoRoundOffD3Q7DescriptorBase
        : public plb::descriptors::D3Q7Constants<T>, public plb::NoOptimizationRoundOffPolicy<T>
{
    typedef NoRoundOffD3Q7DescriptorBase<T> BaseDescriptor;
    enum { numPop=plb::descriptors::D3Q7Constants<T>::q };
};

}

#endif // NOROUNDOFFD3Q7DESCRIPTORBASE_H
