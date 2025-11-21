#ifndef NOROUNDOFFD2Q5DESCRIPTORBASE_H
#define NOROUNDOFFD2Q5DESCRIPTORBASE_H

#include <palabos/latticeBoltzmann/advectionDiffusionLattices.h>

namespace IPP
{

template <typename T> struct NoRoundOffD2Q5DescriptorBase
    : public plb::descriptors::D2Q5Constants<T>, public plb::NoOptimizationRoundOffPolicy<T>
{
    typedef NoRoundOffD2Q5DescriptorBase<T> BaseDescriptor;
    enum { numPop=plb::descriptors::D2Q5Constants<T>::q };
};


}

#endif // NOROUNDOFFD2Q5DESCRIPTORBASE_H
