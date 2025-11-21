#ifndef NOROUNDOFFD2Q9DESCRIPTORBASE_H
#define NOROUNDOFFD2Q9DESCRIPTORBASE_H

#include <palabos/latticeBoltzmann/nearestNeighborLattices2D.h>

namespace IPP
{

template <typename T> struct NoRoundOffD2Q9DescriptorBase
    : public plb::descriptors::D2Q9Constants<T>, public plb::NoOptimizationRoundOffPolicy<T>
{
    typedef NoRoundOffD2Q9DescriptorBase<T> BaseDescriptor;
    enum { numPop=plb::descriptors::D2Q9Constants<T>::q };
};


}

#endif // NOROUNDOFFD2Q9DESCRIPTORBASE_H
