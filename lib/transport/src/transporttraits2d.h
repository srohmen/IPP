#ifndef TRANSPORTTRAITS2D_H
#define TRANSPORTTRAITS2D_H

#include "derivedtransporttraits.h"

#include <palabos2D.h>
#include <palabos2D.hh>

namespace IPP
{

template<typename Scalar,
         template<typename U> class HydrodynamicDescriptorT,
         template<typename U> class DiffusionDescriptorT,
         template<typename U, template<typename V> class DescriptorT> class HydrDynT,
         template<typename U, template<typename V> class DescriptorT> class DiffDynT,
         typename HelperFunctions
         >
struct TransportTraits2D : public DerivedTransportTraits<
        Scalar,
        HydrodynamicDescriptorT,
        DiffusionDescriptorT,
        HydrDynT,
        DiffDynT,
        HelperFunctions,

        plb::MultiBlockLattice2D
        >
{
    static constexpr size_t dim = 2;
};


}

#endif // TRANSPORTTRAITS2D_H
