#ifndef TRANSPORTTRAITS3D_H
#define TRANSPORTTRAITS3D_H

#include <palabos3D.h>
#include <palabos3D.hh>

#include <palabos/atomicBlock/dataProcessingFunctional3D.h>
#include <palabos/atomicBlock/reductiveDataProcessingFunctional3D.h>
#include <palabos/dataProcessors/dataInitializerFunctional3D.h>

#include "derivedtransporttraits.h"


namespace IPP
{

template<typename Scalar,
         template<typename U> class HydrodynamicDescriptorT,
         template<typename U> class DiffusionDescriptorT,
         template<typename U, template<typename V> class DescriptorT> class HydrDynT,
         template<typename U, template<typename V> class DescriptorT> class DiffDynT,
         typename HelperFunctions
         >
struct TransportTraits3D : public DerivedTransportTraits<
        Scalar,
        HydrodynamicDescriptorT,
        DiffusionDescriptorT,
        HydrDynT,
        DiffDynT,
        HelperFunctions,

        plb::MultiBlockLattice3D
        >
{
    static constexpr size_t dim = 3;
};

}

#endif // TRANSPORTTRAITS3D_H
