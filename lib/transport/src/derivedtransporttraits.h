#ifndef DERIVEDTRANSPORTTRAITS_H
#define DERIVEDTRANSPORTTRAITS_H

#include <memory>

#include <palabos/io/vtkDataOutput.hh>

#include "plbtypededuction.h"


namespace IPP
{


template<typename Scalar_,
         template<typename U> class HydrodynamicDescT,
         template<typename U> class DiffusionDescT,
         template<typename U, template<typename V> class DescriptorT> class HydrDynT,
         template<typename U, template<typename V> class DescriptorT> class DiffDynT,
         typename HelperFunctions,
         template<typename U, template<typename V> class DescriptorT> class LatticeT
         >
struct DerivedTransportTraits
{
    using Scalar = Scalar_;

    // for some reason these types must be defined at the very top class...
    template<typename T>
    using HydrodynamicDescriptorT = HydrodynamicDescT<T>;

    template<typename T>
    using DiffusionDescriptorT = DiffusionDescT<T>;


    // hydrodynamic types
    using HydrodynamicDescriptor = HydrodynamicDescriptorT<Scalar>;
    using HydrodynamicLattice = LatticeT<Scalar, HydrodynamicDescriptorT>;
    using HydrodynamicLatticePtr = std::shared_ptr<HydrodynamicLattice>;
    // hydrodynamic dynamic types
    using HydrodynamicDynamic = HydrDynT<Scalar, HydrodynamicDescriptorT>;


    // diffusive types
    using DiffusionDescriptor = DiffusionDescriptorT<Scalar>;
    using DiffusionLattice = LatticeT<Scalar, DiffusionDescriptorT>;
    using DiffusionLatticePtr = std::shared_ptr<DiffusionLattice>;
    // diffusive dynamic types
    using DiffusionDynamic = DiffDynT<Scalar, DiffusionDescriptorT>;


    static constexpr size_t dim = DiffusionDescriptor::d;



    using ScalarFieldPtr = decltype(plb::computeDensity(std::declval<DiffusionLattice&>()));
    using ScalarField = typename ScalarFieldPtr::element_type;
    // TODO: check if this is still needed
    using ScalarFieldSharedPtr = std::shared_ptr<ScalarField>;



    // work around problems with intel compiler
    // using TensorFieldPtr = decltype(plb::computeVelocity(std::declval<DiffusionLattice&>()));
    // using TensorField = typename TensorFieldPtr::element_type;
    using TensorField = typename PlbTypeDeduction::GetMultiTensorField_XD<Scalar, dim, dim>::value;

private:
    template<typename T>
    using PlbPtrType = std::auto_ptr<T>;

public:
    using TensorFieldPtr = PlbPtrType<TensorField>;
    // TODO: check if this is still needed
    using TensorFieldSharedPtr = std::shared_ptr<TensorField>;



    template<typename T>
    using DiffusionOmegaCalc = typename HelperFunctions::template DiffusionOmegaCalc<T>;

    using PorosityFunc = typename HelperFunctions::PorosityFunc;

    static constexpr bool isTRT = HelperFunctions::isTRT;

private:
    using ScalarBlockMap = typename ScalarField::BlockMap;
    using ScalarMappedPtr = typename ScalarBlockMap::mapped_type;
    using ScalarMapped = typename std::remove_pointer<ScalarMappedPtr>::type;

    using TensorBlockMap = typename TensorField::BlockMap;
    using TensorMappedPtr = typename TensorBlockMap::mapped_type;
    using TensorMapped = typename std::remove_pointer<TensorMappedPtr>::type;


public:
    using Vector = plb::Array<Scalar,dim>;
    using Dot = decltype(std::declval<ScalarMapped>().getLocation());
    using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;
    using Box = decltype(std::declval<ScalarMapped>().getBoundingBox());


    using VtkOut = typename PlbTypeDeduction::GetVtkOutXD<Scalar, dim>::value;


};

}

#endif // DERIVEDTRANSPORTTRAITS_H
