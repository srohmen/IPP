#ifndef PALABOSTYPEDEDUCTIONHELPER_H
#define PALABOSTYPEDEDUCTIONHELPER_H

#include <iostream>
#include <type_traits>

#include "palabosfwd.h"

#include <palabos/multiBlock/multiBlock2D.h>
#include <palabos/multiBlock/multiBlock3D.h>

#include <palabos/boundaryCondition/boundaryCondition2D.h>
#include <palabos/boundaryCondition/boundaryCondition3D.h>

#include <palabos/complexDynamics/advectionDiffusionBoundaryCondition2D.h>
#include <palabos/complexDynamics/advectionDiffusionBoundaryCondition3D.h>

#include <palabos/atomicBlock/dataProcessingFunctional2D.h>
#include <palabos/atomicBlock/dataProcessingFunctional3D.h>


namespace IPP
{

namespace PlbTypeDeduction
{

template <typename Descriptor>
class HasTransportScalar
{
private:
    typedef char one;
    typedef long two;

    template <typename C> static one test( decltype(&C::ExternalField::transportScalarBeginsAt) ) ;
    template <typename C> static two test(...);

public:
    enum { value = sizeof(test<Descriptor>(0)) == sizeof(char) };

};


template <typename Descriptor>
class HasPorosity
{
private:
    typedef char one;
    typedef long two;

    template <typename C> static one test( decltype(&C::ExternalField::porosityBeginsAt) ) ;
    template <typename C> static two test(...);

public:
    enum { value = sizeof(test<Descriptor>(0)) == sizeof(char) };

};



template<typename T>
typename std::enable_if<HasTransportScalar<T>::value, void>::type checkTransportScalar()
{
    std::cout << "transport scalar present" << std::endl;
}

template<typename T>
typename std::enable_if<!HasTransportScalar<T>::value, void>::type checkTransportScalar()
{
    std::cout << "NO transport scalar" << std::endl;
}



template <typename Descriptor>
class HasScalar
{
private:
    typedef char one;
    typedef long two;

    template <typename C> static one test( decltype(&C::ExternalField::scalarBeginsAt) ) ;
    template <typename C> static two test(...);

public:
    enum { value = sizeof(test<Descriptor>(0)) == sizeof(char) };

};


template<typename>
struct GetArrayDim;

template<typename T, plb::pluint size>
struct GetArrayDim< plb::Array<T, size> >
{
    static constexpr plb::pluint value = size;
};


template<size_t dim>
struct GetDotXD;

template<>
struct GetDotXD<2>
{
    using value = plb::Dot2D;
};

template<>
struct GetDotXD<3>
{
    using value = plb::Dot3D;
};



template<size_t dim>
struct GetDotListXD;

template<>
struct GetDotListXD<2>
{
    using value = plb::DotList2D;
};

template<>
struct GetDotListXD<3>
{
    using value = plb::DotList3D;
};


template<size_t dim>
struct GetBoxXD;

template<>
struct GetBoxXD<2>
{
    using value = plb::Box2D;
};

template<>
struct GetBoxXD<3>
{
    using value = plb::Box3D;
};


template<size_t dim>
struct GetScalarFieldT_XD;

template<>
struct GetScalarFieldT_XD<2>
{
    template<typename T>
    using value = plb::ScalarField2D<T>;
};

template<>
struct GetScalarFieldT_XD<3>
{
    template<typename T>
    using value = plb::ScalarField3D<T>;
};


template<typename T, size_t dim>
struct GetScalarField_XD;

template<typename T>
struct GetScalarField_XD<T, 2>
{
    using value = plb::ScalarField2D<T>;
};

template<typename T>
struct GetScalarField_XD<T, 3>
{
    using value = plb::ScalarField3D<T>;
};




template<typename T, size_t dim>
struct GetMultiScalarField_XD;

template<typename T>
struct GetMultiScalarField_XD<T, 2>
{
    using value = plb::MultiScalarField2D<T>;
};

template<typename T>
struct GetMultiScalarField_XD<T, 3>
{
    using value = plb::MultiScalarField3D<T>;
};


template<typename T, size_t tensorDim, size_t fieldDim>
struct GetMultiTensorField_XD;

template<typename T, size_t tensorDim>
struct GetMultiTensorField_XD<T, tensorDim, 2>
{
    using value = plb::MultiTensorField2D<T, tensorDim>;
};

template<typename T, size_t tensorDim>
struct GetMultiTensorField_XD<T, tensorDim, 3>
{
    using value = plb::MultiTensorField3D<T, tensorDim>;
};



template<typename T, size_t dim>
struct GetNTensorField_XD;

template<typename T>
struct GetNTensorField_XD<T, 2>
{
    using value = plb::NTensorField2D<T>;
};

template<typename T>
struct GetNTensorField_XD<T, 3>
{
    using value = plb::NTensorField3D<T>;
};


template<typename T, size_t dim>
struct GetMultiNTensorField_XD;

template<typename T>
struct GetMultiNTensorField_XD<T, 2>
{
    using value = plb::MultiNTensorField2D<T>;
};

template<typename T>
struct GetMultiNTensorField_XD<T, 3>
{
    using value = plb::MultiNTensorField3D<T>;
};


template<typename Scalar, template<typename U> class Descriptor>
struct GetBlockLattice_XD
{
private:
    template<typename T, template<typename U> class Desc, size_t dim>
    struct Helper;

    template<typename T, template<typename U> class Desc>
    struct Helper<T, Desc, 2>
    {
        using value = plb::BlockLattice2D<T, Descriptor>;
    };

    template<typename T, template<typename U> class Desc>
    struct Helper<T, Desc, 3>
    {
        using value = plb::BlockLattice3D<T, Descriptor>;
    };

public:
    using value = typename Helper<Scalar, Descriptor,
                                  Descriptor<Scalar>::d>::value;
};





template<size_t dim>
struct GetMultiBlock_XD;

template<>
struct GetMultiBlock_XD<2>
{
    using value = plb::MultiBlock2D;
};

template<>
struct GetMultiBlock_XD<3>
{
    using value = plb::MultiBlock3D;
};


template<size_t dim>
struct GetAtomicBlock_XD;

template<>
struct GetAtomicBlock_XD<2>
{
    using value = plb::AtomicBlock2D;
};

template<>
struct GetAtomicBlock_XD<3>
{
    using value = plb::AtomicBlock3D;
};



template<typename T, size_t dim>
struct GetVtkOutXD;

template<typename T>
struct GetVtkOutXD<T, 2>
{
    using value = plb::VtkImageOutput2D<T>;
};

template<typename T>
struct GetVtkOutXD<T, 3>
{
    using value = plb::VtkImageOutput3D<T>;
};


template<typename Scalar, template<typename U> class DescriptorT>
class GetOneCellIndexedFunctionalXD
{
private:
    template<typename T, template<typename U> class Descriptor, size_t dim>
    struct Helper;

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 2>
    {
        using value = plb::OneCellIndexedFunctional2D<T,Descriptor>;
    };

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 3>
    {
        using value = plb::OneCellIndexedFunctional3D<T,Descriptor>;
    };

public:
    using value = typename Helper<Scalar, DescriptorT, DescriptorT<Scalar>::d>::value;

};


template<typename Scalar, template<typename U> class DescriptorT>
class GetBoxProcessingFunctionalXD_L
{
private:

    template<typename T, template<typename U> class Descriptor, size_t dim>
    struct Helper;

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 2>
    {
        using value = plb::BoxProcessingFunctional2D_L<T,Descriptor>;
    };

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 3>
    {
        using value = plb::BoxProcessingFunctional3D_L<T,Descriptor>;
    };

public:
    using value = typename Helper<Scalar, DescriptorT, DescriptorT<Scalar>::d>::value;

};



template<typename Scalar, template<typename U> class DescriptorT>
class GetOnLatticeBoundaryConditionXD
{
private:

    template<typename T, template<typename U> class Descriptor, size_t dim>
    struct Helper;

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 2>
    {
        using value = plb::OnLatticeBoundaryCondition2D<T,Descriptor>;
    };

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 3>
    {
        using value = plb::OnLatticeBoundaryCondition3D<T,Descriptor>;
    };

public:
    using value = typename Helper<Scalar, DescriptorT, DescriptorT<Scalar>::d>::value;

};



template<typename Scalar, template<typename U> class DescriptorT>
class GetOnLatticeAdvectionDiffusionBoundaryConditionXD
{
private:

    template<typename T, template<typename U> class Descriptor, size_t dim>
    struct Helper;

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 2>
    {
        using value = plb::OnLatticeAdvectionDiffusionBoundaryCondition2D<T,Descriptor>;
    };

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 3>
    {
        using value = plb::OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>;
    };

public:
    using value = typename Helper<Scalar, DescriptorT, DescriptorT<Scalar>::d>::value;

};



template<typename T, size_t dim>
struct GetBoxProcessingFunctionalXD_S;

template<typename T>
struct GetBoxProcessingFunctionalXD_S<T, 2>
{
    using value = plb::BoxProcessingFunctional2D_S<T>;
};

template<typename T>
struct GetBoxProcessingFunctionalXD_S<T, 3>
{
    using value = plb::BoxProcessingFunctional3D_S<T>;
};

template<typename T, size_t tensorDim, size_t dim>
struct GetBoxProcessingFunctionalXD_T;

template<typename T, size_t tensorDim>
struct GetBoxProcessingFunctionalXD_T<T, tensorDim, 2>
{
    using value = plb::BoxProcessingFunctional2D_T<T, tensorDim>;
};

template<typename T, size_t tensorDim>
struct GetBoxProcessingFunctionalXD_T<T, tensorDim, 3>
{
    using value = plb::BoxProcessingFunctional3D_T<T, tensorDim>;
};


template<typename T1, typename T2, size_t dim>
struct GetBoxProcessingFunctionalXD_SS;

template<typename T1, typename T2>
struct GetBoxProcessingFunctionalXD_SS<T1, T2, 2>
{
    using value = plb::BoxProcessingFunctional2D_SS<T1,T2>;
};

template<typename T1, typename T2>
struct GetBoxProcessingFunctionalXD_SS<T1, T2, 3>
{
    using value = plb::BoxProcessingFunctional3D_SS<T1,T2>;
};



template<typename T1, typename T2, size_t tensorDim, size_t dim>
struct GetBoxProcessingFunctionalXD_ST;

template<typename T1, typename T2, size_t tensorDim>
struct GetBoxProcessingFunctionalXD_ST<T1, T2, tensorDim, 2>
{
    using value = plb::BoxProcessingFunctional2D_ST<T1,T2,tensorDim>;
};

template<typename T1, typename T2, size_t tensorDim>
struct GetBoxProcessingFunctionalXD_ST<T1, T2, tensorDim, 3>
{
    using value = plb::BoxProcessingFunctional3D_ST<T1,T2,tensorDim>;
};



template<typename T, size_t dim>
struct GetBoxProcessingFunctionalXD_N;

template<typename T>
struct GetBoxProcessingFunctionalXD_N<T, 2>
{
    using value = plb::BoxProcessingFunctional2D_N<T>;
};

template<typename T>
struct GetBoxProcessingFunctionalXD_N<T, 3>
{
    using value = plb::BoxProcessingFunctional3D_N<T>;
};


template<typename T1, typename T2, size_t dim>
struct GetBoxProcessingFunctionalXD_SN;

template<typename T1, typename T2>
struct GetBoxProcessingFunctionalXD_SN<T1, T2, 2>
{
    using value = plb::BoxProcessingFunctional2D_SN<T1, T2>;
};

template<typename T1, typename T2>
struct GetBoxProcessingFunctionalXD_SN<T1, T2, 3>
{
    using value = plb::BoxProcessingFunctional3D_SN<T1, T2>;
};



template<typename T1, typename T2, size_t dim>
struct GetReductiveBoxProcessingFunctionalXD_SS;

template<typename T1, typename T2>
struct GetReductiveBoxProcessingFunctionalXD_SS<T1, T2, 2>
{
    using value = plb::ReductiveBoxProcessingFunctional2D_SS<T1,T2>;
};

template<typename T1, typename T2>
struct GetReductiveBoxProcessingFunctionalXD_SS<T1, T2, 3>
{
    using value = plb::ReductiveBoxProcessingFunctional3D_SS<T1,T2>;
};





template<size_t dim>
struct GetDotProcessingFunctionalXD;

template<>
struct GetDotProcessingFunctionalXD<2>
{
    using value = plb::DotProcessingFunctional2D;
};

template<>
struct GetDotProcessingFunctionalXD<3>
{
    using value = plb::DotProcessingFunctional3D;
};




template<typename T, size_t dim>
struct GetDotProcessingFunctionalXD_S;

template<typename T>
struct GetDotProcessingFunctionalXD_S<T,2>
{
    using value = plb::DotProcessingFunctional2D_S<T>;
};

template<typename T>
struct GetDotProcessingFunctionalXD_S<T,3>
{
    using value = plb::DotProcessingFunctional3D_S<T>;
};



template<typename T, size_t dim>
struct GetReductiveDotProcessingFunctionalXD_S;

template<typename T>
struct GetReductiveDotProcessingFunctionalXD_S<T,2>
{
    using value = plb::ReductiveDotProcessingFunctional2D_S<T>;
};

template<typename T>
struct GetReductiveDotProcessingFunctionalXD_S<T,3>
{
    using value = plb::ReductiveDotProcessingFunctional3D_S<T>;
};




template<typename T1, typename T2, size_t dim>
struct GetDotProcessingFunctionalXD_SS;

template<typename T1, typename T2>
struct GetDotProcessingFunctionalXD_SS<T1, T2, 2>
{
    using value = plb::DotProcessingFunctional2D_SS<T1,T2>;
};

template<typename T1, typename T2>
struct GetDotProcessingFunctionalXD_SS<T1, T2, 3>
{
    using value = plb::DotProcessingFunctional3D_SS<T1,T2>;
};


template<typename T, size_t dim>
struct GetDotProcessingFunctionalXD_N;

template<typename T>
struct GetDotProcessingFunctionalXD_N<T, 2>
{
    using value = plb::DotProcessingFunctional2D_N<T>;
};

template<typename T>
struct GetDotProcessingFunctionalXD_N<T, 3>
{
    using value = plb::DotProcessingFunctional3D_N<T>;
};



template<size_t dim>
struct GetBoxProcessingFunctionalXD;

template<>
struct GetBoxProcessingFunctionalXD<2>
{
    using value = plb::BoxProcessingFunctional2D;
};

template<>
struct GetBoxProcessingFunctionalXD<3>
{
    using value = plb::BoxProcessingFunctional3D;
};



template<typename T, size_t dim>
struct GetScalarFieldBoxProcessingFunctionalXD;

template<typename T>
struct GetScalarFieldBoxProcessingFunctionalXD<T, 2>
{
    using value = plb::ScalarFieldBoxProcessingFunctional2D<T>;
};

template<typename T>
struct GetScalarFieldBoxProcessingFunctionalXD<T, 3>
{
    using value = plb::ScalarFieldBoxProcessingFunctional3D<T>;
};



template<typename T, size_t dim>
struct GetScalarFieldDotProcessingFunctionalXD;

template<typename T>
struct GetScalarFieldDotProcessingFunctionalXD<T, 2>
{
    using value = plb::ScalarFieldDotProcessingFunctional2D<T>;
};

template<typename T>
struct GetScalarFieldDotProcessingFunctionalXD<T, 3>
{
    using value = plb::ScalarFieldDotProcessingFunctional3D<T>;
};



template<typename Scalar, template<typename U> class Descriptor>
struct GetLatticeBoxProcessingFunctionalXD
{
private:
    template<typename T, template<typename U> class Desc, size_t dim>
    struct Helper;

    template<typename T, template<typename U> class Desc>
    struct Helper<T, Desc, 2>
    {
        using value = plb::LatticeBoxProcessingFunctional2D<T, Desc>;
    };

    template<typename T, template<typename U> class Desc>
    struct Helper<T, Desc, 3>
    {
        using value = plb::LatticeBoxProcessingFunctional3D<T, Desc>;
    };


public:
    using value = typename Helper<Scalar, Descriptor,
                                  Descriptor<Scalar>::d>::value;
};





template<typename Scalar, template<typename U> class DescriptorT>
struct GetDotProcessingFunctionalXD_L
{
private:
    template<typename T, template<typename U> class Descriptor, size_t dim>
    struct Helper;

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 2>
    {
        using value = plb::DotProcessingFunctional2D_L<T,Descriptor>;
    };

    template<typename T, template<typename U> class Descriptor>
    struct Helper<T, Descriptor, 3>
    {
        using value = plb::DotProcessingFunctional3D_L<T,Descriptor>;
    };


public:
    using value = typename Helper<Scalar, DescriptorT,
                                  DescriptorT<Scalar>::d>::value;
};



template<typename LatticeScalar,
         template<typename U> class DescriptorT,
         typename FieldScalar = LatticeScalar>
class GetDotProcessingFunctionalXD_LS
{
private:
    template<typename T1, template<typename U> class Descriptor, typename T2, size_t dim>
    struct Helper;

    template<typename T1, template<typename U> class Descriptor, typename T2>
    struct Helper<T1, Descriptor, T2, 2>
    {
        using value = plb::DotProcessingFunctional2D_LS<T1,Descriptor,T2>;
    };

    template<typename T1, template<typename U> class Descriptor, typename T2>
    struct Helper<T1, Descriptor, T2, 3>
    {
        using value = plb::DotProcessingFunctional3D_LS<T1,Descriptor,T2>;
    };


public:
    using value = typename Helper<LatticeScalar, DescriptorT,
                                  FieldScalar, DescriptorT<LatticeScalar>::d>::value;
};






template<typename LatticeScalar, template<typename U> class DescriptorT, typename FieldScalar = LatticeScalar>
class GetBoxProcessingFunctionalXD_LS
{
private:

    template<typename T1, template<typename U> class Descriptor, typename T2, size_t dim>
    struct Helper;

    template<typename T1, template<typename U> class Descriptor, typename T2>
    struct Helper<T1, Descriptor, T2, 2>
    {
        using value = plb::BoxProcessingFunctional2D_LS<T1, Descriptor, T2>;
    };

    template<typename T1, template<typename U> class Descriptor, typename T2>
    struct Helper<T1, Descriptor, T2, 3>
    {
        using value = plb::BoxProcessingFunctional3D_LS<T1, Descriptor, T2>;
    };

public:
    using value = typename Helper<LatticeScalar, DescriptorT,
                                  FieldScalar, DescriptorT<LatticeScalar>::d>::value;

};


template<class F>
struct function_traits;

template<class R, class... Args>
struct function_traits<R(Args...)>
{
    using return_type = R;

    static constexpr std::size_t arity = sizeof...(Args);

    template <std::size_t N>
    struct argument
    {
        static_assert(N < arity, "error: invalid parameter index.");
        using type = typename std::tuple_element<N,std::tuple<Args...>>::type;
    };
};

// function pointer
template<class R, class... Args>
struct function_traits<R(*)(Args...)> : public function_traits<R(Args...)>
{};

// member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...)> : public function_traits<R(C&,Args...)>
{};

// const member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...) const> : public function_traits<R(C&,Args...)>
{};




template<typename DotList>
struct GetDotDim;

template<>
struct GetDotDim<plb::Dot2D>
{
    static constexpr size_t value = 2;
};

template<>
struct GetDotDim<plb::Dot3D>
{
    static constexpr size_t value = 3;
};


template<typename DotList>
struct GetDotListDim;

template<>
struct GetDotListDim<plb::DotList2D>
{
    static constexpr size_t value = 2;
};

template<>
struct GetDotListDim<plb::DotList3D>
{
    static constexpr size_t value = 3;
};


template<typename Box>
struct GetBoxDim;

template<>
struct GetBoxDim<plb::Box2D>
{
    static constexpr size_t value = 2;
};

template<>
struct GetBoxDim<plb::Box3D>
{
    static constexpr size_t value = 3;
};



template<typename Scalar, typename ScalarField, typename Enable = void>
struct GetFieldDim;

template<typename Scalar, typename ScalarField>
struct GetFieldDim<Scalar, ScalarField, typename std::enable_if<
        std::is_base_of<plb::ScalarFieldBase2D<Scalar>, ScalarField>::value
        >::type>
{
    static constexpr size_t value = 2;
};

template<typename Scalar, typename ScalarField>
struct GetFieldDim<Scalar, ScalarField, typename std::enable_if<
        std::is_base_of<plb::ScalarFieldBase3D<Scalar>, ScalarField>::value
        >::type>
{
    static constexpr size_t value = 3;
};


template<typename T, typename>
struct GetTensorDim;

template<typename T>
struct GetTensorDim<T, plb::MultiScalarField2D<T> >
{
    static const std::size_t value = 1;
};

template<typename T>
struct GetTensorDim<T, plb::MultiScalarField3D<T> >
{
    static const std::size_t value = 1;
};

template<typename T, std::size_t N>
struct GetTensorDim<T, plb::MultiTensorField2D<T, N> >
{
    static const std::size_t value = N;
};

template<typename T, std::size_t N>
struct GetTensorDim<T, plb::MultiTensorField3D<T, N> >
{
    static const std::size_t value = N;
};



template<typename T, typename Field, typename Enable = void>
struct GetBoxProcessingFunc;

template<typename T, typename Field>
struct GetBoxProcessingFunc
        <
            T, Field, typename std::enable_if
            <
                GetTensorDim<T, Field>::value == 1
            >::type
        >
{
    static constexpr size_t dim = PlbTypeDeduction::GetFieldDim<T, Field>::value;
    using value = typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T, dim>::value;
};


template<typename T, typename Field>
struct GetBoxProcessingFunc
        <
            T, Field, typename std::enable_if
            <
                GetTensorDim<T, Field>::value != 1
            >::type
        >
{
    static constexpr size_t tensorDim = GetTensorDim<T, Field>::value;
    static constexpr size_t dim = PlbTypeDeduction::GetFieldDim<T, Field>::value;
    using value = typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_T<T, tensorDim, dim>::value;
};




template<size_t dim>
struct GetSparseBlockStructure_XD;

template<>
struct GetSparseBlockStructure_XD<2>
{
    using value = plb::SparseBlockStructure2D;
};

template<>
struct GetSparseBlockStructure_XD<3>
{
    using value = plb::SparseBlockStructure3D;
};


} // end of namespace PalabosTypeDeductionHelper


}


#endif // PALABOSTYPEDEDUCTIONHELPER_H
