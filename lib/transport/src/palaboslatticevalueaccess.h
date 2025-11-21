#ifndef PALABOSLATTICEVALUEACCESS_H
#define PALABOSLATTICEVALUEACCESS_H

#include <palabos/dataProcessors/dataAnalysisWrapper2D.h>
#include <palabos/dataProcessors/dataInitializerWrapper2D.h>

#include <palabos/dataProcessors/dataAnalysisWrapper3D.hh>
#include <palabos/dataProcessors/dataInitializerWrapper3D.hh>
#include <palabos/dataProcessors/dataInitializerFunctional3D.hh>

#include "plbtypededuction.h"
#include "array_view.h"
#include "ippvector.h"
#include "domainiterator.h"
#include "palabosdataaccess.h"
#include "arrayviewdataaccess.h"

namespace IPP
{

template<typename T, size_t dim>
class SetScalarsFromArray
{
public:
    SetScalarsFromArray(const av::array_view<T, dim>& scalarArray,
                        const IPPVector3DLong& location)
        : m_scalarArray(scalarArray),
          location(location)
    {

    }

    T operator()(plb::plint x, plb::plint y) const
    {
        // std::cerr << x << "\t" << y << std::endl;
        const T& scalar = m_scalarArray[x - location[0]][y - location[1]];
        assert(std::isnan(scalar) == false);
        return scalar;
    }

    T operator()(plb::plint x, plb::plint y, plb::plint z) const
    {
        const T& scalar = m_scalarArray[x - location[0]][y - location[1]][z - location[2]];
        assert(std::isnan(scalar) == false);
        return scalar;
    }

private:
    const av::array_view<T, dim>& m_scalarArray;
    const IPPVector3DLong& location;
};


template<typename T, typename Array, size_t dim>
class GetScalarsToArray
        : public PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T, dim>::value
{
public:
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;


    GetScalarsToArray(Array& result, const Dot& location)
        : m_arr(result),
          location(location)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }

    virtual GetScalarsToArray<T,Array,dim>* clone() const override
    {
        return new GetScalarsToArray<T,Array,dim>(*this);
    }


    using BaseClass = typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T, dim>::value;
    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

    virtual void process(Box domain, InputField& decomposedMultiData) override
    {

        const auto offset = decomposedMultiData.getLocation();

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto posLocal = *it;
            const T val = DataAccess::get(decomposedMultiData, posLocal);

            assert(std::isnan(val) == false);

            const auto posOutput = posLocal + offset - location;
            using ArrType = typename Array::value_type;
            DataAccess::get(m_arr, posOutput) = static_cast<ArrType>(val);

        }
    }

private:
    Array& m_arr;
    const Dot& location;
};

template<typename T, typename Array, size_t tensorDim, size_t dim>
class GetTensorsToArray
        : public PlbTypeDeduction::GetBoxProcessingFunctionalXD_T<T, tensorDim, dim>::value
{
public:
};

namespace PalabosLatticeValueAccess
{

    // set transport scalar to constant
    template <typename Descriptor, typename LatticeType, typename BoxType, typename T>
    static typename std::enable_if<PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setTransportScalarIfPossible(LatticeType& lattice, const BoxType& domain, const T& scalar)
    {
        plb::setExternalScalar(lattice, domain, Descriptor::ExternalField::transportScalarBeginsAt, scalar);
    }

    template <typename Descriptor, typename LatticeType, typename BoxType, typename T>
    static typename std::enable_if<!PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setTransportScalarIfPossible(LatticeType& lattice, const BoxType& domain, const T& scalar )
    {
        // do nothing
    }


    // set transport scalar to constant
    template <typename Descriptor, typename LatticeType, typename T>
    static typename std::enable_if<PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setTransportScalarIfPossible(LatticeType& lattice, const T& scalar)
    {
        plb::setExternalScalar(lattice, lattice.getBoundingBox(),
                               Descriptor::ExternalField::transportScalarBeginsAt, scalar);
    }

    template <typename Descriptor, typename LatticeType, typename T>
    static typename std::enable_if<!PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setTransportScalarIfPossible(LatticeType&, const T& )
    {
        // do nothing
    }


    // set transport scalar to field
    template <typename Descriptor, typename LatticeType, typename ScalarFieldType>
    static typename std::enable_if<PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setTransportScalarIfPossible(LatticeType& lattice, ScalarFieldType& scalarField)
    {
        plb::setExternalScalar(lattice, lattice.getBoundingBox(),
                               Descriptor::ExternalField::transportScalarBeginsAt,
                               scalarField);
    }

    template <typename Descriptor, typename LatticeType, typename ScalarFieldType>
    static typename std::enable_if<!PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setTransportScalarIfPossible(LatticeType&, ScalarFieldType&)
    {
        // do nothing
    }


    // set porosity to constant value
    template <typename Descriptor, typename LatticeType, typename T>
    static typename std::enable_if<PlbTypeDeduction::HasPorosity<Descriptor>::value, void>::type
    setPorosityIfPossible(LatticeType& lattice, const T& scalar)
    {
        plb::setExternalScalar(lattice, lattice.getBoundingBox(),
                               Descriptor::ExternalField::porosityBeginsAt, scalar);
    }

    template <typename Descriptor, typename LatticeType, typename T>
    static typename std::enable_if<!PlbTypeDeduction::HasPorosity<Descriptor>::value, void>::type
    setPorosityIfPossible(LatticeType& /*lattice*/, const T& /*scalar*/)
    {
        // do nothing
    }
    /////////////////////////


    // set porosity to field
    template <typename Descriptor, typename LatticeType, typename ScalarFieldType>
    static typename std::enable_if<PlbTypeDeduction::HasPorosity<Descriptor>::value, void>::type
    setPorosityIfPossible(LatticeType& lattice, ScalarFieldType& scalarField)
    {
        plb::setExternalScalar(lattice, lattice.getBoundingBox(),
                               Descriptor::ExternalField::porosityBeginsAt,
                               scalarField);
    }

    template <typename Descriptor, typename LatticeType, typename ScalarFieldType>
    static typename std::enable_if<!PlbTypeDeduction::HasPorosity<Descriptor>::value, void>::type
    setPorosityIfPossible(LatticeType& /*lattice*/, ScalarFieldType& /*scalarField*/)
    {
        // do nothing
    }
    ////////////////////////


    // get porosity
    template <typename Descriptor, typename LatticeType, typename ScalarFieldTypePtr>
    static typename std::enable_if<PlbTypeDeduction::HasPorosity<Descriptor>::value, void>::type
    computePorosityIfPossible(LatticeType& lattice, ScalarFieldTypePtr& porosityField)
    {
        porosityField = plb::computeExternalScalar(lattice, Descriptor::ExternalField::porosityBeginsAt);
    }

    template <typename Descriptor, typename LatticeType, typename ScalarFieldTypePtr>
    static typename std::enable_if<!PlbTypeDeduction::HasPorosity<Descriptor>::value, void>::type
    computePorosityIfPossible(LatticeType & /*lattice*/, ScalarFieldTypePtr& /*scalarField*/)
    {
        // do nothing
    }
    /////////////////////////////////


    // get transport scalar field
    template <typename Descriptor, typename LatticeType, typename ScalarFieldTypePtr>
    static typename std::enable_if<PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    computeTransportScalarIfPossible(LatticeType &lattice, ScalarFieldTypePtr& scalarField)
    {
        scalarField = plb::computeExternalScalar(lattice, Descriptor::ExternalField::transportScalarBeginsAt);
    }

    template <typename Descriptor, typename LatticeType, typename ScalarFieldTypePtr>
    static typename std::enable_if<!PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    computeTransportScalarIfPossible(LatticeType &lattice, ScalarFieldTypePtr& scalarField)
    {
        // do nothing
    }





    // set constant scalar
    template <typename Descriptor, typename LatticeType, typename BoxType, typename T>
    static typename std::enable_if<PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setScalarIfPossible(LatticeType& lattice, const BoxType& domain, const T& scalar )
    {
        plb::setExternalScalar(lattice, domain, Descriptor::ExternalField::scalarBeginsAt, scalar);
    }

    template <typename Descriptor, typename LatticeType, typename BoxType, typename T>
    static typename std::enable_if<!PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setScalarIfPossible(LatticeType& /*lattice*/, const BoxType& /*domain*/, const T& /*scalar*/ )
    {
        // do nothing
    }


    // set constant scalar
    template <typename Descriptor, typename LatticeType, typename T>
    static typename std::enable_if<PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setScalarIfPossible(LatticeType& lattice, const T& scalar )
    {
        setScalarIfPossible<Descriptor>(lattice, lattice.getBoundingBox(), scalar);
    }

    template <typename Descriptor, typename LatticeType, typename T>
    static typename std::enable_if<!PlbTypeDeduction::HasTransportScalar<Descriptor>::value, void>::type
    setScalarIfPossible(LatticeType& /*lattice*/, const T& /*scalar*/ )
    {
        // do nothing
    }



    // set scalar field
    template <typename Descriptor, typename LatticeType, typename ScalarFieldType>
    static typename std::enable_if<PlbTypeDeduction::HasScalar<Descriptor>::value, void>::type
    setScalarIfPossible(LatticeType& lattice, ScalarFieldType& scalarField)
    {
        plb::setExternalScalar(lattice, lattice.getBoundingBox(),
                               Descriptor::ExternalField::scalarBeginsAt,
                               scalarField);
    }

    template <typename Descriptor, typename LatticeType, typename ScalarFieldType>
    static typename std::enable_if<!PlbTypeDeduction::HasScalar<Descriptor>::value, void>::type
    setScalarIfPossible(LatticeType& lattice, ScalarFieldType& scalarField)
    {
        // do nothing
    }



} // end of namespace PalabosLatticeValueAccess

}

#endif // PALABOSLATTICEVALUEACCESS_H
