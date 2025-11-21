#ifndef ISZEROFUNC_H
#define ISZEROFUNC_H

#include <palabos/multiBlock/multiDataField2D.h>
#include <dataProcessors/dataInitializerFunctional2D.h>
#include <palabos/multiBlock/multiDataField3D.h>
#include <dataProcessors/dataInitializerFunctional3D.h>

#include "mathtools.h"

namespace IPP
{

template<typename Scalar, size_t dim>
class IsZeroFuncT;

template<typename Scalar>
class IsZeroFuncT<Scalar, 2> : public plb::DomainFunctional2D
{
public:
    IsZeroFuncT(const plb::MultiScalarField2D<Scalar>& field, const Scalar& minVal)
        : field(field),
          minVal(minVal)
    {

    }

    virtual bool operator() (plb::plint x, plb::plint y) const
    {
        const plb::plint wrappedX = MathTools::wrapToRange(x, field.getNx());
        const plb::plint wrappedY = MathTools::wrapToRange(y, field.getNy());
        const Scalar& val = field.get(wrappedX, wrappedY);
        bool isZero = val <= minVal;
        return isZero;
    }

    virtual DomainFunctional2D* clone() const
    {
        return new IsZeroFuncT<Scalar, 2>(*this);
    }

private:
    const plb::MultiScalarField2D<Scalar>& field;
    const Scalar minVal;
};

template<typename Scalar>
class IsZeroFuncT<Scalar, 3> : public plb::DomainFunctional3D
{
public:
    IsZeroFuncT(const plb::MultiScalarField3D<Scalar>& field, const Scalar& minVal)
        : field(field),
          minVal(minVal)
    {

    }

    virtual bool operator() (plb::plint x, plb::plint y, plb::plint z) const
    {
        const plb::plint wrappedX = MathTools::wrapToRange(x, field.getNx());
        const plb::plint wrappedY = MathTools::wrapToRange(y, field.getNy());
        const plb::plint wrappedZ = MathTools::wrapToRange(z, field.getNz());
        const Scalar& val = field.get(wrappedX, wrappedY, wrappedZ);
        bool isZero = val <= minVal;
        return isZero;
    }

    virtual DomainFunctional3D* clone() const
    {
        return new IsZeroFuncT<Scalar, 3>(*this);
    }

private:
    const plb::MultiScalarField3D<Scalar>& field;
    const Scalar minVal;
};

}

#endif // ISZEROFUNC_H
