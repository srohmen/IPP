#ifndef ARRAYVIEWSCALARFIELD_H
#define ARRAYVIEWSCALARFIELD_H

#include "abstractscalarfield.h"
#include "array_view.h"
#include "fieldvalueaccess.h"

namespace IPP
{


template<typename ArrType>
class ArrayViewScalarField : public AbstractScalarField
{
public:
    ArrayViewScalarField(const ArrType& arr)
        : m_arr(arr)
    {

    }

    virtual size_t getNx() const
    {
        return m_arr.bounds()[0];
    }

    virtual size_t getNy() const
    {
        return m_arr.bounds()[1];
    }


    virtual size_t getNz() const
    {
        return helper<ArrType>(m_arr);
    }

    virtual double get(int x, int y, int z) const
    {
        return FieldValueAccess<ArrType::rank>::get(m_arr, x, y, z);
    }

    virtual void set(int x, int y, int z, const double& newVal)
    {
        FieldValueAccess<ArrType::rank>::get(m_arr, x, y, z) = newVal;
    }

private:

    template<typename Arr>
    typename std::enable_if<Arr::rank == 3, size_t>::type
    helper(const Arr &arr) const
    {
         return arr.bounds()[2];
    }

    template<typename Arr>
    typename std::enable_if<Arr::rank != 3, size_t>::type
    helper(const Arr &) const
    {
        return 1;
    }

    ArrType m_arr;
};


template<typename ArrType, bool takeOwnership = false>
class ConstArrayViewScalarField : public ConstAbstractScalarField
{
public:
    ConstArrayViewScalarField(const ArrType& arr)
        : m_arr(arr)
    {

    }

    virtual size_t getNx() const
    {
        return m_arr.bounds()[0];
    }

    virtual size_t getNy() const
    {
        return m_arr.bounds()[1];
    }


    virtual size_t getNz() const
    {
        return helper<ArrType>(m_arr);
    }    

    virtual double get(int x, int y, int z) const
    {
        return FieldValueAccess<ArrType::rank>::get(m_arr, x, y, z);
    }


private:

    template<typename Arr>
    typename std::enable_if<Arr::rank == 3, size_t>::type
    helper(const Arr &arr) const
    {
         return arr.bounds()[2];
    }

    template<typename Arr>
    typename std::enable_if<Arr::rank != 3, size_t>::type
    helper(const Arr &) const
    {
        return 1;
    }

    using ArrTypeRef = ArrType&;

    using ArrTypeContainer = typename
    std::conditional<takeOwnership, ArrType, const ArrType&>::type;

    const ArrTypeContainer m_arr;
};

}

#endif // ARRAYVIEWSCALARFIELD_H
