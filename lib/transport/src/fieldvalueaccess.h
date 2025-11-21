#ifndef FIELDVALUEACCESS_H
#define FIELDVALUEACCESS_H

#include <palabos/atomicBlock/dataField2D.h>
#include <palabos/atomicBlock/dataField3D.h>

#include "array_view.h"

namespace IPP
{

template<size_t dim>
struct FieldValueAccess;

template<>
struct FieldValueAccess<2>
{
    template<typename T, template<typename U> class Field, typename I>
    static T& get(Field<T>& field, const I x, const I y)
    {
        return field.get(x, y);
    }

    template<typename T, template<typename U> class Field, typename I>
    static const T& get(const Field<T>& field, const I x, const I y)
    {
        return field.get(x, y);
    }


    template<typename T, typename I>
    static T& get(av::array_view<T, 2>& field, const I x, const I y)
    {
        return field[x][y];
    }

    template<typename T, typename I>
    static const T& get(const av::array_view<T, 2>& field, const I x, const I y)
    {
        return field[x][y];
    }



    template<typename T, template<typename U> class Field, typename I>
    static T& get(Field<T>& field, const I x, const I y, const I /*z*/)
    {
        return field.get(x, y);
    }

    template<typename T, template<typename U> class Field, typename I>
    static const T& get(const Field<T>& field, const I x, const I y, const I /*z*/)
    {
        return field.get(x, y);
    }

    template<typename T, typename I>
    static T& get(av::array_view<T, 2>& field, const I x, const I y, const I /*z*/)
    {
        return field[x][y];
    }

    template<typename T, typename I>
    static const T& get(const av::array_view<T, 2>& field, const I x, const I y, const I /*z*/)
    {
        return field[x][y];
    }


    template<typename T, typename Array>
    static T& get(av::array_view<T, 2>& field, const Array& arr)
    {
        return get(field, arr[0], arr[1]);
    }

    template<typename T, typename Array>
    static const T& get(const av::array_view<T, 2>& field, const Array& arr)
    {
        return get(field, arr[0], arr[1]);
    }


};

template<>
struct FieldValueAccess<3>
{
    template<typename T, template<typename U> class Field, typename I>
    static T& get(Field<T>& field, const I x, const I y, const I z)
    {
        return field.get(x, y, z);
    }

    template<typename T, template<typename U> class Field, typename I>
    static const T& get(const Field<T>& field, const I x, const I y, const I z)
    {
        return field.get(x, y, z);
    }


    template<typename T, typename I>
    static T& get(av::array_view<T, 3>& field, const I x, const I y, const I z)
    {
        return field[x][y][z];
    }

    template<typename T, typename I>
    static const T& get(const av::array_view<T, 3>& field, const I x, const I y, const I z)
    {
        return field[x][y][z];
    }


    template<typename T, typename Array>
    static T& get(av::array_view<T, 3>& field, const Array& arr)
    {
        return get(field, arr[0], arr[1], arr[2]);
    }

    template<typename T, typename Array>
    static const T& get(const av::array_view<T, 3>& field, const Array& arr)
    {
        return get(field, arr[0], arr[1], arr[2]);
    }

};

}

#endif // FIELDVALUEACCESS_H
