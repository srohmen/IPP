#ifndef PALABOSFIELDWRAPPER_H
#define PALABOSFIELDWRAPPER_H

#include <cstddef>
#include <type_traits>
#include "is_equal.h"

namespace IPP
{

template<typename FieldType, size_t dim>
class PalabosFieldWrapper;

template<typename FieldType>
class PalabosFieldWrapper<FieldType, 2>
{
public:
    PalabosFieldWrapper(FieldType& scalarField)
        : m_field(scalarField)
    {

    }

    const FieldType& getWrapped() const
    {
        return m_field;
    }

    size_t getNx() const
    {
        return m_field.getNx();
    }

    size_t getNy() const
    {
        return m_field.getNy();
    }

    size_t getNz() const
    {
        return 1;
    }


    decltype(auto) get(int x, int y, int)
    {
        return m_field.get(x, y);
    }

    decltype(auto) get(int x, int y, int) const
    {
        return m_field.get(x, y);
    }


private:
    FieldType& m_field;
};


template<typename FieldType>
class PalabosFieldWrapper<FieldType, 3>
{
public:
    PalabosFieldWrapper(FieldType& scalarField)
        : m_field(scalarField)
    {

    }

    const FieldType& getWrapped() const
    {
        return m_field;
    }

    size_t getNx() const
    {
        return m_field.getNx();
    }

    size_t getNy() const
    {
        return m_field.getNy();
    }

    size_t getNz() const
    {
        return m_field.getNz();
    }

    decltype(auto) get(int x, int y, int z)
    {
        return m_field.get(x, y, z);
    }

    decltype(auto) get(int x, int y, int z) const
    {
        return m_field.get(x, y, z);
    }

private:
    FieldType& m_field;
};


}


#endif // PALABOSFIELDWRAPPER_H
