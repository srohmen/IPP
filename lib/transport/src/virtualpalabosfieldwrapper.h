#ifndef VIRTUALPALABOSSCALARFIELDWRAPPER_H
#define VIRTUALPALABOSSCALARFIELDWRAPPER_H

#include "abstractscalarfield.h"
#include "palabosfieldwrapper.h"

namespace IPP
{

template<typename T, typename ScalarField, size_t dim>
class VirtualPalabosFieldWrapper : public AbstractScalarField
{
public:
    VirtualPalabosFieldWrapper(ScalarField& field)
        : m_field(field)
    {

    }

    virtual ~VirtualPalabosFieldWrapper()
    {

    }

    virtual size_t getNx() const
    {
        return m_field.getNx();
    }

    virtual size_t getNy() const
    {
        return m_field.getNy();
    }

    virtual size_t getNz() const
    {
        return m_field.getNz();
    }

    virtual double get(int x, int y, int z) const
    {
        const T& entry = m_field.get(x, y, z);
        return static_cast<double>(entry);
    }

    virtual void set(int x, int y, int z, const double& val)
    {
        T& entry = m_field.get(x, y, z);
        entry = val;
    }

private:
    PalabosFieldWrapper<ScalarField, dim> m_field;
};


}

#endif // VIRTUALPALABOSSCALARFIELDWRAPPER_H
