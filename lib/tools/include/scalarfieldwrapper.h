#ifndef RAWSCALARFIELDWRAPPER_H
#define RAWSCALARFIELDWRAPPER_H

#include <assert.h>
#include "abstractscalarfield.h"
#include "rawscalarfield.h"

namespace IPP
{

template<typename T>
class RawScalarFieldWrapper2D : public AbstractScalarField
{
public:
    RawScalarFieldWrapper2D(RawScalarField2D<T>& toWrap)
        : m_toWrap(toWrap)
    {

    }

    virtual size_t getNx() const
    {
        return m_toWrap.getNx();
    }

    virtual size_t getNy() const
    {
        return m_toWrap.getNy();
    }

    virtual size_t getNz() const
    {
        // z axis not defined for 2D
        return 1;
    }

    virtual T get(int x, int y, int z) const
    {
        assert(x < (int)this->getNx());
        assert(y < (int)this->getNy());
        assert(z < (int)this->getNz());
        // z axis not defined for 2D
        assert(z == 0);
        return m_toWrap.get(x, y);
    }

    virtual void set(int x, int y, int z, const T& val)
    {
        assert(x < (int)this->getNx());
        assert(y < (int)this->getNy());
        assert(z < (int)this->getNz());
        // z axis not defined for 2D
        assert(z == 0);

        m_toWrap[x][y] = val;
    }

private:
    RawScalarField2D<T>& m_toWrap;
};


template<typename T>
class RawScalarFieldWrapper3D : public AbstractScalarField
{
public:
    RawScalarFieldWrapper3D(RawScalarField3D<T>& toWrap)
        : m_toWrap(toWrap)
    {

    }

    virtual size_t getNx() const
    {
        return m_toWrap.getNx();
    }

    virtual size_t getNy() const
    {
        return m_toWrap.getNy();
    }

    virtual size_t getNz() const
    {
        return m_toWrap.getNz();
    }

    virtual T get(int x, int y, int z) const
    {
        assert(x < (int)this->getNx());
        assert(y < (int)this->getNy());
        assert(z < (int)this->getNz());

        return m_toWrap.get(x, y, z);
    }

    virtual void set(int x, int y, int z, const T& val)
    {
        assert(x < (int)this->getNx());
        assert(y < (int)this->getNy());
        assert(z < (int)this->getNz());

        m_toWrap[x][y][z] = val;
    }

private:
    RawScalarField3D<T>& m_toWrap;
};

}

#endif // RAWSCALARFIELDWRAPPER_H
