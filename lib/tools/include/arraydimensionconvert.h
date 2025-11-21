#ifndef ARRAYDIMENSIONCONVERT_H
#define ARRAYDIMENSIONCONVERT_H

#include <cassert>
#include "ippvector.h"

namespace IPP
{

class ArrayDimensionConvert
{
public:
    ArrayDimensionConvert()
        : m_nx(0),
          m_ny(0),
          m_nz(0),
          m_nxy(0),
          m_nyz(0)
    {

    }

    ArrayDimensionConvert(const size_t nx, const size_t ny, const size_t nz = 1)
        : m_nx(nx),
          m_ny(ny),
          m_nz(nz),
          m_nxy(nx * ny),
          m_nyz(ny * nz)
    {

    }

    template<typename Vector>
    ArrayDimensionConvert(const Vector& size)
        : m_nx(size[0]),
          m_ny(size[1]),
          m_nz(size[2]),
          m_nxy(size[0] * size[1]),
          m_nyz(size[1] * size[2])
    {

    }

    inline bool isInBounds(const int x, const int y, const int z) const
    {
        return isInRange(x, 0, (int)m_nx) &&
                isInRange(y, 0, (int)m_ny) &&
                isInRange(z, 0, (int)m_nz);
    }

    inline size_t getNxyz() const
    {
        return m_nx * m_nyz;
    }

    inline size_t calcIndex(const int x, const int y, const int z) const
    {
        assert(isInBounds(x,y,z));
        return x * m_nyz + y * m_nz + z;
    }    

    template<typename CoordType>
    inline size_t calcIndex(const CoordType& coords) const
    {
        return calcIndex(coords[0], coords[1], coords[2]);
    }

    inline IPPVector3DInt calcCoordinate(const size_t index) const
    {
        IPPVector3DInt coord;

        const std::div_t xDiv = std::div((int)index, (int)m_nyz);
        coord[0] = xDiv.quot;

        const std::div_t yDiv = std::div(xDiv.rem, (int)m_nz);
        coord[1] = yDiv.quot;

        coord[2] = yDiv.rem;

        return coord;
    }

    size_t getNx() const
    {
        return m_nx;
    }

    size_t getNy() const
    {
        return m_ny;
    }

    size_t getNz() const
    {
        return m_nz;
    }

    void getSize(size_t& nx, size_t& ny, size_t& nz) const
    {
        nx = m_nx;
        ny = m_ny;
        nz = m_nz;
    }

    inline IPPVector3DInt colMajorToRowMajorCoord(const size_t index) const
    {
        const std::div_t zDiv = std::div((int)index, (int)m_nxy);
        const size_t z = zDiv.quot;

        const std::div_t yDiv = std::div(zDiv.rem, (int)m_nx);
        const size_t y = yDiv.quot;

        const size_t x = yDiv.rem;

        return { {(int)x, (int)y, (int)z} };
    }

private:

    template <typename T>
    bool isInRange(const T& value, const T& low, const T& high) const
    {
        return !(value < low) && (value < high);
    }


    size_t m_nx;
    size_t m_ny;
    size_t m_nz;

    size_t m_nxy;
    size_t m_nyz;
};

}

#endif // ARRAYDIMENSIONCONVERT_H
