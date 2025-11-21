#ifndef INPUTGEOMETRYSCALE_H
#define INPUTGEOMETRYSCALE_H

#include <cstddef>
#include <vector>

namespace IPP
{

namespace InputGeometryScale
{


class ColumnMajorArrayDimensionConvert
{
public:
    ColumnMajorArrayDimensionConvert(const size_t nx, const size_t ny, const size_t nz = 1)
        : m_nx(nx),
          m_ny(ny),
          m_nz(nz),
          m_nxy(nx*ny)

    {

    }

    inline size_t getNxyz() const
    {
        return m_nxy * m_nz;
    }

    inline size_t calcIndex(const size_t x, const size_t y, const size_t z) const
    {
        return z * m_nxy + y * m_nx + x;
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

private:
    size_t m_nx;
    size_t m_ny;
    size_t m_nz;
    size_t m_nxy;
};


size_t calcUpscaledSize(const size_t n, const size_t scalePasses);
size_t calcDownscaledSize(const size_t n, const size_t scalePasses);

void scale(const ColumnMajorArrayDimensionConvert& srcIndexConv,
           const ColumnMajorArrayDimensionConvert& dstIndexConv,
           const std::vector<size_t>& input,
           std::vector<size_t>& output);
}

}

#endif // INPUTGEOMETRYSCALE_H
