#include "geometrytools.h"

#include <cassert>

namespace IPP
{

void GeometryTools::transpose(std::vector<double> &arr2D, const size_t nx, const size_t ny)
{
    assert(nx * ny == arr2D.size());

    std::vector<double> transposed(arr2D.size());

    for (size_t x = 0; x < nx; x++)
    {
        for (size_t y = 0; y < ny; y++)
        {
            const size_t iSrc = y * nx + x;
            const size_t iDst = x * ny + y;
            transposed[iDst] = arr2D[iSrc];
        }
    }

    arr2D.swap(transposed);
}



}


