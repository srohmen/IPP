#ifndef DIMSIZECALC_H
#define DIMSIZECALC_H

#include <vector>
#include <palabos/core/globalDefs.h>
#include "is_equal.h"

namespace IPP
{

namespace detail
{

template<size_t dim>
struct DimSizeCalcHelper;

template<>
struct DimSizeCalcHelper<2>
{
    static size_t calc(const std::vector<plb::plint>& allGlobalIndices,
                       const size_t startIndex)
    {
        const size_t x0 = allGlobalIndices[startIndex];
        const size_t x1 = allGlobalIndices[startIndex + 1];
        const size_t y0 = allGlobalIndices[startIndex + 2];
        const size_t y1 = allGlobalIndices[startIndex + 3];

        const size_t nx = x1 - x0 + 1;
        const size_t ny = y1 - y0 + 1;
        const size_t size = nx*ny;
        return size;
    }
};


template<>
struct DimSizeCalcHelper<3>
{
    static size_t calc(const std::vector<plb::plint>& allGlobalIndices,
                       const size_t startIndex)
    {
        const size_t x0 = allGlobalIndices[startIndex];
        const size_t x1 = allGlobalIndices[startIndex + 1];
        const size_t y0 = allGlobalIndices[startIndex + 2];
        const size_t y1 = allGlobalIndices[startIndex + 3];
        const size_t z0 = allGlobalIndices[startIndex + 4];
        const size_t z1 = allGlobalIndices[startIndex + 5];

        const size_t nx = x1 - x0 + 1;
        const size_t ny = y1 - y0 + 1;
        const size_t nz = z1 - z0 + 1;
        const size_t size = nx*ny*nz;
        return size;
    }
};

}

template<size_t dimension>
struct DimSizeCalc
{
    static constexpr size_t dim = dimension;
    static constexpr size_t nComp = 2 * dim;


    static size_t calc(const std::vector<plb::plint>& allGlobalIndices,
                       const size_t startIndex)
    {
        return detail::DimSizeCalcHelper<dim>::calc(allGlobalIndices, startIndex);
    }

};

typedef DimSizeCalc<2> DimSizeCalc2D;
typedef DimSizeCalc<3> DimSizeCalc3D;


}

#endif // DIMSIZECALC_H
