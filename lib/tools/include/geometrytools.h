#ifndef GEOMETRYTOOLS_H
#define GEOMETRYTOOLS_H

#include <vector>
#include <cstddef>

#include "ippvector.h"

namespace IPP
{

namespace GeometryTools
{

template<typename Vector>
bool isWithin(const Vector& pos, const Vector& bounds)
{
    for(size_t dim = 0; dim < pos.size(); ++dim)
    {
        if(pos[dim] < 0 || pos[dim] >= bounds[dim])
        {
            return false;
        }
    }

    return true;
}

void transpose(std::vector<double>& arr2D,
               const size_t nx, const size_t ny);
}

}

#endif // GEOMETRYTOOLS_H
