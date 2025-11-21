#ifndef ARRAYVIEWDATAACCESS_H
#define ARRAYVIEWDATAACCESS_H

#include <palabos/core/geometry2D.h>
#include <palabos/core/geometry3D.h>
#include "array_view.h"

namespace IPP
{

namespace DataAccess
{

template<typename T>
auto& get(av::array_view<T, 2>& container, const plb::Dot2D& pos)
{
    return container[pos.x][pos.y];
}

template<typename T>
const auto& get(const av::array_view<T, 2>& container, const plb::Dot2D& pos)
{
    return container[pos.x][pos.y];
}

template<typename T, typename I>
const auto& get(const av::array_view<T, 2>& container, const I x, const I y, const I /*z*/)
{
    return container[x][y];
}




template<typename T>
auto& get(av::array_view<T,3>& container, const plb::Dot3D& pos)
{
    return container[pos.x][pos.y][pos.z];
}

template<typename T>
const auto& get(const av::array_view<T,3>& container, const plb::Dot3D& pos)
{
    return container[pos.x][pos.y][pos.z];
}

template<typename T, typename I>
const auto& get(const av::array_view<T, 3>& container, const I x, const I y, const I z)
{
    return container[x][y][z];
}

}


}


#endif // ARRAYVIEWDATAACCESS_H
