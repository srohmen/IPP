#ifndef PALABOSDATAACCESS2D_H
#define PALABOSDATAACCESS2D_H

#include <palabos/core/geometry2D.h>
#include <palabos/atomicBlock/dataField2D.h>


namespace IPP
{

namespace DataAccess
{


template<typename Container>
auto& get(Container& container, const plb::Dot2D& pos)
{
    return container.get(pos.x, pos.y);
}

template<typename Container>
const auto& get(const Container& container, const plb::Dot2D& pos)
{
    return container.get(pos.x, pos.y);
}

template<typename T>
T* get(plb::NTensorField2D<T>& container, const plb::Dot2D& pos)
{
    return container.get(pos.x, pos.y);
}

template<typename T>
const T* get(const plb::NTensorField2D<T>& container, const plb::Dot2D& pos)
{
    return container.get(pos.x, pos.y);
}


inline bool isEqual(const plb::Dot2D& a, const plb::Dot2D& b)
{
    return a.x == b.x && a.y == b.y;
}

inline bool isEqual(const plb::Box2D& a, const plb::Box2D& b)
{
    return a.x0 == b.x0
            && a.x1 == b.x1
            && a.y0 == b.y0
            && a.y1 == b.y1;
}

inline bool isLower(const plb::Dot2D& a, const plb::Dot2D& b)
{
    return a.y < b.y && a.x < b.x;

}

}

}

#endif // PALABOSDATAACCESS2D_H
