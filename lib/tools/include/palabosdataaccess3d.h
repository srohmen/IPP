#ifndef PALABOSDATAACCESS3D_H
#define PALABOSDATAACCESS3D_H

#include <palabos/core/geometry3D.h>
#include <palabos/atomicBlock/dataField3D.h>


namespace IPP
{

namespace DataAccess
{

template<typename Container>
auto& get(Container& container, const plb::Dot3D& pos)
{
    return container.get(pos.x, pos.y, pos.z);
}

template<typename Container>
const auto& get(const Container& container, const plb::Dot3D& pos)
{
    return container.get(pos.x, pos.y, pos.z);
}

template<typename T>
T* get(plb::NTensorField3D<T>& container, const plb::Dot3D& pos)
{
    return container.get(pos.x, pos.y, pos.z);
}

template<typename T>
const T* get(const plb::NTensorField3D<T>& container, const plb::Dot3D& pos)
{
    return container.get(pos.x, pos.y, pos.z);
}


inline bool isEqual(const plb::Dot3D& a, const plb::Dot3D& b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool isEqual(const plb::Box3D& a, const plb::Box3D& b)
{
    return a.x0 == b.x0
            && a.x1 == b.x1
            && a.y0 == b.y0
            && a.y1 == b.y1
            && a.z0 == b.z0
            && a.z1 == b.z1;
}

inline bool isLower(const plb::Dot3D& a, const plb::Dot3D& b)
{
    return a.z < b.z && a.y < b.y && a.x < b.x;
}


}

}

#endif // PALABOSDATAACCESS3D_H
