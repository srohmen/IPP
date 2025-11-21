#ifndef DOTFINDNEIGHBORS_H
#define DOTFINDNEIGHBORS_H

#include <palabos/core/geometry2D.h>
#include <palabos/core/geometry3D.h>

namespace DotFindNeighbors
{

inline plb::DotList2D find(const plb::Dot2D &pos)
{
    plb::DotList2D list;

    static const plb::Dot2D dx(1, 0);
    list.addDot(pos - dx);
    list.addDot(pos + dx);

    static const plb::Dot2D dy(0, 1);
    list.addDot(pos - dy);
    list.addDot(pos + dy);

    return list;
}

inline plb::DotList3D find(const plb::Dot3D &pos)
{
    plb::DotList3D list;

    static const plb::Dot3D dx(1, 0, 0);
    list.addDot(pos - dx);
    list.addDot(pos + dx);

    static const plb::Dot3D dy(0, 1, 0);
    list.addDot(pos - dy);
    list.addDot(pos + dy);

    static const plb::Dot3D dz(0, 0, 1);
    list.addDot(pos - dz);
    list.addDot(pos + dz);

    return list;
}

}

#endif // DOTFINDNEIGHBORS_H
