#ifndef PALABOSIO_H
#define PALABOSIO_H

#include <iostream>

#include <palabos/core/geometry2D.h>
#include <palabos/core/geometry3D.h>

namespace IPP
{

template<typename T, size_t size>
inline std::ostream& operator<<(std::ostream& ostr, const plb::Array<T,size>& arr)
{
    ostr << "(";
    for(size_t i = 0; i < size; ++i)
    {
        ostr << arr[i];

        if(i+1 != size)
        {
            ostr << ",";
        }
    }

    ostr << ")";

    return ostr;
}

inline std::ostream& operator<<(std::ostream& ostr, const plb::Dot2D& dot)
{
    ostr << "(" << dot.x << ","<< dot.y << ")";
    return ostr;
}

inline std::ostream& operator<<(std::ostream& ostr, const plb::Dot3D& dot)
{
    ostr << "(" << dot.x << ","<< dot.y << "," << dot.z << ")";
    return ostr;
}

inline std::ostream& operator<<(std::ostream& ostr, const plb::Box2D& box)
{
    ostr << box.x0
       << "\t "<< box.x1
       << "\t" << box.y0
       << "\t "<< box.y1;
    return ostr;
}

inline std::ostream& operator<<(std::ostream& ostr, const plb::Box3D& box)
{
    ostr << box.x0
       << "\t "<< box.x1
       << "\t" << box.y0
       << "\t "<< box.y1
       << "\t" << box.z0
       << "\t "<< box.z1;
    return ostr;
}

inline std::ostream& operator<<(std::ostream& ostr, const plb::DotList2D& dots)
{
    for(plb::plint i = 0; i < dots.getN(); ++i)
    {
        const auto& dot = dots.getDot(i);
        ostr << dot;

        if(i+1 < dots.getN())
        {
            ostr << ",";
        }
    }

    return ostr;
}

inline std::ostream& operator<<(std::ostream& ostr, const plb::DotList3D& dots)
{
    for(plb::plint i = 0; i < dots.getN(); ++i)
    {
        const auto& dot = dots.getDot(i);
        ostr << dot << ",";

        if(i+1 < dots.getN())
        {
            ostr << ",";
        }
    }

    return ostr;
}

}

#endif // PALABOSIO_H
