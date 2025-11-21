#ifndef IPPBOX_H
#define IPPBOX_H

#include <iostream>
#include "ippvector.h"

namespace boost
{
namespace serialization
{
class access;
}
}


namespace IPP
{

template<typename VecType>
struct IPPBox
{
    typedef typename VecType::value_type value_type;

    IPPBox()
    {
    }

    IPPBox(const VecType& lower, const VecType& upper)
        : lower(lower),
          upper(upper)
    {

    }

    value_type& operator[](const size_t index)
    {
        const std::div_t dv = std::div(index, 2);
        if(dv.rem == 0)
        {
            return lower[dv.quot];
        }
        else
        {
            return upper[dv.quot];
        }
    }

    const value_type& operator[](const size_t index) const
    {
        const std::div_t dv = std::div(index, 2);
        if(dv.rem == 0)
        {
            return lower[dv.quot];
        }
        else
        {
            return upper[dv.quot];
        }
    }

    bool operator==(const IPPBox<VecType>& other) const
    {
        const bool result
                = lower == other.lower
                && upper == other.upper;
        return result;
    }

    bool operator!=(const IPPBox<VecType>& other) const
    {
        return !(*this == other);
    }


    value_type getDiagonalSize() const
    {
        const VecType diff = getSize();
        size_t size = diff[0];
        for(size_t i = 1; i < diff.size(); ++i)
        {
            size *= diff[i];
        }
        return size;
    }

    VecType getSize() const
    {
        VecType unitVec;
        unitVec.fill(1);

        const VecType diff = upper - lower + unitVec;
        return diff;
    }


    bool contains(const IPPBox& other) const
    {
        return other.lower >= lower && other.upper <= upper;
    }

    template<typename Vector>
    bool contains(const Vector& pos) const
    {
        return pos >= lower && pos <= upper;
    }



    VecType lower;
    VecType upper;

    friend std::ostream& operator <<(std::ostream& ostr, const IPPBox& box)
    {
        ostr << "( (" << box.lower << "), ("<< box.upper<< ") )";
        return ostr;
    }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & lower;
        ar & upper;
    }
};

template<typename VecType>
inline void intersection(const IPPBox<VecType>& a, const IPPBox<VecType>& b, IPPBox<VecType>& intersectBox)
{
    for(size_t i = 0; i < a.lower.size(); ++i)
    {
        intersectBox.lower[i] = std::max(a.lower[i], b.lower[i]);
    }

    for(size_t i = 0; i < a.lower.size(); ++i)
    {
        intersectBox.upper[i] = std::min(a.upper[i], b.upper[i]);
    }

}

template<typename VecType>
inline IPPBox<VecType> intersection(const IPPBox<VecType>& a, const IPPBox<VecType>& b)
{
    IPPBox<VecType> intersectBox;
    intersection(a, b, intersectBox);
    return intersectBox;
}

template<typename VecType>
inline bool isIntersecting(const IPPBox<VecType>& a, const IPPBox<VecType>& b)
{
    bool result = true;

    for(size_t i = 0; i < a.lower.size(); ++i)
    {
        result = result && (std::min(a.upper[i], b.upper[i]) >= std::max(a.lower[i], b.lower[i]));
    }

    return result;
}



typedef IPPBox<IPPVector2DInt> IPPBox2DInt;
typedef IPPBox<IPPVector3DInt> IPPBox3DInt;

typedef IPPBox<IPPVector2DLong> IPPBox2DLong;
typedef IPPBox<IPPVector3DLong> IPPBox3DLong;

typedef IPPBox<IPPVector2DUInt> IPPBox2DUInt;
typedef IPPBox<IPPVector3DUInt> IPPBox3DUInt;

}

#endif // IPPBOX_H
