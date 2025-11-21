#ifndef SETNTENSORFROMARRAY_H
#define SETNTENSORFROMARRAY_H

#include <palabos/core/geometry2D.h>
#include <palabos/core/geometry3D.h>
#include "ippvector.h"

namespace IPP
{

template<typename T, typename ArrayViewVec, size_t dim>
class SetNTensorFromArray
{
public:
    SetNTensorFromArray(const ArrayViewVec& scalarArrays,
                        const IPPVector3DLong& location,
                        const size_t tensorDim)
        : m_scalarArrays(scalarArrays)
        , m_location(location)
        , m_tensorDim(tensorDim)
    {

    }

    void operator()(const plb::Dot2D& pos, T* arr) const
    {
        (*this)(pos.x, pos.y, arr);
    }

    void operator()(plb::plint x, plb::plint y, T* arr) const
    {
        for(size_t i = 0; i < m_tensorDim; ++i)
        {
            const auto& arrView = m_scalarArrays[i];
            const auto& scalar = arrView[x - m_location[0]][y - m_location[1]];
            arr[i] = scalar;
        }
    }

    void operator()(const plb::Dot3D& pos, T* arr) const
    {
        (*this)(pos.x, pos.y, pos.z, arr);
    }

    void operator()(plb::plint x, plb::plint y, plb::plint z, T* arr) const
    {
        for(size_t i = 0; i < m_tensorDim; ++i)
        {
            const auto& arrView = m_scalarArrays[i];
            const auto& scalar = arrView[x - m_location[0]][y - m_location[1]][z - m_location[2]];
            arr[i] = scalar;
        }
    }

private:
    const ArrayViewVec& m_scalarArrays;
    const IPPVector3DLong& m_location;
    const size_t m_tensorDim;
};


}


#endif // SETNTENSORFROMARRAY_H
