#ifndef ARRAYVIEWITERATORIMPL_H
#define ARRAYVIEWITERATORIMPL_H

#include "array_view.h"

namespace IPP
{


class ArrayViewIteratorImpl
{
public:
    typedef av::strided_array_view<const double, 2> array_view;
    typedef av::bounds_iterator<2> bounds_iterator;

    ArrayViewIteratorImpl(const array_view& view, bounds_iterator it)
        : m_view(view),
          m_it(it)
    {

    }

    ArrayViewIteratorImpl &operator ++()
    {
        ++m_it;
        return *this;
    }

    bool operator ==(const ArrayViewIteratorImpl &other) const
    {
        return m_it == other.m_it;
    }

    bool operator !=(const ArrayViewIteratorImpl &other) const
    {
        return m_it != other.m_it;
    }

    const double &operator *()
    {
        return m_view[*m_it];
    }

private:
    const array_view& m_view;
    bounds_iterator m_it;
};


}

#endif // ARRAYVIEWITERATORIMPL_H
