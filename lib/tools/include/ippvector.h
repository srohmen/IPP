#ifndef IPPVECTOR_H
#define IPPVECTOR_H

#include <iterator>
#include <array>
#include <boost/serialization/array.hpp>
#include <initializer_list>

namespace IPP
{

template<typename T, size_t dim>
class Array
{
public:

    typedef T value_type;

    Array() = default;

    Array(std::initializer_list<T> il)
    {
        std::copy(il.begin(), il.end(), m_data.begin());
    }


    void fill(const T& v)
    {
        m_data.fill(v);
    }

    size_t size() const
    {
        return m_data.size();
    }

    inline size_t calcLength2() const
    {
        T l = m_data[0];
        for(size_t i = 1; i < dim; ++i)
        {
            l *= m_data[i];
        }
        return l;
    }



    inline T& operator[](const size_t i)
    {
        return m_data[i];
    }

    inline const T& operator[](const size_t i) const
    {
        return m_data[i];
    }

    template<typename U>
    inline void operator+=(const Array<U, dim>& o)
    {
        for(size_t i = 0; i < dim; ++i)
        {
            m_data[i] += o[i];
        }
    }

    template<typename U>
    inline void operator-=(const Array<U, dim>& o)
    {
        for(size_t i = 0; i < dim; ++i)
        {
            m_data[i] -= o[i];
        }
    }

    template<typename U>
    inline bool operator==(const Array<U, dim>& o) const
    {
        return m_data == o.m_data;
    }

    template<typename U>
    inline bool operator!=(const Array<U, dim>& o) const
    {
        return (*this == o) == false;
    }

    template<typename U>
    inline bool operator<(const Array<U, dim>& o) const
    {
        for(size_t i = 0; i < dim; ++i)
        {
            if(m_data[i] >= o[i])
            {
                return false;
            }
        }
        return true;
    }

    template<typename U>
    inline bool operator<=(const Array<U, dim>& o) const
    {
        for(size_t i = 0; i < dim; ++i)
        {
            if(m_data[i] > o[i])
            {
                return false;
            }
        }
        return true;
    }

    template<typename U>
    inline bool operator>(const Array<U, dim>& o) const
    {
        return (*this <= o) == false;
    }

    template<typename U>
    inline bool operator>=(const Array<U, dim>& o) const
    {
        return (*this < o) == false;
    }

    friend inline std::ostream& operator<<(std::ostream& o, const Array& arr)
    {
        std::copy(arr.m_data.cbegin(), arr.m_data.cend(), std::ostream_iterator<T>(o, " "));
        return o;
    }


private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
        ar & boost::serialization::make_nvp(
            "elems",
            *static_cast<T (*)[dim]>(static_cast<void *>(m_data.data()))
            );
        // ar & m_data;
    }

    std::array<T,dim> m_data;
};


template<typename T, typename U, size_t dim>
inline Array<T, dim> operator+(const Array<T, dim>& a,
                               const Array<U, dim>& b)
{
    Array<T, dim> tmp = a;
    tmp += b;
    return tmp;
}

template<typename T, typename U, size_t dim>
inline Array<T, dim> operator-(const Array<T, dim>& a,
                               const Array<U, dim>& b)
{
    Array<T, dim> tmp = a;
    tmp -= b;
    return tmp;
}


template<typename T, typename U, size_t N>
inline Array<T, N> operator *(const Array<T, N>& arr, const U& scalar)
{
    Array<T, N> tmp = arr;
    for(size_t i = 0; i < N; ++i)
    {
        tmp[i] *= scalar;
    }
    return tmp;
}

template <typename T, typename A, std::size_t N>
inline Array<T, N> convert(const Array<A, N>& src)
{
    Array<T, N> dst;
    for(size_t i = 0; i < N; ++i)
    {
        dst[i] = src[i];
    }
    return dst;
}


typedef Array<double, 2> IPPVector2D;
typedef Array<double, 3> IPPVector3D;

typedef Array<int, 2> IPPVector2DInt;
typedef Array<int, 3> IPPVector3DInt;

typedef Array<long, 2> IPPVector2DLong;
typedef Array<long, 3> IPPVector3DLong;

typedef Array<size_t, 2> IPPVector2DUInt;
typedef Array<size_t, 3> IPPVector3DUInt;

} // end of namespace IPP

#endif // IPPVECTOR_H
