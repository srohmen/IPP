#ifndef DERIV3D_H
#define DERIV3D_H

#include <cmath>
#include <limits>

namespace IPP
{

template<typename T>
struct Deriv3D
{
private:
    static const T highInt;

public:
    Deriv3D()
        : dx(highInt),
          dy(highInt),
          dz(highInt)
    {

    }

    Deriv3D(const T dx, const T dy, const T dz)
        : dx(dx),
          dy(dy),
          dz(dz)
    {

    }

    inline Deriv3D& operator+=(const Deriv3D& other)
    {
        dx += other.dx;
        dy += other.dy;
        dz += other.dz;
        return *this;
    }

    inline Deriv3D operator+(const Deriv3D& other) const
    {
        Deriv3D result = *this;
        result += other;
        return result;
    }

    inline T length2() const
    {
        return dx*dx + dy*dy + dz*dz;
    }

    T dx, dy, dz;
};


template<typename T>
const T Deriv3D<T>::highInt = std::sqrt(std::numeric_limits<T>::max()) / 2;


}


#endif // DERIV3D_H
