#ifndef DERIV2D_H
#define DERIV2D_H

#include <cmath>
#include <limits>

namespace IPP
{

template<typename T>
struct Deriv2D
{
private:
    static const T highInt;

public:
    Deriv2D()
        : dx(highInt),
          dy(highInt)
    {

    }

    Deriv2D(const T dx, const T dy)
        : dx(dx),
          dy(dy)
    {

    }

    T length2() const
    {
        return dx*dx + dy*dy;
    }


    T dx, dy;
};


template<typename T>
const T Deriv2D<T>::highInt = std::sqrt(std::numeric_limits<T>::max()) / 2;


}

#endif // DERIV2D_H
