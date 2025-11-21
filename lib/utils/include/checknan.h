#ifndef CHECKNAN_H
#define CHECKNAN_H

#include <cmath>

namespace IPP
{

template<typename T>
class CheckNaN
{
public:
    CheckNaN(bool& hasNan)
        : m_hasNaN(hasNan)
    {

    }

    T operator()(const T& input)
    {
        if(std::isnan(input) || std::isnan(-input))
        {
            m_hasNaN = true;
            assert(false);
        }

        return input;
    }

private:
    bool& m_hasNaN;
};


}

#endif // CHECKNAN_H
