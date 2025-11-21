#ifndef LBGEOCHEMMATHTOOLS_H
#define LBGEOCHEMMATHTOOLS_H

namespace IPP
{

namespace MathTools
{

template<typename T, typename U>
T wrapToRange(const T& toWrap, const U& maxVal)
{
    T mod = toWrap % maxVal;

    if(mod < 0)
    {
        mod = maxVal + mod;
    }

    return mod;
}

}

}


#endif // LBGEOCHEMMATHTOOLS_H
