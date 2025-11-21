#ifndef ARRAY_VIEW_STRM_H
#define ARRAY_VIEW_STRM_H

#include <iostream>

namespace av
{

template <size_t Rank> class bounds;
template <size_t Rank> class offset;

}

template<size_t rank>
inline std::ostream& operator<<(std::ostream& ostr, const av::bounds<rank>& b)
{
    ostr << "(" << b[0];
    for(size_t iRank = 2; iRank <= rank; ++iRank)
    {
        ostr << ","<< b[iRank-1];
    }

    ostr << ")";

    return ostr;
}


template<size_t rank>
inline std::ostream& operator<<(std::ostream& ostr, const av::offset<rank>& b)
{
    ostr << "(" << b[0];
    for(size_t iRank = 2; iRank <= rank; ++iRank)
    {
        ostr << ","<< b[iRank-1];
    }

    ostr << ")";

    return ostr;
}

#endif // ARRAY_VIEW_STRM_H
