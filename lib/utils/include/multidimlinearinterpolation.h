#ifndef MULTIDIMLINEARINTERPOLATION_H
#define MULTIDIMLINEARINTERPOLATION_H

#include "multidiminterpolation.h"


namespace IPP
{

struct LinearInterpolation
{
    template<typename Element, typename CoordScalar>
    static void interpolate(Element &a,
                            Element &b,
                            const CoordScalar& interp,
                            Element& result)
    {
        assert(a.size() == b.size());

        const CoordScalar inv = 1.0 - interp;

        result.swap(a);
        for(size_t i = 0; i < result.size(); ++i)
        {
            result[i] *= inv;
            result[i] += interp * b[i];
        }
    }
};

template<typename ElementScalar, typename CoordScalar = ElementScalar>
using MultiDimLinearInterpolation =
MultiDimInterpolation<LinearInterpolation, ElementScalar, CoordScalar>;


} // end of namespace IPP



#endif // MULTIDIMLINEARINTERPOLATION_H
