#ifndef MULTIDIMNEARESTNEIGHBORINTERPOLATION_H
#define MULTIDIMNEARESTNEIGHBORINTERPOLATION_H

#include "multidiminterpolation.h"


namespace IPP
{


struct NearestNeighborInterpolation
{
    template<typename Element, typename CoordScalar>
    static void interpolate(Element &a,
                            Element &b,
                            const CoordScalar& interp,
                            Element& result)
    {
        assert(a.size() == b.size());

        if(interp < 0.5)
        {
            result.swap(a);
        }
        else
        {
            result.swap(b);
        }

    }
};


template<typename ElementScalar, typename CoordScalar = ElementScalar>
using MultiDimNearestNeighborInterpolation =
MultiDimInterpolation<NearestNeighborInterpolation, ElementScalar, CoordScalar>;


} // end of namespace IPP


#endif // MULTIDIMNEARESTNEIGHBORINTERPOLATION_H
