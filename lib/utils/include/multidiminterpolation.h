#ifndef MULTIDIMINTERPOLATION_H
#define MULTIDIMINTERPOLATION_H

#include "multidimarray.h"
#include <cmath>
#include <limits>

#include "stringutils.h"

namespace boost
{
namespace serialization
{
class access;
}
}


namespace IPP
{



template<typename Interpolation, typename ElementScalar, typename CoordScalar = ElementScalar>
class MultiDimInterpolation
{
public:
    static_assert(std::is_floating_point<CoordScalar>::value,
                  "Only floating point type is supported");

    typedef typename MultiDimArray<ElementScalar>::Element Element;
    typedef std::vector<CoordScalar> Coord;
    typedef std::vector<size_t> IndexCoord;

    MultiDimInterpolation() = default;

    template<typename... Args>
    MultiDimInterpolation(Args&&... args);


    const IndexCoord& getDimensions() const;
    const Coord& getOrigin() const;
    const Coord& getSpacing() const;
    size_t getElementSize() const;

    bool contains(const Coord &coord) const;


    void setValue(const IndexCoord& coord, const Element &element);

    void evaluate(const Coord& coord, Element& result) const;
    Element& evaluate(const IndexCoord& coord);
    const Element& evaluate(const IndexCoord& coord) const;

private:
    // TODO: move to static compile time function
    void init();
    void initNeigh();

    void findIndexCoord(const Coord& coord, IndexCoord& indexCoord, Coord& interp) const;


    static void interpolate(std::vector<Element> &elements, const Coord& interp,
                            const size_t currDim, Element& result);

    static inline void interpolateElement(Element& a, Element& b,
                                          const CoordScalar& interp, Element& result);

    void checkContains(const Coord& coord) const;


    struct Data
    {
        Data() = default;


        Data(const MultiDimArray<ElementScalar> &array,
             const Coord &spacing,
             const Coord &origin = {})
            : array(array)
            , spacing(spacing)
            , origin(origin)
        {
            IPPCheck::assertCheck(array.getDimensions().size() == spacing.size());
            init();
        }

        Data(const std::vector<size_t> &dimSizes,
             const size_t elementDim,
             const Coord& spacing,
             const Coord &origin = {})
            : array(dimSizes, elementDim)
            , spacing(spacing)
            , origin(origin)
        {
            IPPCheck::assertCheck(array.getDimensions().size() == spacing.size());
            init();
        }


        void init()
        {
            if(origin.size() != spacing.size())
            {
                origin.resize(spacing.size(), 0.0);
            }

            initNeigh();
        }



        void initNeigh()
        {
            // TODO: move to static compile time function

            // finding neighbours with help of binary numbers
            // sorted by x, y, z, etc. (increment in big endian fashion)
            const size_t nDim = array.getDimensions().size();
            const size_t nNeigh = std::exp2(nDim);
            neighOffset.resize(nNeigh);
            for(size_t iNeigh = 0; iNeigh < nNeigh; ++iNeigh)
            {
                IndexCoord& offset = neighOffset[iNeigh];
                offset.resize(nDim);

                for(size_t d = 0; d < nDim; ++d)
                {
                    const size_t mask = (1 << d);
                    const bool digit = iNeigh & mask;
                    offset[d] = digit;
                }
            }
        }


        MultiDimArray<ElementScalar> array;
        Coord spacing;
        Coord origin;


        // TODO: move to static compile time function
        std::vector<IndexCoord> neighOffset;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & array;
            ar & spacing;
            ar & origin;

            // TODO: remove from serialization if determined at compile time
            ar & neighOffset;
        }

    };


    Data m_data;



    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & m_data;
    }


};


template<typename Interpolation, typename ElementScalar, typename CoordScalar>
template<typename... Args>
MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::MultiDimInterpolation(Args&&... args)
    : m_data(std::forward<Args>(args)...)

{

}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
const typename MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::IndexCoord&
MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::getDimensions() const
{
    return m_data.array.getDimensions();
}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
const typename MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::Coord&
MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::getOrigin() const
{
    return m_data.origin;
}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
const typename MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::Coord &
MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::getSpacing() const
{
    return m_data.spacing;
}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
size_t MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::getElementSize() const
{
    return m_data.array.getElementSize();
}


template<typename Interpolation, typename ElementScalar, typename CoordScalar>
bool MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::contains(const Coord &coord) const
{
    assert(coord.size() == m_data.spacing.size());

    const IndexCoord& dimSizes = m_data.array.getDimensions();

    for(size_t i = 0; i < coord.size(); ++i)
    {
        const CoordScalar val = coord[i] - m_data.origin[i];
        const CoordScalar max = dimSizes[i] * m_data.spacing[i];
        if(val < 0.0 || val >= max)
        {
            return false;
        }
    }

    return true;

}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
void MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::setValue(const IndexCoord &coord,
                                                                                const Element &element)
{
    m_data.array.setValue(coord, element);
}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
void MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::checkContains(const Coord &coord) const
{
    if(contains(coord) == false)
    {
        std::string errStr = "linear interpolation does not contain coordinate: "
                             + StringUtils::toString(coord)
                             + " boundarys are:\n";

        const IndexCoord& dimSizes = m_data.array.getDimensions();
        for(size_t i = 0; i < coord.size(); ++i)
        {
            const CoordScalar min = m_data.origin[i];
            const CoordScalar max = dimSizes[i] * m_data.spacing[i] + m_data.origin[i];
            errStr += std::to_string(i) + ": [ " + std::to_string(min) + ", " + std::to_string(max) + " ]\n";
        }

        throw std::runtime_error(errStr);
    }
}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
void MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::interpolate(std::vector<Element> &elements,
                                                                                   const Coord& interp,
                                                                                   const size_t currDim,
                                                                                   Element& result)
{
    const CoordScalar& axisInterp = interp[currDim];

    const size_t nNewElements = elements.size() / 2;
    std::vector<Element> reducedElements(nNewElements);

    assert(elements.size() == 2 * nNewElements );


    for(size_t i = 0; i < nNewElements; ++i)
    {
        const size_t index = i * 2;
        Element& a = elements[index];
        Element& b = elements[index+1];

        assert(reducedElements.size() > i);
        Element& interpolated = reducedElements[i];
        Interpolation::interpolate(a, b, axisInterp, interpolated);
    }


    if(reducedElements.size() == 1)
    {
        Element& final = reducedElements.front();
        result.swap(final);
        return;
    }
    else
    {
        interpolate(reducedElements, interp, currDim+1, result);
    }

}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
void MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::evaluate(const Coord &coord,
                                                                                Element &result) const
{
    checkContains(coord);

    IndexCoord index;
    Coord interp;
    findIndexCoord(coord, index, interp);



    std::vector<Element> elements(m_data.neighOffset.size());
    for(size_t iNeigh = 0; iNeigh < m_data.neighOffset.size(); ++iNeigh)
    {
        const IndexCoord& localCoord = m_data.neighOffset[iNeigh];

        IndexCoord absIndex(index.size());
        for(size_t d = 0; d < index.size(); ++d)
        {
            static const CoordScalar tolerance = 1.0E-12;
            if(interp[d] <= tolerance)
            {
                // snap down
                absIndex[d] = index[d];
            }
            else
            {
                absIndex[d] = index[d] + localCoord[d];

            }


        }

        assert(m_data.array.contains(absIndex));
        const Element& neighElement = m_data.array.evaluate(absIndex);
        elements[iNeigh] = neighElement;
    }

    interpolate(elements, interp, 0, result);



}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
typename MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::Element&
MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::evaluate(const IndexCoord &coord)
{
    return m_data.array.evaluate(coord);
}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
const typename MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::Element&
MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::evaluate(const IndexCoord &coord) const
{
    return m_data.array.evaluate(coord);
}

template<typename Interpolation, typename ElementScalar, typename CoordScalar>
void MultiDimInterpolation<Interpolation, ElementScalar, CoordScalar>::findIndexCoord(const Coord& coord,
                                                                                      IndexCoord& indexCoord,
                                                                                      Coord& interp) const
{
    assert(coord.size() == m_data.spacing.size());

    indexCoord.resize(coord.size());
    interp.resize(coord.size());

    for(size_t i = 0; i < coord.size(); ++i)
    {
        const CoordScalar& val = coord[i] - m_data.origin[i];
        const CoordScalar& space = m_data.spacing[i];

        const CoordScalar diffToNode = std::fmod(val, space);
        const CoordScalar a = val - diffToNode + space * 0.5;
        const CoordScalar div = a / space;
        const size_t index = div;

        const CoordScalar t = diffToNode / space;

        assert(t >= 0.0);
        assert(t < 1.0);

        interp[i] = t;
        indexCoord[i] = index;
    }

    assert(m_data.array.contains(indexCoord));

}


} // end of namespace IPP

#endif // MULTIDIMINTERPOLATION_H
