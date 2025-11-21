#ifndef MULTIDIMARRAY_H
#define MULTIDIMARRAY_H

#include <vector>
#include <cassert>

#include "ippexception.h"

namespace boost
{
namespace serialization
{
class access;
}
}


namespace IPP
{

template<typename T>
class MultiDimArray
{
public:
    typedef std::vector<T> Element;
    typedef std::vector<size_t> IndexCoord;

    MultiDimArray();

    MultiDimArray(const IndexCoord& dimSizes, const size_t elementDim);

    const IndexCoord& getDimensions() const;
    size_t getElementSize() const;

    bool contains(const IndexCoord& coord) const;


    void setValue(const IndexCoord& coord, const Element &element);
    Element& evaluate(const IndexCoord& coord);
    const Element& evaluate(const IndexCoord& coord) const;

private:
    size_t calcIndex(const IndexCoord& coord) const;


    std::vector<Element> m_array;
    IndexCoord dimSizes;
    size_t elementDim;
    size_t dim;
    size_t totalSize;
    IndexCoord sizeConvFactors;


    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & m_array;
        ar & dimSizes;
        ar & elementDim;
        ar & dim;
        ar & totalSize;
        ar & sizeConvFactors;        
    }

};


template<typename T>
MultiDimArray<T>::MultiDimArray()
{

}

template<typename T>
MultiDimArray<T>::MultiDimArray(const std::vector<size_t> &dimSizes,
                                const size_t elementDim)
    : dimSizes(dimSizes)
    , elementDim(elementDim)
    , dim(dimSizes.size())
    , totalSize(0)
{
    IPPCheck::assertCheck(dimSizes.empty() == false);
    IPPCheck::assertCheck(elementDim > 0);

    totalSize = dimSizes.front();
    for(size_t i = 1; i < dimSizes.size(); ++i)
    {
        totalSize *= dimSizes[i];
    }

    sizeConvFactors.resize(dimSizes.size());
    size_t tmp = 1;
    for(size_t i = dimSizes.size(); i --> 0 ;)
    {
        sizeConvFactors[i] = tmp;
        tmp *= dimSizes[i];
    }

    const Element emptyElem(elementDim);
    m_array.resize(totalSize, emptyElem);
}

template<typename T>
const std::vector<size_t> &MultiDimArray<T>::getDimensions() const
{
    return dimSizes;
}

template<typename T>
size_t MultiDimArray<T>::getElementSize() const
{
    return elementDim;
}

template<typename T>
bool MultiDimArray<T>::contains(const IndexCoord& coord) const
{
    assert(dimSizes.size() == coord.size());

    for(size_t i = 0; i < coord.size(); ++i)
    {
        if(dimSizes[i] <= coord[i])
        {
            return false;
        }
    }

    return true;
}

template<typename T>
void MultiDimArray<T>::setValue(const std::vector<size_t> &coord, const Element& element)
{
    const size_t index = calcIndex(coord);
    m_array[index] = element;
}

template<typename T>
typename MultiDimArray<T>::Element &MultiDimArray<T>::evaluate(const IndexCoord &coord)
{
    const size_t index = calcIndex(coord);
    return m_array[index];
}

template<typename T>
const typename MultiDimArray<T>::Element &MultiDimArray<T>::evaluate(const IndexCoord &coord) const
{
    const size_t index = calcIndex(coord);
    return m_array[index];
}

template<typename T>
size_t MultiDimArray<T>::calcIndex(const IndexCoord& coord) const
{
    assert(coord.size() == dimSizes.size());
    // clamp coord to boudaries ?
    assert(contains(coord));
    assert(sizeConvFactors.size() == coord.size());

    // index = x * m_ny + y;
    // index = x * m_nyz + y * m_nz + z;
    // etc....
    size_t index = 0;
    for(size_t i = 0; i < coord.size(); ++i)
    {
        index += coord[i] * sizeConvFactors[i];
    }

    assert(index < m_array.size());

    return index;
}


}


#endif // MULTIDIMARRAY_H
