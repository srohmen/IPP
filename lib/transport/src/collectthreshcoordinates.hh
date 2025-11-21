#ifndef COLLECTTHRESHCOORDINATES_HH
#define COLLECTTHRESHCOORDINATES_HH

#include "collectthreshcoordinates.h"

#include "domainiterator.h"
#include "palabosdataaccess.h"

namespace IPP
{


template<typename T, size_t dim>
CollectThreshCoordinates<T,dim>::CollectThreshCoordinates(const T &thresh,
                                                        DotList &above, DotList &below)
    : m_thresh(thresh)
    , m_above(above)
    , m_below(below)
{

}

template<typename T, size_t dim>
void CollectThreshCoordinates<T, dim>::getTypeOfModification(std::vector<plb::modif::ModifT> &modified) const
{
    modified[0] = plb::modif::nothing;
}

template<typename T, size_t dim>
plb::BlockDomain::DomainT CollectThreshCoordinates<T,dim>::appliesTo() const
{
    return plb::BlockDomain::bulkAndEnvelope;
}

template<typename T, size_t dim>
CollectThreshCoordinates<T, dim> *CollectThreshCoordinates<T, dim>::clone() const
{
    return new CollectThreshCoordinates<T,dim>(*this);
}

template<typename T, size_t dim>
void CollectThreshCoordinates<T,dim>::process(Box domain, InputField& field)
{
    const auto offset = field.getLocation();

    const auto end = DataAccess::end(domain);
    for(auto it = DataAccess::begin(domain); it < end; ++it)
    {
        const auto& pos = *it;
        const T& value = DataAccess::get(field, pos);


        const auto posGlobal = pos + offset;

        if(value < m_thresh)
        {            
            m_below.addDot(posGlobal);
        }
        else
        {
            m_above.addDot(posGlobal);
        }
    }
}

}

#endif // COLLECTTHRESHCOORDINATES_HH
