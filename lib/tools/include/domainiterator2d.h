#ifndef DOMAINITERATOR2D_H
#define DOMAINITERATOR2D_H

#include "domainiterator.h"
#include <palabos/core/geometry2D.h>

#include "palabosdataaccess2d.h"

namespace IPP
{

namespace DataAccess
{

template<>
class DomainIterator<2>
{
public:
    DomainIterator(const plb::Box2D& domain, const plb::Dot2D& currPos)
        : domain(domain),
          currPos(currPos)
    {

    }

    DomainIterator& operator++()
    {
        ++currPos.y;
        if(currPos.y > domain.y1)
        {
            // end of column reached
            currPos.y = domain.y0;
            ++currPos.x;
        }

        return *this;
    }

    DomainIterator operator++(int)
    {
        DomainIterator tmp(*this);
        operator ++();
        return tmp;
    }

    const plb::Dot2D& operator*() const
    {
        return currPos;
    }

    bool operator!=(const DomainIterator& other) const
    {
        return !isEqual(domain, other.domain)
                || !isEqual(currPos, other.currPos);
    }

    bool operator<(const DomainIterator& other) const
    {
        return isLower(currPos, other.currPos);
    }

private:
    const plb::Box2D& domain;
    plb::Dot2D currPos;
};

inline DomainIterator<2> begin(const plb::Box2D& domain)
{
    return DomainIterator<2>(domain, plb::Dot2D(domain.x0, domain.y0));
}

inline DomainIterator<2> end(const plb::Box2D& domain)
{
    return DomainIterator<2>(domain, plb::Dot2D(domain.x1+1, domain.y1+1));
}


}
}
#endif // DOMAINITERATOR2D_H
