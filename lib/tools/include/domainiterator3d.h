#ifndef DOMAINITERATOR3D_H
#define DOMAINITERATOR3D_H

#include "domainiterator.h"

#include <palabos/core/geometry3D.h>

#include "palabosdataaccess3d.h"

namespace IPP
{
namespace DataAccess
{

template<>
class DomainIterator<3>
{
public:
    DomainIterator(const plb::Box3D& domain, const plb::Dot3D& currPos)
        : domain(domain),
          currPos(currPos)
    {

    }

    DomainIterator& operator++()
    {
        ++currPos.z;
        if(currPos.z > domain.z1)
        {
            // end of column reached
            currPos.z = domain.z0;
            ++currPos.y;

            if(currPos.y > domain.y1)
            {
                // end of column reached
                currPos.y = domain.y0;
                ++currPos.x;
            }
        }

        return *this;
    }

    DomainIterator operator++(int)
    {
        DomainIterator tmp(*this);
        operator ++();
        return tmp;
    }

    const plb::Dot3D& operator*() const
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
    const plb::Box3D& domain;
    plb::Dot3D currPos;
};

inline DomainIterator<3> begin(const plb::Box3D& domain)
{
    return DomainIterator<3>(domain, plb::Dot3D(domain.x0, domain.y0, domain.z0));
}

inline DomainIterator<3> end(const plb::Box3D& domain)
{
    return DomainIterator<3>(domain, plb::Dot3D(domain.x1+1, domain.y1+1, domain.z1+1));
}

}
}

#endif // DOMAINITERATOR3D_H
