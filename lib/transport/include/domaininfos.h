#ifndef DOMAININFOS_H
#define DOMAININFOS_H

#include <unordered_set>

#include "boundaryconditiondomain.h"
#include "ippvector.h"

namespace IPP
{

struct DomainInfos
{
    // enums need special hashing...
    typedef std::unordered_set<BoundaryConditionDomain::BoundaryPosition, std::hash<int>> BoundaryEnvelops;

    BoundaryEnvelops boundaryEnvelops;
    IPPVector3DInt totalSize;
    IPPVector3DInt origin;
};


}

#endif // DOMAININFOS_H
