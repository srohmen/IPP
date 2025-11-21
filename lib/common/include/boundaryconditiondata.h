#ifndef BOUNDARYCONDITIONDATA_H
#define BOUNDARYCONDITIONDATA_H

#include "boundaryconditiontype.h"
#include "boundaryconditiondomain.h"

namespace IPP
{

struct BoundaryConditionData
{
    BoundaryConditionDomain domain;
    BoundaryConditionType type = IPP::BCT_Closed;
};

}

#endif // BOUNDARYCONDITIONDATA_H
