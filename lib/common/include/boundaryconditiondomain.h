#ifndef BOUNDARYCONDITIONDOMAIN_H
#define BOUNDARYCONDITIONDOMAIN_H

#include "ippbox.h"

namespace IPP
{

struct BoundaryConditionDomain
{
    enum BoundaryPosition
    {
        BP_left,    // low x
        BP_right,   // high x
        BP_bottom,  // low y
        BP_top,     // high y
        BP_back,    // low z
        BP_front    // high z
    };

    BoundaryPosition position;
    IPPBox3DInt range;
};

}

#endif // BOUNDARYCONDITIONDOMAIN_H
