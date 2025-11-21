#include "nonreturningbcgen.h"

#include "configboundaryconditions.h"

namespace IPP
{


NonReturningBCGen::NonReturningBCGen()
{

}

ConfigBoundaryConditions* NonReturningBCGen::generate() const
{
    return new ConfigBoundaryConditions;
}


}
