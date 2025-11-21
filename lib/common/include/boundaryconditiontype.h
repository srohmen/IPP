#ifndef BOUNDARYCONDITIONTYPE_H
#define BOUNDARYCONDITIONTYPE_H

namespace IPP
{

enum BoundaryConditionType
{
    BCT_Closed,     // bounce back (density -> neumann)
    BCT_DensityDirichlet,  // density -> dirichlet (density -> constant)
    BCT_VelocityDirichlet,    // velocity -> dirichlet (velocity -> constant)
    BCT_DensityCauchy,     // density -> dirichlet + velocity -> dirichlet (density -> neumann)
    BCT_VelocityOutflow // velocity -> neumann (open outflow)
};

} // end of namespace

#endif // BOUNDARYCONDITIONTYPE_H
