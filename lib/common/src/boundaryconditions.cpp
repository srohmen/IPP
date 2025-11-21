#include "boundaryconditions.h"

namespace IPP
{


BoundaryConditions::BoundaryConditions()
{

}

const std::vector<int>& BoundaryConditions::getPeriodicBoundaryDimensions() const
{
    return periodicDims;
}

bool BoundaryConditions::contains(const std::string& elementName) const
{
    return diffusiveBC.find(elementName) != diffusiveBC.end();
}


const std::vector<AbstractBoundaryConditions::BoundaryConditionData>&
BoundaryConditions::getAdvectiveBoundaryConditions() const
{
    return advectiveBC;
}

const AbstractBoundaryConditions::BCVec&
BoundaryConditions::getDiffusiveBoundaryConditions(const std::string& elementName) const
{
    return diffusiveBC.at(elementName);
}

}
