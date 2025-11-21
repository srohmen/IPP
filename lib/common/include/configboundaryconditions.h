#ifndef CONFIGBOUNDARYCONDITIONS_H
#define CONFIGBOUNDARYCONDITIONS_H

#include <memory>
#include <vector>

#include "boundaryconditiondata.h"

namespace IPP
{

class Composition;

class ConfigBoundaryConditions
{
public:
    struct AdvectiveBoundaryConditionData : public BoundaryConditionData
    {
        double density = 0.0;
        IPPVector3D velocity = {{0.0,0.0,0.0}};
    };
    typedef std::vector<AdvectiveBoundaryConditionData> AdvectiveBCVec;

    struct DiffusiveBoundaryConditionData : public BoundaryConditionData
    {
        const Composition* composition;
    };
    typedef std::vector<DiffusiveBoundaryConditionData> DiffusiveBCVec;

    std::vector<int> periodicDims;
    AdvectiveBCVec advectiveBC;
    DiffusiveBCVec diffusiveBC;
};

typedef std::shared_ptr<ConfigBoundaryConditions> ConfigBoundaryConditionsPtr;

} // end of namespace

#endif // CONFIGBOUNDARYCONDITIONS_H

