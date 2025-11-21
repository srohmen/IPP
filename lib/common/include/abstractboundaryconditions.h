#ifndef ABSTRACTBOUNDARYCONDITIONS_H
#define ABSTRACTBOUNDARYCONDITIONS_H

#include "abstractboundaryconditionsfwd.h"
#include <vector>
#include <map>

#include "boundaryconditiontype.h"
#include "boundaryconditiondomain.h"

namespace IPP
{

class AbstractBoundaryConditions
{
public:

    struct BoundaryConditionData
    {
        BoundaryConditionDomain domain;
        BoundaryConditionType type = BCT_Closed;
        double density = 0.0;
        IPPVector3D velocity = {{0.0,0.0,0.0}};
    };

    typedef std::vector<BoundaryConditionData> BCVec;


    AbstractBoundaryConditions()
    {

    }
    virtual ~AbstractBoundaryConditions()
    {

    }

    virtual const std::vector<int>& getPeriodicBoundaryDimensions() const = 0;
    virtual const BCVec& getAdvectiveBoundaryConditions() const = 0;

    virtual bool contains(const std::string& elementName) const = 0;
    virtual const BCVec& getDiffusiveBoundaryConditions(const std::string& elementName) const = 0;
};



}

#endif // ABSTRACTBOUNDARYCONDITIONS_H
