#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "abstractboundaryconditions.h"

namespace IPP
{

class BoundaryConditions : public AbstractBoundaryConditions
{
public:
    BoundaryConditions();

    virtual const std::vector<int>& getPeriodicBoundaryDimensions() const override;
    virtual const BCVec& getAdvectiveBoundaryConditions() const override;

    virtual bool contains(const std::string& elementName) const override;
    virtual const BCVec& getDiffusiveBoundaryConditions(const std::string& elementName) const override;

    std::vector<int> periodicDims;
    BCVec advectiveBC;

    typedef std::map<std::string, BCVec> ElementNameToBCVec;
    ElementNameToBCVec diffusiveBC;
};

typedef std::shared_ptr<BoundaryConditions> BoundaryConditionsPtr;

}

#endif // BOUNDARYCONDITIONS_H
