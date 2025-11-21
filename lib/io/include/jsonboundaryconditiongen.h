#ifndef JSONBOUNDARYCONDITIONGEN_H
#define JSONBOUNDARYCONDITIONGEN_H

#include "bcgenerator.h"
#include "ippconfig.h"

namespace IPP
{

class Composition;

class JSONBoundaryConditionGen : public BCGenerator
{
public:
    JSONBoundaryConditionGen();

    virtual void setDim(const size_t dim);

    virtual ConfigBoundaryConditions* generate() const;

    virtual void setCompositionCookie(const void* cookie);

    typedef std::map<std::string, Composition*> CompNameToComp;

private:
    size_t m_dim;
    const CompNameToComp* m_compNameToComp;

};

}
#endif // JSONBOUNDARYCONDITIONGEN_H
