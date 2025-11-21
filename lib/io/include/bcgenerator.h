#ifndef BCGENERATOR_H
#define BCGENERATOR_H

#include "configfactorybase.h"
#include "configboundaryconditionsfwd.h"

namespace IPP
{

class BCGenerator : public ConfigFactoryBase
{
public:
    BCGenerator()
    {

    }

    virtual void setDim(const size_t dim)
    {
        // do nothing
        (void)dim;
    }

    virtual ConfigBoundaryConditions* generate() const = 0;

    virtual void setCompositionCookie(const void* cookie)
    {
        // do nothing
        (void)cookie;
    }

};

} // end of namespace

#endif // BCGENERATOR_H
