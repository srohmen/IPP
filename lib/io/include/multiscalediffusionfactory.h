#ifndef MULTISCALEDIFFUSIONFACTORY_H
#define MULTISCALEDIFFUSIONFACTORY_H

#include "configfactorybase.h"

namespace IPP
{

class AbstractMultiScaleDiffusionCalc;

class MultiScaleDiffusionFactory : public ConfigFactoryBase
{
public:
    MultiScaleDiffusionFactory();

    virtual ~MultiScaleDiffusionFactory();

    virtual AbstractMultiScaleDiffusionCalc* generate() const = 0;
};


}

#endif // MULTISCALEDIFFUSIONFACTORY_H
