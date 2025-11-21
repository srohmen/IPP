#ifndef POROSITYMODELFACTORY_H
#define POROSITYMODELFACTORY_H

#include "configfactorybase.h"

namespace IPP
{

class AbstractPorosityCalc;

class PorosityModelFactory : public ConfigFactoryBase
{
public:
    PorosityModelFactory() = default;
    virtual ~PorosityModelFactory() = default;

    virtual AbstractPorosityCalc* generate() const = 0;
};

}

#endif // POROSITYMODELFACTORY_H
