#ifndef JSONPOROSITYMODELFACTORY_H
#define JSONPOROSITYMODELFACTORY_H

#include "porositymodelfactory.h"

namespace IPP
{

class JSONPorosityModelFactory : public PorosityModelFactory
{
public:
    JSONPorosityModelFactory();

    virtual AbstractPorosityCalc* generate() const override;
};

}

#endif // JSONPOROSITYMODELFACTORY_H
