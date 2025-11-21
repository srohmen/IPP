#ifndef NONRETURNINGPOROSITYMODELFACTORY_H
#define NONRETURNINGPOROSITYMODELFACTORY_H

#include "porositymodelfactory.h"

namespace IPP
{

class NonReturningPorosityModelFactory : public PorosityModelFactory
{
public:
    NonReturningPorosityModelFactory();

    virtual AbstractPorosityCalc* generate() const override;
};


}

#endif // NONRETURNINGPOROSITYMODELFACTORY_H
