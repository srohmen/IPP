#ifndef NONRETURNINGMULTISCALEDIFFUSIONFACTORY_H
#define NONRETURNINGMULTISCALEDIFFUSIONFACTORY_H

#include "multiscalediffusionfactory.h"

namespace IPP
{

class NonReturningMultiScaleDiffusionFactory : public MultiScaleDiffusionFactory
{
    virtual AbstractMultiScaleDiffusionCalc* generate() const;
};

}


#endif // NONRETURNINGMULTISCALEDIFFUSIONFACTORY_H
