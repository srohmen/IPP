#ifndef NONRETURNINGBCGEN_H
#define NONRETURNINGBCGEN_H

#include "bcgenerator.h"

namespace IPP
{

class NonReturningBCGen : public BCGenerator
{
public:
    NonReturningBCGen();

    virtual ConfigBoundaryConditions *generate() const override;
};


}

#endif // NONRETURNINGBCGEN_H
