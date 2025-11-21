#ifndef NONRETURNINGSICALCFACTORY_H
#define NONRETURNINGSICALCFACTORY_H

#include "saturationindexcalcfactory.h"

namespace IPP
{

class NonReturningSICalcFactory : public SaturationIndexCalcFactory
{
public:
    NonReturningSICalcFactory();

    virtual AbstractSICalc* generate() const override;
};

}

#endif // NONRETURNINGSICALCFACTORY_H
