#ifndef SIMPLEPOROSITYCALC_H
#define SIMPLEPOROSITYCALC_H

#include "abstractporositycalc.h"

namespace IPP
{

class SimplePorosityCalc : public AbstractPorosityCalc
{
public:
    SimplePorosityCalc();


    virtual SimplePorosityInfosPtr calc(const PorosityCalcInput& input) const override;

private:

};

}

#endif // SIMPLEPOROSITYCALC_H
