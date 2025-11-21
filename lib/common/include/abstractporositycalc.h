#ifndef ABSTRACTPOROSITYCALC_H
#define ABSTRACTPOROSITYCALC_H

#include <memory>
#include "arrayviewiterator.h"
#include "phasenametoinfos.h"

namespace IPP
{

struct SimplePorosityInfos;

class AbstractPorosityCalc
{
public:

    typedef ArrayViewIterator ValueIt;
    struct BeginEndIt
    {
        ValueIt begin;
        ValueIt end;
    };

    struct PorosityCalcInput
    {
        PorosityCalcInput()
            : inertFraction(0.0)
        {

        }

        double inertFraction;

        BeginEndIt volFractions;
        BeginEndIt phaseAmounts;
    };

    AbstractPorosityCalc() = default;
    virtual ~AbstractPorosityCalc() = default;

    virtual void setCompNames(const std::vector<std::string>* /*names*/)
    {

    }

    virtual void setPhaseNames(const std::vector<std::string>* /*names*/)
    {

    }

    virtual void setPhaseInfos(const PhaseNameToInfos* /*phaseInfos*/)
    {

    }


    using SimplePorosityInfosPtr = std::unique_ptr<SimplePorosityInfos>;
    virtual SimplePorosityInfosPtr calc(const PorosityCalcInput& input) const = 0;

};


}

#endif // ABSTRACTPOROSITYCALC_H
