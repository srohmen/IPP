#include "simpleporositycalc.h"

#include <cassert>

#include "simpleporosityinfos.h"

namespace IPP
{

SimplePorosityCalc::SimplePorosityCalc()
{

}

AbstractPorosityCalc::SimplePorosityInfosPtr SimplePorosityCalc::calc(const PorosityCalcInput &input) const
{

    const BeginEndIt& range = input.volFractions;

    double solidFrac = 0.0;
    for(ValueIt it = range.begin; it != range.end; ++it)
    {
        const double val = *it;
        solidFrac += val;
    }

    solidFrac += input.inertFraction;

    const double porosity = 1.0 - solidFrac;

    std::unique_ptr<SimplePorosityInfos> output(new SimplePorosityInfos);

    output->porosityTotal = porosity;
    output->porosityCapillary = porosity;

    return output;
}



}
