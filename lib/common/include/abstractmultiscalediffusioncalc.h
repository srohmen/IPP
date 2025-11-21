#ifndef ABSTRACTMULTISCALEDIFFUSIONCALC_H
#define ABSTRACTMULTISCALEDIFFUSIONCALC_H

#include <vector>
#include <string>

namespace IPP
{

struct SimplePorosityInfos;

class AbstractMultiScaleDiffusionCalc
{
public:
    AbstractMultiScaleDiffusionCalc();

    virtual ~AbstractMultiScaleDiffusionCalc();


    virtual double calc(const SimplePorosityInfos& input) const = 0;
};

}


#endif // ABSTRACTMULTISCALEDIFFUSIONCALC_H
