#ifndef ABSTRACTCALCINTERFACEPROPERTIES_H
#define ABSTRACTCALCINTERFACEPROPERTIES_H

#include <vector>

namespace IPP
{

class AbstractCalcInterfaceProperties
{
public:
    virtual ~AbstractCalcInterfaceProperties() = default;

    virtual void setnComps(const size_t nComps) = 0;

    virtual void execute(const std::vector<double>& conc,
                         const std::vector<double>& porosity,
                         std::vector<double>& avgConc,
                         std::vector<double>& avgPoros) const = 0;

};

}

#endif // ABSTRACTCALCINTERFACEPROPERTIES_H
