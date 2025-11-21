#ifndef SATURATIONINDEXCALCFACTORY_H
#define SATURATIONINDEXCALCFACTORY_H

#include "configfactorybase.h"

namespace IPP
{

class AbstractSICalc;
class AbstractDissPrecipOnlyInfo;

class SaturationIndexCalcFactory : public ConfigFactoryBase
{
public:
    SaturationIndexCalcFactory()
    {

    }

    virtual ~SaturationIndexCalcFactory()
    {

    }

    virtual void setSpatialResolution(const double&)
    {
    }

    virtual void setDissPrecipitationBehaviour(const AbstractDissPrecipOnlyInfo*)
    {
    }

    virtual AbstractSICalc* generate() const = 0;
};

} // end of namespace

#endif // SATURATIONINDEXCALCFACTORY_H
