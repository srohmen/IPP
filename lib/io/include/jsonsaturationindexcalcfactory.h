#ifndef JSONSATURATIONINDEXCALCFACTORY_H
#define JSONSATURATIONINDEXCALCFACTORY_H

#include "saturationindexcalcfactory.h"
#include "configfactorybase.h"

namespace IPP
{

class JSONSaturationIndexCalcFactory : public SaturationIndexCalcFactory
{
public:
    JSONSaturationIndexCalcFactory();

    virtual void setSpatialResolution(const double& resolution) override;
    virtual void setDissPrecipitationBehaviour(const AbstractDissPrecipOnlyInfo* dissPrec) override;

    virtual AbstractSICalc* generate() const override;

private:
    double m_spatialResolution;
    const AbstractDissPrecipOnlyInfo* m_dissPrec;
};


}

#endif // JSONSATURATIONINDEXCALCFACTORY_H
