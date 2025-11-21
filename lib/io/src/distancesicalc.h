#ifndef DISTANCESICALC_H
#define DISTANCESICALC_H

#include <map>
#include <vector>
#include "abstractgeometricsicalc.h"


namespace IPP
{
class AbstractDissPrecipOnlyInfo;

class DistanceSICalc : public AbstractGeometricSICalc
{
public:
    enum SI_Calc_DistanceTechnique
    {
        SICT_Threshold,
        SICT_Linear,
        SICT_Exponential
    };

    struct SIData
    {
        SI_Calc_DistanceTechnique technique;
        double distanceThresh;
        double maxSI;
        double molarVolume;
    };

    DistanceSICalc(const double& spatialResolution);
    virtual ~DistanceSICalc();

    virtual void setDistanceField(const ConstAbstractScalarFieldPtr& distField);
    virtual void setPorosityField(const ConstAbstractScalarFieldPtr& porosField);
    virtual void setPhaseInfos(const PhaseNameToInfos* phaseInfos);
    virtual void setConcentrations(const size_t iComp, const AbstractScalarFieldPtr& concField);
    virtual void setSaturationIndices(const std::string& phaseName, const AbstractScalarFieldPtr& saturationIndices);
    virtual void setPrecipAmount(const std::string& phaseName, const AbstractScalarFieldPtr& precipAmount);

    virtual void evaluate(const IPPVector3DInt& coord, const iterator& begin, const iterator& end) const;

    void setSIFunc(const std::string& phaseName, const SIData& data);

    void setDissPrecipBehaviour(const AbstractDissPrecipOnlyInfo *dissPrec);


private:
    double evaluate(const IPPVector3DInt& coord, const std::string& phaseName) const;

    const double m_spatialResolution;
    typedef std::map<std::string, SIData> PhaseNameToData;
    PhaseNameToData m_phaseNameToData;
    const AbstractDissPrecipOnlyInfo* m_dissPrec;

    ConstAbstractScalarFieldPtr m_distField;
    ConstAbstractScalarFieldPtr m_porosField;
    const PhaseNameToInfos* m_phaseInfos;
    std::vector<AbstractScalarFieldPtr> m_concVec;
    std::map<std::string, AbstractScalarFieldPtr> m_phaseNameToSaturationIndices;
    std::map<std::string, AbstractScalarFieldPtr> m_phaseNameToPrecipAmount;

};

} // end of namespace LBGeoChem

#endif // DISTANCESICALC_H
