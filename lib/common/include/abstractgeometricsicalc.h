#ifndef ABSTRACTGEOMETRICSICALC_H
#define ABSTRACTGEOMETRICSICALC_H

#include <vector>
#include <string>

#include "ippvector.h"
#include "abstractscalarfieldfwd.h"

// TODO: find abstraction
#include "phasenametoinfos.h"

namespace IPP
{

class AbstractGeometricSICalc
{
public:
    typedef double return_type;

    AbstractGeometricSICalc()
    {

    }

    virtual ~AbstractGeometricSICalc()
    {

    }

    virtual void setDistanceField(const ConstAbstractScalarFieldPtr& distField) = 0;
    virtual void setPorosityField(const ConstAbstractScalarFieldPtr& porosField) = 0;
    virtual void setPhaseInfos(const PhaseNameToInfos* phaseInfos) = 0;
    virtual void setConcentrations(const size_t iComp, const AbstractScalarFieldPtr& concField) = 0;
    virtual void setSaturationIndices(const std::string& phaseName, const AbstractScalarFieldPtr& saturationIndices) = 0;
    virtual void setPrecipAmount(const std::string& phaseName, const AbstractScalarFieldPtr& precipAmount) = 0;

    typedef std::vector<return_type>::iterator iterator;
    virtual void evaluate(const IPPVector3DInt& coord, const iterator& begin, const iterator& end) const = 0;
};

} // end of namespace

#endif // ABSTRACTGEOMETRICSICALC_H
