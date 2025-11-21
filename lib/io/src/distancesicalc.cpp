#include <iostream>
#include <cmath>
#include <assert.h>
#include <limits>

#include "distancesicalc.h"
#include "abstractscalarfield.h"
#include "abstractdisspreciponlyinfo.h"
#include "ippexception.h"
#include "phasenametoinfos.h"


namespace IPP
{

static double applySIData(const DistanceSICalc::SIData& data,
                          const double& physicalDistance)
{
    double newSI = 0.0;

    if(physicalDistance <= 0.0)
    {
        // within solid grain
        newSI = 0.0;
    }
    else if(physicalDistance <= data.distanceThresh)
    {
        // within grain neigborhood

        switch(data.technique)
        {
            case(DistanceSICalc::SICT_Threshold):
            {
                newSI = 0.0;
                break;
            }
            case(DistanceSICalc::SICT_Linear):
            {
                const double slope = data.maxSI / data.distanceThresh;
                newSI = physicalDistance * slope;
                break;
            }
            case(DistanceSICalc::SICT_Exponential):
            {
                const double k = -std::log(0.01) / data.distanceThresh;
                newSI = (1.0 - std::exp(-k*physicalDistance)) * data.maxSI;
                break;
            }
            default:
            {
                throw std::runtime_error("Unknown SI calc technique");
                break;
            }
        }

    }
    else
    {
        // far away from grain
        newSI = data.maxSI;
    }

    return newSI;
}

static const double s_highSI = 100.0;
static const double s_ultraLowPoros = 1.0E-3;
static const double s_veryLowPoros = 1.0E-2;
static const double s_lowPoros = 1.0E-1;

class PorosToSIcalc
{
public:
    PorosToSIcalc(const double& x, const double& c0, const double& Vc)
        : x(x),
          c0(c0),
          Vc(Vc)
    {

    }

    double calcPorosDiff(const double dSI) const
    {
        const double deLog = std::pow(10.0, dSI/x);
        const double foo = deLog - 1;
        const double result = foo * c0 * Vc;
        return result;
    }

    double calcSIDiff(const double dP) const
    {
        const double c0Vc = Vc*c0;
        const double clamped_dP = std::max(dP, -(c0Vc - 1.0E-6));
        const double result = x * std::log10(1 + clamped_dP / c0Vc);
        return result;
    }

private:
    const double x;
    const double c0;
    const double Vc;
};


DistanceSICalc::DistanceSICalc(const double& spatialResolution)
    : m_spatialResolution(spatialResolution),
      m_distField(),
      m_porosField(),
      m_phaseInfos(nullptr)
{

}

DistanceSICalc::~DistanceSICalc()
{

}

void DistanceSICalc::setDistanceField(const ConstAbstractScalarFieldPtr& distField)
{
    m_distField = distField;
}

void DistanceSICalc::setPorosityField(const ConstAbstractScalarFieldPtr& porosField)
{
    m_porosField = porosField;
}

double DistanceSICalc::evaluate(const IPPVector3DInt& coord, const std::string& phaseName) const
{
    // check for dissolve only phase
    const DissPrecipBehaviour& behaviour = m_dissPrec->getBehaviour(phaseName);
    if(behaviour == DPB_DissOnly)
    {
        return 1000.0;
    }
    else if(behaviour == DPB_PrecOnly)
    {
        return -1000.0;
    }
    else if(behaviour == DPB_Normal)
    {
        const double latticeDistance = m_distField->get(coord[0], coord[1], coord[2]);
        // std::cout << coord[0] << "\t" << coord[1] << "\t" << coord[2] << "\t" << distance << std::endl;

        const double physicalDistance = latticeDistance * m_spatialResolution - m_spatialResolution;

        double newSI;

        PhaseNameToData::const_iterator it = m_phaseNameToData.find(phaseName);
        if(it == m_phaseNameToData.end())
        {
            // try to find generic phase SI calc definition
            it = m_phaseNameToData.find("generic");
        }

        if(it != m_phaseNameToData.end())
        {
            const SIData& data = it->second;
            newSI = applySIData(data, physicalDistance);
        }
        else
        {
            // not defined SI func.
            // assume precipitation on grain surface only
            if(physicalDistance <= m_spatialResolution)
            {
                newSI = 0.0;
            }
            else
            {
                // incredible high SI value....
                newSI = 100.0;
            }
        }

        const double poros = m_porosField->get(coord[0], coord[1], coord[2]);


        if(poros < 0.0)
        {
            newSI = poros * -10.0;
            newSI = std::min(newSI, 1.0);
        }

        assert(newSI >= 0.0);
        return newSI;

    }
    else
    {
        throw std::runtime_error("Unknown dissolution/precipitation behaviour");
    }

}

void DistanceSICalc::evaluate(const IPPVector3DInt& coord, const iterator& begin, const iterator& end) const
{
    IPPCheck::assertCheck(end > begin);
    IPPCheck::assertCheck((size_t)(end - begin) == m_phaseNameToSaturationIndices.size());

    iterator it = begin;

    for(auto phaseIt = m_phaseNameToSaturationIndices.begin();
        phaseIt != m_phaseNameToSaturationIndices.end(); ++phaseIt, ++it)
    {
        const std::string& phaseName = phaseIt->first;

        const double newSI = evaluate(coord, phaseName);

        *it = newSI;
    }
}



void DistanceSICalc::setSIFunc(const std::string& phaseName, const SIData& data)
{
    m_phaseNameToData[phaseName] = data;
}

void DistanceSICalc::setDissPrecipBehaviour(const AbstractDissPrecipOnlyInfo *dissPrec)
{
    m_dissPrec = dissPrec;
}

void DistanceSICalc::setPhaseInfos(const PhaseNameToInfos* phaseInfos)
{
    m_phaseInfos = phaseInfos;
}

void DistanceSICalc::setConcentrations(const size_t iComp, const AbstractScalarFieldPtr& concField)
{
    if(iComp >= m_concVec.size())
    {
        m_concVec.resize(iComp+1);
    }
    m_concVec[iComp] = concField;
}

void DistanceSICalc::setSaturationIndices(const std::string& phaseName, const AbstractScalarFieldPtr& saturationIndices)
{
    m_phaseNameToSaturationIndices[phaseName] = saturationIndices;
}

void DistanceSICalc::setPrecipAmount(const std::string& phaseName, const AbstractScalarFieldPtr& precipAmount)
{
    m_phaseNameToPrecipAmount[phaseName] = precipAmount;
}

} // end of namespace LBGeoChem
