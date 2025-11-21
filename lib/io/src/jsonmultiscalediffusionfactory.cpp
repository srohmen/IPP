#include "jsonmultiscalediffusionfactory.h"

#include "configdatawrapper.h"
#include "ippexception.h"

#include "degradedcshdiffusioncalc.h"
#include "weightedaveragediffusioncalc.h"
#include "moritanakadiffusioncalc.h"
#include "cshgeldiffusioncoefcalc.h"
#include "archieslawdiffusioncalc.h"
#include "lowporosdegradedcshdiffusioncalc.h"
#include "lowporosthreshdegradedcshdiffusioncalc.hpp"

namespace IPP
{

static double parseFreePoreDiffCoef(const Json::Value& node)
{
    double value = 1.0E-9; // m2/s

    if(node.isMember("free_pore_diff_coef"))
    {
        value = node["free_pore_diff_coef"].asDouble();
    }
    else
    {
        std::cout << "WARNING: \"free_pore_diff_coef\" not defined in diff calc, assuming: "
                  << value << " m2/s" << std::endl;
    }

    return value;
}

static AbstractMultiScaleDiffusionCalc* parseDegradedCSHInterpolation(const Json::Value& node)
{
    const double D0 = parseFreePoreDiffCoef(node);

    DegradedCSHDiffusionCalc* multiScaleCalc =
            new DegradedCSHDiffusionCalc(new CSHGelDiffusionCoefCalc,
                                         new MoriTanakaDiffusionCalc,
                                         D0);

    return multiScaleCalc;
}



static double parseArchiesExpo(const Json::Value& node)
{
    double value = 2.0;

    if(node.isMember("expo"))
    {
        value = node["expo"].asDouble();
    }
    else
    {
        std::cout << "WARNING: \"expo\" not defined in diff calc, assuming: " << value << std::endl;
    }

    return value;
}


static AbstractMultiScaleDiffusionCalc* parseLowPorosDegradedCSHInterpolation(const Json::Value& node)
{
    const double D0 = parseFreePoreDiffCoef(node);
    const double expo = parseArchiesExpo(node);

    double porosThreshArchies = 0.0;
    if(node.isMember("poros_thresh_archies"))
    {
        porosThreshArchies = node["poros_thresh_archies"].asDouble();
    }
    else
    {
        std::cout << "WARNING: \"poros_thresh_archies\" not defined in diff calc, assuming: "
                  << porosThreshArchies << std::endl;
    }


    double porosThreshClog = 0.0;
    if(node.isMember("poros_thresh_clog"))
    {
        porosThreshClog = node["poros_thresh_clog"].asDouble();
    }
    else
    {
        std::cout << "WARNING: \"poros_thresh_clog\" not defined in diff calc, assuming: "
                  << porosThreshClog << std::endl;
    }

    LowPorosDegradedCSHDiffusionCalc* multiScaleCalc =
            new LowPorosDegradedCSHDiffusionCalc(new CSHGelDiffusionCoefCalc,
                                                 new MoriTanakaDiffusionCalc,
                                                 D0, expo, porosThreshArchies,
                                                 porosThreshClog);
    return multiScaleCalc;
}


static AbstractMultiScaleDiffusionCalc* parseLowPorosThreshDegradedCSHInterpolation(const Json::Value& node)
{
    const double D0 = parseFreePoreDiffCoef(node);

    double porosThreshClog = 0.0;
    if(node.isMember("poros_thresh_clog"))
    {
        porosThreshClog = node["poros_thresh_clog"].asDouble();
    }
    else
    {
        std::cout << "WARNING: \"poros_thresh_clog\" not defined in diff calc, assuming: "
                  << porosThreshClog << std::endl;
    }

    LowPorosThreshDegradedCSHDiffusionCalc* multiScaleCalc =
            new LowPorosThreshDegradedCSHDiffusionCalc(new CSHGelDiffusionCoefCalc,
                                                       new MoriTanakaDiffusionCalc,
                                                       D0, porosThreshClog);
    return multiScaleCalc;
}


static AbstractMultiScaleDiffusionCalc* parseArchiesLaw(const Json::Value& node)
{
    const double D0 = parseFreePoreDiffCoef(node);
    const double expo = parseArchiesExpo(node);

    return new ArchiesLawDiffusionCalc(D0, expo);
}


JSONMultiScaleDiffusionFactory::JSONMultiScaleDiffusionFactory()
{

}

JSONMultiScaleDiffusionFactory::~JSONMultiScaleDiffusionFactory()
{

}

AbstractMultiScaleDiffusionCalc *JSONMultiScaleDiffusionFactory::generate() const
{
    AbstractMultiScaleDiffusionCalc* multiScaleCalc = nullptr;

    const Json::Value& node = m_configData->node;

    const std::string typeStr = node.get("type", "").asString();

    if(typeStr == "degraded_CSH")
    {
        multiScaleCalc = parseDegradedCSHInterpolation(node);
    }
    else if(typeStr == "archies")
    {
        multiScaleCalc = parseArchiesLaw(node);
    }
    else if(typeStr == "low_poros_degraded_CSH")
    {
        multiScaleCalc = parseLowPorosDegradedCSHInterpolation(node);
    }
    else if(typeStr == "low_poros_thresh_degraded_CSH")
    {
        multiScaleCalc = parseLowPorosThreshDegradedCSHInterpolation(node);
    }
    else
    {
        throw std::runtime_error("Unknown multi scale diffusion type: \"" + typeStr + "\"");
    }


    return multiScaleCalc;
}


}
