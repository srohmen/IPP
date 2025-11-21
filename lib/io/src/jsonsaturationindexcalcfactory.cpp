#include "jsonsaturationindexcalcfactory.h"

#include "ippexception.h"
#include "distancesicalc.h"
#include "configdatawrapper.h"
#include "pcs_si_func.h"

namespace IPP
{

JSONSaturationIndexCalcFactory::JSONSaturationIndexCalcFactory()
    : m_spatialResolution()
{

}

void JSONSaturationIndexCalcFactory::setSpatialResolution(const double& resolution)
{
    m_spatialResolution = resolution;
}

void JSONSaturationIndexCalcFactory::setDissPrecipitationBehaviour(const AbstractDissPrecipOnlyInfo *dissPrec)
{
    m_dissPrec = dissPrec;
}

AbstractSICalc* JSONSaturationIndexCalcFactory::generate() const
{
    const Json::Value& node = m_configData->node;
    const std::string technique = node.get("technique", "").asString();
    if(technique == "pcs")
    {
        const std::string poreGeom = node["pore_geometry"].asString();
        if(poreGeom == "cylinder")
        {
            const double temp = node.get("T", 298).asDouble();
            const double porosInit = node.get("poros_init", -1.0).asDouble();
            const double rInit = node.get("r_init", -1.0).asDouble();
            // mol/m3
            const double surfTensFallbackConc = node.get("surface_tension_fallback_conc", 1.0e-3).asDouble() * 1000;

            std::map<std::string, double> surfaceTensions;
            const Json::Value& surfTensionsNode = node["surface_tensions"];
            for ( size_t iEntry = 0; iEntry < surfTensionsNode.size(); ++iEntry )
            {
                const Json::Value& value = surfTensionsNode[(int)iEntry];
                const Json::Value::Members keys = value.getMemberNames();
                const std::string phaseName = keys.front();
                const double surfTens = value[phaseName].asDouble();

                // std::cout << "phaseName: " << phaseName << std::endl;

                IPPCheck::assertCheck(surfaceTensions.find(phaseName) == surfaceTensions.end(), "double entry of: " + phaseName);
                surfaceTensions[phaseName] = surfTens;
            }

            return new PCS_FuncCylinder(temp, porosInit, rInit, surfTensFallbackConc, surfaceTensions);
        }
        else
        {
            throw std::runtime_error("unknown pore geometry: " + poreGeom);
        }

    }
    else
    {
        throw std::runtime_error("unknown PCS implementation: " + technique);
    }

}

}
