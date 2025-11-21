#include "jsonporositymodelfactory.h"

#include "configdatawrapper.h"

#include "simpleporositycalc.h"
#include "tennisjenningscshporositycalc.h"

namespace IPP
{

static AbstractPorosityCalc* parseSimple(const Json::Value& node)
{
    return new SimplePorosityCalc;
}

static AbstractPorosityCalc* parseTenJen(const Json::Value& node)
{
    const double wc = node.get("wc", -1.0).asDouble();
    const double degreeHydration = node.get("degree_hydration", -1.0).asDouble();

    if(wc <= 0.0 || wc >= 1.0 || degreeHydration < 0.0 || degreeHydration > 1.0)
    {
        throw std::runtime_error("invalid cement hydration parameter passed: wc = "
                                 + std::to_string(wc) + " ; degree_hydration = "
                                 + std::to_string(degreeHydration));
    }

    TennisJenningsCSHPorosityCalc* calc = new TennisJenningsCSHPorosityCalc(wc, degreeHydration);

    const Json::Value& phasesNode = node["phases"];
    for(size_t i = 0; i < phasesNode.size(); ++i)
    {
        const Json::Value& entry = phasesNode[(int)i];
        const std::string phaseName = entry.asString();
        calc->add_CSH_Phase(phaseName);
    }


    return calc;
}

JSONPorosityModelFactory::JSONPorosityModelFactory()
{

}

AbstractPorosityCalc* JSONPorosityModelFactory::generate() const
{
    AbstractPorosityCalc* porosCalc = nullptr;

    const Json::Value& node = m_configData->node;

    const std::string typeStr = node.get("type", "").asString();

    if(typeStr == "simple")
    {
        porosCalc = parseSimple(node);
    }
    else if(typeStr == "TenJen")
    {
        porosCalc = parseTenJen(node);
    }
    else
    {
        throw std::runtime_error("Unknown porosity model type: \"" + typeStr + "\"");
    }


    return porosCalc;
}

}
