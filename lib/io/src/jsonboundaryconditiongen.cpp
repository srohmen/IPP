#include "jsonboundaryconditiongen.h"

#include <assert.h>
#include <iostream>
#include "configboundaryconditions.h"
#include "ippexception.h"
#include "configdatawrapper.h"


namespace IPP
{

namespace
{

typedef ConfigBoundaryConditions::AdvectiveBoundaryConditionData AdvectiveBoundaryConditionData;
typedef ConfigBoundaryConditions::DiffusiveBoundaryConditionData DiffusiveBoundaryConditionData;
typedef std::map<std::string, BoundaryConditionDomain> DomainNameToDomain;

static IPP::BoundaryConditionType convertToBCType(const std::string& bcTypeName)
{
    if(bcTypeName == "closed")
    {
        return IPP::BCT_Closed;
    }
    else if(bcTypeName == "constant")
    {
        return IPP::BCT_DensityDirichlet;
    }
    else if(bcTypeName == "flux")
    {
        return IPP::BCT_VelocityDirichlet;
    }
    else if(bcTypeName == "cauchy")
    {
        return IPP::BCT_DensityCauchy;
    }
    else if(bcTypeName == "open")
    {
        return IPP::BCT_VelocityOutflow;
    }
    else
    {
        throw std::runtime_error("Unknown Boundary Condition type: " + bcTypeName);
    }
}

struct ReadAdvectionValues
{
    void operator()(const Json::Value& condition, AdvectiveBoundaryConditionData& bc) const
    {
    // read density + velocity

    bc.density = condition.get("density", -1.0).asDouble();

    const Json::Value& velocityNode = condition["velocity"];
    bc.velocity = {{ velocityNode.get(0u, 0.0).asDouble(),
                     velocityNode.get(1u, 0.0).asDouble(),
                     velocityNode.get(2u, 0.0).asDouble()
                   }};
    }
};

typedef JSONBoundaryConditionGen::CompNameToComp CompNameToComp;
struct ReadDiffusiveValues
{
    ReadDiffusiveValues(const CompNameToComp& compNameToComp)
        : compNameToComp(compNameToComp)
    {

    }

    void operator()(const Json::Value& condition, DiffusiveBoundaryConditionData& bc) const
    {
        const std::string typeName = condition["type"].asString();
        if(typeName == "constant")
        {
            // read composition
            IPPCheck::assertCheck(condition.isMember("composition"), "Composition name not defined for BC");
            const std::string compositionName = condition["composition"].asString();

            const CompNameToComp::const_iterator it = compNameToComp.find(compositionName);
            if(it == compNameToComp.end())
            {
                throw std::runtime_error("composition not defined: " + compositionName);
            }
            else
            {
                bc.composition = it->second;
            }
        }
        else if(typeName == "open" || typeName == "closed" || typeName == "flux")
        {
            // do nothing
        }
        else
        {
            throw std::runtime_error("unsupported diffusive boundary type: " + typeName);
        }
    }

private:
    const CompNameToComp& compNameToComp;
};

template<typename ReadFunc, typename BCType>
static void read(const ReadFunc& func, const DomainNameToDomain &domainNameToBox,
                 const Json::Value& condition, BCType& bc)
{
    // read type
    const std::string typeName = condition["type"].asString();
    bc.type = convertToBCType(typeName);

    // read domain
    const std::string domainName = condition["domain"].asString();
    const DomainNameToDomain::const_iterator it = domainNameToBox.find(domainName);
    IPPCheck::assertCheck(it != domainNameToBox.end(),
                                "Domain not defined: " + domainName);
    const BoundaryConditionDomain& domain = it->second;
    bc.domain = domain;

    func(condition, bc);
}

template<typename ReadFunc, typename BCType>
static void readBC(const ReadFunc& func,
                   const DomainNameToDomain &domainNameToBox,
                   const Json::Value& conditions,
                   std::vector<BCType>& bcVec)
{
    for(size_t i = 0; i < conditions.size(); ++i)
    {
        const Json::Value condition = conditions[(int)i];
        // insert into vector and change property afterwards
        bcVec.push_back(BCType());
        BCType& bc = bcVec.back();

        read(func, domainNameToBox, condition, bc);
    }
}

static BoundaryConditionDomain::BoundaryPosition readPosition(const Json::Value& node)
{
    const std::string positionStr = node.asString();

    if("left" == positionStr)
    {
        return BoundaryConditionDomain::BP_left;
    }
    else if("right" == positionStr)
    {
        return BoundaryConditionDomain::BP_right;
    }
    else if("bottom" == positionStr)
    {
        return BoundaryConditionDomain::BP_bottom;
    }
    else if("top" == positionStr)
    {
        return BoundaryConditionDomain::BP_top;
    }
    else if("back" == positionStr)
    {
        return BoundaryConditionDomain::BP_back;
    }
    else if("front" == positionStr)
    {
        return BoundaryConditionDomain::BP_front;
    }
    else
    {
        throw std::runtime_error("Unknown Boundary Condition position: " + positionStr);
    }
}

static void readBoundaryConditions(const CompNameToComp& compNameToComp,
                                   const Json::Value& bcNode,
                                   const size_t dims,
                                   ConfigBoundaryConditions& bcOut)
{
    const Json::Value& periodic = bcNode["periodic"];

    for(int iDim = 0; iDim < (int)periodic.size(); ++iDim)
    {
        const Json::Value& node = periodic[iDim];
        if(node.asBool())
        {
            bcOut.periodicDims.push_back(iDim);
        }
    }




    DomainNameToDomain domainNameToData;

    const Json::Value& domainsNode = bcNode["domains"];
    for(size_t i = 0; i < domainsNode.size(); ++i)
    {
        const Json::Value domainNode = domainsNode[(int)i];

        IPPCheck::assertCheck(domainNode.isMember("name"), "name of BC domain not defined");
        const std::string name = domainNode["name"].asString();

        IPPCheck::assertCheck(domainNameToData.find(name) == domainNameToData.end(),
                                    "Boundary Condition domain defined multiple times: " + name);

        domainNameToData[name] = BoundaryConditionDomain();
        BoundaryConditionDomain& domain = domainNameToData[name];

        IPPCheck::assertCheck(domainNode.isMember("position"),
                              "position of BC domain \"" + name + "\" not defined");
        const Json::Value& positionNode = domainNode["position"];
        BoundaryConditionDomain::BoundaryPosition& pos = domain.position;
        pos = readPosition(positionNode);

        IPPCheck::assertCheck(domainNode.isMember("range"), "range of BC domain not defined");
        const Json::Value& rangeNode = domainNode["range"];
        const size_t nElements = (dims-1) * 2;
        IPPCheck::assertCheck(rangeNode.size() == nElements,
                              "Boundary condition domain \"" + name
                              + "\" has invalid range. Must have nElements: "
                              + std::to_string( nElements ) );

        // TODO: optimize that since it is static lookup...
        std::vector<size_t> indices;
        switch(pos)
        {
        case BoundaryConditionDomain::BP_left:
        case BoundaryConditionDomain::BP_right:
            indices.push_back(2);
            indices.push_back(3);
            indices.push_back(4);
            indices.push_back(5);
            break;

        case BoundaryConditionDomain::BP_bottom:
        case BoundaryConditionDomain::BP_top:
            indices.push_back(0);
            indices.push_back(1);
            indices.push_back(4);
            indices.push_back(5);
            break;

        case BoundaryConditionDomain::BP_back:
        case BoundaryConditionDomain::BP_front:
            indices.push_back(0);
            indices.push_back(1);
            indices.push_back(2);
            indices.push_back(3);
            break;

        default:
            throw std::runtime_error("Unknown Boundary Condition position: " + std::to_string(pos));
            break;
        }



        IPPBox3DInt& range = domain.range;
        for(size_t i = 0; i < rangeNode.size(); ++i)
        {
            const Json::Value& entry = rangeNode[(int)i];
            const int coord = entry.asInt();

            const size_t index = indices[i];
            range[index] = coord;
        }
    }

    ConfigBoundaryConditions::AdvectiveBCVec& advectiveBCs = bcOut.advectiveBC;
    const Json::Value& advBC = bcNode["advective"];
    readBC(ReadAdvectionValues(), domainNameToData, advBC, advectiveBCs);

    ConfigBoundaryConditions::DiffusiveBCVec& diffusiveBCs = bcOut.diffusiveBC;
    const Json::Value& diffBC = bcNode["diffusive"];
    readBC(ReadDiffusiveValues(compNameToComp), domainNameToData, diffBC, diffusiveBCs);

}

}

JSONBoundaryConditionGen::JSONBoundaryConditionGen()
    : m_dim(-1),
      m_compNameToComp(nullptr)
{

}

void JSONBoundaryConditionGen::setDim(const size_t dim)
{
    assert(dim <= 3);
    m_dim = dim;
}

IPP::ConfigBoundaryConditions* JSONBoundaryConditionGen::generate() const
{
    IPPCheck::assertCheck(m_dim <= 3);

    IPPCheck::assertCheck(m_compNameToComp);
    IPPCheck::assertCheck(m_configData);

    ConfigBoundaryConditions* out = new ConfigBoundaryConditions;

    readBoundaryConditions(*m_compNameToComp, m_configData->node, m_dim, *out);

    return out;
}

void JSONBoundaryConditionGen::setCompositionCookie(const void* cookie)
{
    m_compNameToComp = static_cast<const CompNameToComp*>(cookie);
    IPPCheck::assertCheck(m_compNameToComp);
}


} // end of namespace LBGeoChem
