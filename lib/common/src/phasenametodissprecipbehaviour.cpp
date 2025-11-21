#include "phasenametodissprecipbehaviour.h"

namespace IPP
{

PhaseNameToDissPrecipBehaviour::PhaseNameToDissPrecipBehaviour()
{

}

DissPrecipBehaviour PhaseNameToDissPrecipBehaviour::getBehaviour(const std::string &phaseName) const
{
    const auto it = m_dissPrecData.find(phaseName);
    if(it == m_dissPrecData.end())
    {
        return DPB_Normal;
    }
    else
    {
        return it->second;
    }
}

void PhaseNameToDissPrecipBehaviour::addData(const std::string &phaseName, const DissPrecipBehaviour &data)
{
    m_dissPrecData[phaseName] = data;
}

}
