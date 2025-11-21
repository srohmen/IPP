#ifndef PHASENAMETODISSPRECIPBEHAVIOUR_H
#define PHASENAMETODISSPRECIPBEHAVIOUR_H

#include "abstractdisspreciponlyinfo.h"
#include <map>

namespace IPP
{

class PhaseNameToDissPrecipBehaviour : public AbstractDissPrecipOnlyInfo
{
public:
    PhaseNameToDissPrecipBehaviour();

    virtual DissPrecipBehaviour getBehaviour(const std::string& phaseName) const override;

    void addData(const std::string& phaseName, const DissPrecipBehaviour& data);

private:
    std::map<std::string, DissPrecipBehaviour> m_dissPrecData;
};




}

#endif // PHASENAMETODISSPRECIPBEHAVIOUR_H
