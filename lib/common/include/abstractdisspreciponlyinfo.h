#ifndef ABSTRACTDISSPRECIPONLYINFO_H
#define ABSTRACTDISSPRECIPONLYINFO_H

#include <string>
#include "dissprecipbehaviour.h"

namespace IPP
{

class AbstractDissPrecipOnlyInfo
{
public:

    AbstractDissPrecipOnlyInfo()
    {

    }

    virtual ~AbstractDissPrecipOnlyInfo()
    {

    }

    virtual DissPrecipBehaviour getBehaviour(const std::string& phaseName) const = 0;
};

}

#endif // ABSTRACTDISSPRECIPONLYINFO_H
