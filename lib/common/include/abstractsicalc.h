#ifndef ABSTRACTSICALC_H
#define ABSTRACTSICALC_H

#include <vector>
#include <string>

#include "phasenametoinfos.h"

namespace IPP
{

class AbstractSICalc
{
public:
    typedef double return_type;

    AbstractSICalc()
    {

    }

    virtual ~AbstractSICalc()
    {

    }

    virtual void init(const std::vector<std::string> &phaseNames,
                      const PhaseNameToInfos& phaseInfos,
                      const std::vector<double>* poros) = 0;

    typedef std::vector<return_type>::iterator iterator;
    virtual void evaluate(const size_t cellId, const iterator& begin, const iterator& end) const = 0;
};

} // end of namespace

#endif // ABSTRACTGEOMETRICSICALC_H
