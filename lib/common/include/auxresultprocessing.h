#ifndef AUXRESULTPROCESSING_H
#define AUXRESULTPROCESSING_H

#include <string>
#include <vector>

#include "array_view.h"

namespace IPP
{

class AuxResultProcessing
{
public:
    AuxResultProcessing() = default;
    virtual ~AuxResultProcessing() = default;

    virtual void setPorosity(const std::vector<double>& porosity)
    {

    }

    virtual void begin(const size_t iteration, const double& time,
                       const std::vector<std::string> &dataNames) = 0;

    virtual void process(const std::string& name,
                         const std::vector<double>& data) = 0;

    virtual void end() = 0;
};

}

#endif // AUXRESULTPROCESSING_H
