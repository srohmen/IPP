#ifndef AUXDATA_H
#define AUXDATA_H

#include <vector>
#include <string>

namespace IPP
{

typedef std::vector<double> AuxData;

struct AuxDataName
{
    std::string name;
    AuxData data;
};

typedef std::vector<AuxDataName> AuxDataVec;

}

#endif // AUXDATA_H
