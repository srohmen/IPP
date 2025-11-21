#ifndef PHREEQCINITIALDATA_H
#define PHREEQCINITIALDATA_H

#include <vector>

namespace IPP
{

struct PhreeqcInitialData
{
    std::vector<double> conc;
    std::vector<double> porosity;
};


}

#endif // PHREEQCINITIALDATA_H
