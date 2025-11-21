#ifndef PHREEQCGLOBALINITDATA_H
#define PHREEQCGLOBALINITDATA_H

#include "phreeqcfeatures.h"

namespace IPP
{

struct PhreeqcGlobalInitData
{
    PhreeqcGlobalInitData()
        :  genericCellID(-1)
    {

    }

    size_t genericCellID;
    CellIdToPhreeqcFeatureMask featureMask;
    std::vector<double> tmpGlobalPoros;
    std::vector<double> tmpGlobalPorosCapillary;
    std::vector<double> tmpGlobalDiff;
    std::vector<double> tmpGlobalInertVolFrac;
};


}

#endif // PHREEQCGLOBALINITDATA_H
