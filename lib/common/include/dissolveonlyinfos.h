#ifndef DISSOLVEONLYINFOS_H
#define DISSOLVEONLYINFOS_H

#include <vector>
#include "cellneighborinfo.h"

namespace IPP
{

struct DissolveOnlyInfos
{
    const double& porosity;
    const CellNeighborInfo& neighInfo;

    using ArrIterator = std::vector<double>::const_iterator;
    const ArrIterator& SI_begin;

    const std::vector<size_t>& phaseToNuclPhaseID;
    const std::vector<size_t>& nuclPhaseToMonomerID;
    const ArrIterator& monomerConc_begin;
};

}

#endif // DISSOLVEONLYINFOS_H
