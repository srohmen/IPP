#ifndef ABSTRACTDISTRIBUTEEXCESSTOTALS_H
#define ABSTRACTDISTRIBUTEEXCESSTOTALS_H

#include "celltotalsdiff.h"

namespace IPP
{

class AbstractDistributeExcessTotals
{
public:
    virtual ~AbstractDistributeExcessTotals() = default;

    virtual void init() = 0;

    virtual void setnComps(const std::size_t nComps) = 0;
    virtual void updatePorosity(const std::vector<double>& porosity) = 0;

    virtual void execute(const std::vector<CellTotalsDiff>& diffPerCell,
                         std::vector<double>& distributedSource) = 0;
};

}

#endif // ABSTRACTDISTRIBUTEEXCESSTOTALS_H
