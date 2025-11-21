#include "cellid_porositsdepedenddissolveonly_func.h"

#include "ippexception.h"

namespace IPP
{


CellID_PorositsDepedendDissolveOnly_Func::CellID_PorositsDepedendDissolveOnly_Func(const std::vector<double>& porosity,
                                                                                   const std::vector<double>& threshold)
    : m_porosity(porosity)
    , m_threshold(threshold)
{
    IPPCheck::assertCheck(porosity.size() == threshold.size());
}

bool CellID_PorositsDepedendDissolveOnly_Func::evaluate(const size_t iCell,
                                                        const iterator& begin,
                                                        const iterator& end) const
{
    const double& porosity = m_porosity[iCell];
    const double& threshold = m_threshold[iCell];
    const bool result = porosity <= threshold;

    for(iterator it = begin; it != end; ++it)
    {
        *it = result;
    }

    return result;
}



}
