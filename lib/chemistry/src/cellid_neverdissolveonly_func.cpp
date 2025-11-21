#include "cellid_neverdissolveonly_func.h"

namespace IPP
{

CellID_NeverDissolveOnly_Func::CellID_NeverDissolveOnly_Func(const size_t nCells)
    : m_nCells(nCells)
{

}

size_t CellID_NeverDissolveOnly_Func::nCells() const
{
    return m_nCells;
}

bool CellID_NeverDissolveOnly_Func::evaluate(const size_t /*iCell*/,
                                             const iterator& begin,
                                             const iterator& end) const
{
    for(iterator it = begin; it != end; ++it)
    {
        *it = false;
    }

    return false;
}


}
