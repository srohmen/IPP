#ifndef CELLID_NEVERDISSOLVEONLY_FUNC_H
#define CELLID_NEVERDISSOLVEONLY_FUNC_H

#include "cellid_dissolveonly_func.h"

namespace IPP
{

class CellID_NeverDissolveOnly_Func : public CellID_DissolveOnly_Func
{
public:
    CellID_NeverDissolveOnly_Func(const size_t nCells);

    virtual size_t nCells() const;
    virtual bool evaluate(const size_t iCell, const iterator& begin, const iterator& end) const;

private:
    const size_t m_nCells;

};

}

#endif // CELLID_NEVERDISSOLVEONLY_FUNC_H
