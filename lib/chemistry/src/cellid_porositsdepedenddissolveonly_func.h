#ifndef CELLID_POROSITSDEPEDENDDISSOLVEONLY_FUNC_H
#define CELLID_POROSITSDEPEDENDDISSOLVEONLY_FUNC_H

#include "cellid_dissolveonly_func.h"

#include <vector>

namespace IPP
{

class CellID_PorositsDepedendDissolveOnly_Func : public CellID_DissolveOnly_Func
{
public:
    CellID_PorositsDepedendDissolveOnly_Func(const std::vector<double>& porosity,
                                             const std::vector<double>& threshold);

    virtual bool evaluate(const size_t iCell, const iterator& begin, const iterator& end) const;


private:
    const std::vector<double>& m_porosity;
    const std::vector<double>& m_threshold;
};

}

#endif // CELLID_POROSITSDEPEDENDDISSOLVEONLY_FUNC_H
