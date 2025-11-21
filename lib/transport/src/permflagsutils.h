#ifndef PERMFLAGSUTILS_H
#define PERMFLAGSUTILS_H

#include "permeabiltyflag.h"
#include "setpermmaskflag.h"

namespace IPP
{

namespace PermFlagUtils
{

template<typename MaskType, typename DotList, typename Field>
void updatePermMask(const DotList& permNodes, const DotList& nonPermNodes, Field& permMask)
{
    static constexpr size_t dim = PlbTypeDeduction::GetDotListDim<DotList>::value;
    using Enable = DotSetPermMaskFlag<MaskType, dim, EnableFlag>;
    using Disable = DotSetPermMaskFlag<MaskType, dim, DisableFlag>;

    // permable     = PF_isPermeable == true    -> 1
    // not permable = PF_isPermeable == false   -> 0
    plb::applyProcessingFunctional(new Enable(PF_isPermeable), permNodes, permMask);
    plb::applyProcessingFunctional(new Disable(PF_isPermeable), nonPermNodes, permMask);
}


} // end of namespace PermFlagUtils

} // end of namespace IPP

#endif // PERMFLAGSUTILS_H
