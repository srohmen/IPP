#ifndef CELLINDEXCORRECTION_H
#define CELLINDEXCORRECTION_H

#include <vector>
#include "ippbox.h"
#include "arraydimensionconvert.h"

namespace IPP
{

class CellIndexCorrection
{
public:
    CellIndexCorrection(const IPPVector3DInt &oldSize,
                        const IPPVector3DInt &newSize,
                        const IPPVector3DInt &oldOrigin,
                        const IPPBox3DLong& localDomain,
                        const std::vector<size_t>& decompCorr);

    size_t correctForBCs(const size_t orig) const;
    size_t correctForDecomp(const size_t orig) const;
    size_t convertToGlobal(const size_t iCellLocal) const;
    size_t convertToLocal(const size_t iCellGlobal) const;

private:
    const ArrayDimensionConvert m_oldGlobal;
    const ArrayDimensionConvert m_newGlobal;
    const IPPVector3DInt m_oldOriginGlobal;

    const ArrayDimensionConvert m_localConv;
    const IPPVector3DLong m_ownDomainOrigin;

    const std::vector<size_t> m_decompCorr;

};


}


#endif // CELLINDEXCORRECTION_H

