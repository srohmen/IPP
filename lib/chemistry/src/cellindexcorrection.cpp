#include "cellindexcorrection.h"

namespace IPP
{

CellIndexCorrection::CellIndexCorrection(const IPPVector3DInt &oldSize,
                                         const IPPVector3DInt &newSize,
                                         const IPPVector3DInt &oldOrigin,
                                         const IPPBox3DLong& localDomain,
                                         const std::vector<size_t> &decompCorr)
    : m_oldGlobal(oldSize)
    , m_newGlobal(newSize)
    , m_oldOriginGlobal(oldOrigin)
    , m_localConv(localDomain.getSize())
    , m_ownDomainOrigin(localDomain.lower)
    , m_decompCorr(decompCorr)
{

}

size_t CellIndexCorrection::correctForBCs(const size_t orig) const
{
    IPPVector3DInt coord = m_oldGlobal.calcCoordinate(orig);
    coord = coord - m_oldOriginGlobal;
    const size_t newIndex = m_newGlobal.calcIndex(coord);
    return newIndex;
}

size_t CellIndexCorrection::correctForDecomp(const size_t orig) const
{
    const size_t corrected = m_decompCorr.at(orig);
    return corrected;
}

size_t CellIndexCorrection::convertToGlobal(const size_t iCellLocal) const
{
    const IPPVector3DInt localCoord = m_localConv.calcCoordinate(iCellLocal);
    const IPPVector3DInt globalCoord = localCoord + m_ownDomainOrigin;
    const size_t globalIndex = m_newGlobal.calcIndex(globalCoord);
    return globalIndex;
}

size_t CellIndexCorrection::convertToLocal(const size_t iCellGlobal) const
{
    const IPPVector3DInt globalCoord = m_newGlobal.calcCoordinate(iCellGlobal);
    const IPPVector3DInt localCoord = globalCoord - m_ownDomainOrigin;
    const size_t localIndex = m_localConv.calcIndex(localCoord);
    return localIndex;
}

}
