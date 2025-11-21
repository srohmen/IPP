#ifndef PHREEQCGEOMETRYSYNC_H
#define PHREEQCGEOMETRYSYNC_H

#include "abstractgeometrysync.h"

namespace IPP
{

class ArrayDimensionConvert;
struct PhreeqcInitialData;
struct NodeInfos;

class PhreeqcGeometrySync : public AbstractGeometrySync
{
public:
    PhreeqcGeometrySync(const IPPVector3DLong& localOrigin,
                        const ArrayDimensionConvert& indexConvLocal,
                        const PhreeqcInitialData& initialData,
                        NodeInfos& nodeInfos,
                        NeighborInfos& neighInfos);


    virtual void run(const NodeList& newPermNodes,
                     const NodeList& newInterfaceNodes,
                     const NodeList& newNonPermNonInterfaceNodes,
                     std::vector<double>& conc) override;

    virtual bool needDistField() const override;
    virtual bool needNeighMinPoros() const override;
    virtual NeighborInfos &getNeighInfos() override;

    virtual void setPrintDebug(const bool printDebug) override;

private:
    bool m_printDebug;

    const IPPVector3DLong& m_localOrigin;
    const ArrayDimensionConvert& m_indexConvLocal;
    const PhreeqcInitialData& m_initialData;

    std::vector<char>& m_interfaceCells;
    std::vector<char>& m_nonPermNonInterfaceCells;
    size_t &m_nInterfaceNodes;

    NeighborInfos m_neighInfos;
};


}

#endif // PHREEQCGEOMETRYSYNC_H
