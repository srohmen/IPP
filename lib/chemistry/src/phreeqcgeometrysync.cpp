#include "phreeqcgeometrysync.h"

#include "arraydimensionconvert.h"
#include "phreeqcinitialdata.h"
#include "phreeqcconstants.h"
#include "indexhelper.h"
#include "calcporositychangefactor.h"

#include "mpi_types.h"
#include <algorithm>
#include "mpimanager.h"

#include "nodeinfos.h"


namespace IPP
{

PhreeqcGeometrySync::PhreeqcGeometrySync(const IPPVector3DLong &localOrigin,
                                         const ArrayDimensionConvert &indexConvLocal,
                                         const PhreeqcInitialData &initialData,
                                         NodeInfos &nodeInfos,
                                         NeighborInfos &neighInfos)
    : m_printDebug(false)
    , m_localOrigin(localOrigin)
    , m_indexConvLocal(indexConvLocal)
    , m_initialData(initialData)
    , m_interfaceCells(nodeInfos.interfaceCells)
    , m_nonPermNonInterfaceCells(nodeInfos.nonPermNonInterfaceCells)
    , m_nInterfaceNodes(nodeInfos.nInterfaceNodesGlobal)
    , m_neighInfos(neighInfos)
{

}

void PhreeqcGeometrySync::run(const NodeList &newPermNodes,
                              const NodeList &newInterfaceNodes,
                              const NodeList &newNonPermNonInterfaceNodes,
                              std::vector<double> &conc)
{

    for(const IPPVector3DLong& node : newNonPermNonInterfaceNodes)
    {
        if(m_printDebug)
        {
            std::cout << MPIManager::getInstance().getRank() <<  " newNonPermNonInterfaceNodes: "
                      << node[0] << "\t" << node[1] << "\t" << node[2] << std::endl;
        }

        const IPPVector3DLong localCoord = node - m_localOrigin;
        const size_t iCell = m_indexConvLocal.calcIndex(localCoord);

        assert(m_nonPermNonInterfaceCells.size() > iCell);
        m_nonPermNonInterfaceCells[iCell] = true;

        assert(m_interfaceCells.size() > iCell);
        m_interfaceCells[iCell] = false;
    }

    for(const IPPVector3DLong& node : newInterfaceNodes)
    {
        if(m_printDebug)
        {
            std::cout << MPIManager::getInstance().getRank() <<  " newInterfaceNodes: "
                      << node[0] << "\t" << node[1] << "\t" << node[2] << std::endl;
        }

        const IPPVector3DLong localCoord = node - m_localOrigin;
        const size_t iCell = m_indexConvLocal.calcIndex(localCoord);

        assert(m_nonPermNonInterfaceCells.size() > iCell);
        m_nonPermNonInterfaceCells[iCell] = false;

        assert(m_interfaceCells.size() > iCell);
        m_interfaceCells[iCell] = true;
    }


    const size_t nCells = m_indexConvLocal.getNxyz();

    for(const IPPVector3DLong& node : newPermNodes)
    {
        if(m_printDebug)
        {
            std::cout << MPIManager::getInstance().getRank() <<  " newPermNodes: "
                      << node[0] << "\t" << node[1] << "\t" << node[2] << std::endl;
        }

        const IPPVector3DLong localCoord = node - m_localOrigin;
        const size_t iCell = m_indexConvLocal.calcIndex(localCoord);

        assert(m_nonPermNonInterfaceCells.size() > iCell);
        m_nonPermNonInterfaceCells[iCell] = false;

        assert(m_interfaceCells.size() > iCell);
        m_interfaceCells[iCell] = false;

        for(size_t iComp = 0; iComp < PhreeqcConstants::s_primaryCompsBegin; ++iComp)
        {
            const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, iCell, iComp);
            assert(m_initialData.conc.size() > index);
            const double& initConc = m_initialData.conc[index];

            assert(conc.size() > index);
            conc[index] = initConc;
        }
    }


    size_t nInterfaceNodesLocal = std::count(m_interfaceCells.begin(), m_interfaceCells.end(), true);

    MPI_Allreduce(&nInterfaceNodesLocal, &m_nInterfaceNodes, 1,
                  MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

}

bool PhreeqcGeometrySync::needDistField() const
{
    // TODO: make depending on SI calc
    return false;
}

bool PhreeqcGeometrySync::needNeighMinPoros() const
{
    // TODO: find better check
    return m_neighInfos.results.empty() == false;
}

AbstractGeometrySync::NeighborInfos& PhreeqcGeometrySync::getNeighInfos()
{
    return m_neighInfos;
}

void PhreeqcGeometrySync::setPrintDebug(const bool printDebug)
{
    m_printDebug = printDebug;
}

}
