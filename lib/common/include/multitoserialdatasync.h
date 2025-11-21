#ifndef MULTITOSERIALDATASYNC_H
#define MULTITOSERIALDATASYNC_H

#include "mpimanager.h"
#include "fielddecomposition.h"
#include "arraydimensionconvert.h"

namespace IPP
{

template<typename T>
class MultiToSerialDataSync
{
public:
    MultiToSerialDataSync(std::vector<T>& globalData,
                          const FieldDecomposition& decomp,
                          const size_t vectorDim = 1)
        : m_globalData(globalData)
        , m_decomp(decomp)
        , m_vectorDim(vectorDim)
    {

    }

    void pullToRoot(const std::vector<T>& localRawData)
    {
        const std::vector<int>& displs = m_decomp.getDispls();
        const std::vector<int>& sndCounts = m_decomp.getSndCounts();

        std::vector<T> globalRawData;
        if(MPIManager::getInstance().isMainProc())
        {
            const size_t nCellsTotal = m_decomp.getGlobalNumberCells();
            globalRawData.resize(nCellsTotal * m_vectorDim);
        }

        // copy data to root proc
        const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<T>();

        MPI_Gatherv(localRawData.data(), localRawData.size(), mpiType,
                    globalRawData.data(), sndCounts.data(), displs.data(), mpiType,
                    0, MPI_COMM_WORLD );

        // copy consecutive format to global arr
        if(MPIManager::getInstance().isMainProc())
        {
            copyToShape(globalRawData, m_globalData);
        }
    }

private:
    void copyToShape(const std::vector<T>& inputArr,
                     std::vector<T>& outputArr)
    {
        outputArr.resize(inputArr.size());
        typename std::vector<T>::const_iterator srcIt = inputArr.begin();

        const IPPVector3DLong& globalSize = m_decomp.getGlobalSize();
        const ArrayDimensionConvert indexConv(globalSize);

        const size_t nCellsGlobal = m_decomp.getGlobalNumberCells();


        const std::vector<IPPBox3DLong>& localDomains = m_decomp.getDomains();

        for(size_t iDomain = 0; iDomain < localDomains.size(); ++iDomain)
        {
            const IPPBox3DLong& domain = localDomains[iDomain];
            for(size_t iVecDim = 0; iVecDim < m_vectorDim; ++iVecDim)
            {
                const size_t vecDimOff = iVecDim * nCellsGlobal;

                for(long xDst = domain.lower[0]; xDst <= domain.upper[0]; ++xDst)
                {
                    for(long yDst = domain.lower[1]; yDst <= domain.upper[1]; ++yDst)
                    {
                        for(long zDst = domain.lower[2]; zDst <= domain.upper[2]; ++zDst)
                        {
                            const size_t dstCellIndex = indexConv.calcIndex(xDst, yDst, zDst);
                            const size_t dstIndex = dstCellIndex + vecDimOff;
                            outputArr[dstIndex] = *srcIt;
                            ++srcIt;
                        }
                    }
                }
            }
        }
    }


    std::vector<T>& m_globalData;
    const FieldDecomposition& m_decomp;
    const size_t m_vectorDim; // elements per cell


};

}


#endif // MULTITOSERIALDATASYNC_H
