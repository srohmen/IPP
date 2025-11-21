#ifndef PALABOSMULTIDATAWRAPPER_H
#define PALABOSMULTIDATAWRAPPER_H

#include "mpimanager.h"
//#include "palabosmpiutils.h"
#include "fielddecomposition.h"
#include "arraydimensionconvert.h"
#include "indexhelper.h"

namespace IPP
{
/**
 * upon creation the wrapper class does determine the local bounds.
 * when a raw data field is set it is decomposed accordingly and transfered to the
 * single processes.
 */



template<typename T>
class SerialToMultiDataSync
{
public:
    SerialToMultiDataSync(std::vector<T>& localData,
                          const FieldDecomposition& decomp,
                          const size_t vectorDim = 1)
        : m_localData(localData)
        , m_decomp(decomp)
        , m_vectorDim(vectorDim)
    {

    }

    void pushFromRoot(const std::vector<T>& globalRawData)
    {
        const std::vector<int>& displs = m_decomp.getDispls();
        const std::vector<int>& sndCounts = m_decomp.getSndCounts();
        const size_t localSize = m_decomp.getLocalNumberOfCells();

        // std::cout << "localSize: " << localSize << std::endl;

        // copy global data to consecutive format
        std::vector<T> globalArr;
        if(MPIManager::getInstance().isMainProc())
        {
            copyToConsecutive(globalRawData, globalArr);
        }

        // copy data to child procs
        const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<T>();
        m_localData.resize(localSize * m_vectorDim);

        MPI_Scatterv(globalArr.data(), sndCounts.data(), displs.data(), mpiType,
                     m_localData.data(), m_localData.size(), mpiType,
                     0, MPI_COMM_WORLD);

    }

private:
    void copyToConsecutive(const std::vector<T>& inputArr,
                           std::vector<T>& outputArr)
    {
        const size_t nCellsGlobal = m_decomp.getGlobalNumberCells();
        const size_t nValues = nCellsGlobal * m_vectorDim;
        assert(inputArr.size() == nValues);

        outputArr.resize(nValues);
        typename std::vector<T>::iterator out = outputArr.begin();

        const IPPVector3DLong& globalSize = m_decomp.getGlobalSize();
        const ArrayDimensionConvert indexConv(globalSize);


        const std::vector<IPPBox3DLong>& localDomains = m_decomp.getDomains();

        for(size_t iDomain = 0; iDomain < localDomains.size(); ++iDomain)
        {
            const IPPBox3DLong& domain = localDomains[iDomain];

            for(size_t iComp = 0; iComp < m_vectorDim; ++iComp)
            {
                for(long xSrc = domain.lower[0]; xSrc <= domain.upper[0]; ++xSrc)
                {
                    for(long ySrc = domain.lower[1]; ySrc <= domain.upper[1]; ++ySrc)
                    {
                        for(long zSrc = domain.lower[2]; zSrc <= domain.upper[2]; ++zSrc)
                        {
                            const size_t iCell = indexConv.calcIndex(xSrc, ySrc, zSrc);
                            const size_t index = IndexHelper::toPhreeqcRawIndex(nCellsGlobal, iCell, iComp);

                            assert(inputArr.size() > index);
                            const T& val = inputArr[index];
                            *out = val;

                            ++out;
                        }
                    }
                }
            }
        }
    }


    std::vector<T>& m_localData;
    const FieldDecomposition& m_decomp;
    const size_t m_vectorDim; // elements per cell
};

}


#endif // PALABOSMULTIDATAWRAPPER_H
