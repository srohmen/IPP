#include "mpitools.h"

#include <iostream>

#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include "mpi_types.h"
#include "ippexception.h"

#include <sstream>

namespace IPP
{


namespace MpiTools
{



template<typename T>
void syncFromRootVec(const boost::mpi::communicator& comm, const std::vector<T> &vec)
{
    std::vector<T>& nonConstVec = const_cast<std::vector<T>&>(vec);

    size_t nElements = nonConstVec.size();
    MPI_Bcast(&nElements, 1, MPI_SIZE_T, 0, comm);

    if(comm.rank() != 0)
    {
        nonConstVec.resize(nElements);
    }

    const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<T>();

    MPI_Bcast(nonConstVec.data(), nElements, mpiType, 0, comm);
}

void syncFromRoot(const boost::mpi::communicator &comm, const std::vector<char> &vec)
{
    syncFromRootVec(comm, vec);
}

void syncFromRoot(const boost::mpi::communicator &comm, const std::vector<double> &vec)
{
    syncFromRootVec(comm, vec);
}

void syncFromRoot(const boost::mpi::communicator& comm,
                  const std::vector<std::string>& strVec)
{
    std::vector<std::string>& nonConstVec = const_cast<std::vector<std::string>&>(strVec);

    size_t nElements = nonConstVec.size();
    MPI_Bcast(&nElements, 1, MPI_SIZE_T, 0, comm);

    if(comm.rank() != 0)
    {
        nonConstVec.resize(nElements);

        for(size_t i = 0; i < nElements; ++i)
        {
            std::string& str = nonConstVec[i];

            size_t nChar;
            MPI_Bcast(&nChar, 1, MPI_SIZE_T, 0, comm);
            char *buf = new char[nChar];

            MPI_Bcast(buf, nChar, MPI_CHAR, 0, comm);
            str = std::string(buf, nChar);

            delete [] buf;
        }
    }
    else
    {
        for(size_t i = 0; i < nElements; ++i)
        {
            const std::string& str = strVec[i];

            size_t nChar = str.size();
            MPI_Bcast(&nChar, 1, MPI_SIZE_T, 0, comm);

            char* buf = const_cast<char*>(str.c_str());
            MPI_Bcast(buf, nChar, MPI_CHAR, 0, comm);
        }
    }


}

void syncFromRoot(const boost::mpi::communicator &comm,
                  const std::vector<std::pair<size_t, size_t>> &map)
{
    std::vector<std::pair<size_t, size_t>>& nonConstObj = const_cast< std::vector<std::pair<size_t, size_t>>& >(map);

    size_t nElements = nonConstObj.size();
    MPI_Bcast(&nElements, 1, MPI_SIZE_T, 0, comm);

    if(comm.rank() != 0)
    {
        nonConstObj.resize(nElements);
    }

    for(size_t i = 0; i < nElements; ++i)
    {
        std::pair<size_t, size_t>& element = nonConstObj[i];
        MPI_Bcast(&element.first, 1, MPI_SIZE_T, 0, comm);
        MPI_Bcast(&element.second, 1, MPI_SIZE_T, 0, comm);
    }

}

void checkSync(MPI_Comm rawComm,
               const std::string& fileName,
               const std::string& functionName, const int line)
{
    boost::mpi::communicator comm(rawComm, boost::mpi::comm_attach);
    comm.barrier();

    if(comm.rank() == 0)
    {
        std::vector<std::string> allFiles;
        boost::mpi::gather(comm, fileName, allFiles, 0);

        std::vector<std::string> allFunctions;
        boost::mpi::gather(comm, functionName, allFunctions, 0);

        std::vector<int> all_numbers;
        boost::mpi::gather(comm, line, all_numbers, 0);

        for (int proc = 0; proc < comm.size(); ++proc)
        {
            const std::string& otherFileName = allFiles[proc];
            const std::string& otherFunctionName = allFunctions[proc];
            const int otherLine = all_numbers[proc];

            if(otherFileName != fileName || otherFunctionName != functionName || otherLine != line)
            {
                std::string errStr = "process not in sync:\n"
                        + std::to_string(0) + " vs. " + std::to_string(proc) + "\n"
                        + fileName + " vs. " + otherFileName + "\n"
                        + functionName + " vs. " + otherFunctionName + "\n"
                        + std::to_string(line) + " vs. " + std::to_string(otherLine) + "\n";
                std::cerr << errStr;
                //throw std::runtime_error(errStr);

                while(true)
                {

                }
            }
        }
    }
    else
    {
        boost::mpi::gather(comm, fileName, 0);
        boost::mpi::gather(comm, functionName, 0);
        boost::mpi::gather(comm, line, 0);
    }

    comm.barrier();
}

void printVector(const int myRank, const std::vector<double>& conc, const std::string& prefix)
{
    std::stringstream ss;
    ss << prefix << "\t" << myRank << "\t";
    for(size_t i = 64; i < conc.size(); ++i)
    {
        ss << conc[i] << " ";
    }
    ss << std::endl;
    std::cerr << ss.str();
}

} // end of namespace MpiTools
} // end of namespace IPP
