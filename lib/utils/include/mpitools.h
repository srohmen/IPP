#ifndef MPITOOLS_H
#define MPITOOLS_H

#include <boost/mpi/communicator.hpp>
#include <vector>
#include <string>

namespace IPP
{

namespace MpiTools
{

template<typename T>
void syncFromRank(const boost::mpi::communicator &comm,
                  int srcRank,
                  std::vector<T> &vec)
{
    size_t nElements = vec.size();
    MPI_Datatype mpiSizeType = boost::mpi::get_mpi_datatype(nElements);
    MPI_Bcast(&nElements, 1, mpiSizeType, srcRank, comm);

    if(comm.rank() != srcRank)
    {
        vec.resize(nElements);
    }

    MPI_Datatype mpiContainerType = boost::mpi::get_mpi_datatype<T>();
    MPI_Bcast(vec.data(), nElements, mpiContainerType, srcRank, comm);
}


void syncFromRoot(const boost::mpi::communicator& comm,
                  const std::vector<char>& vec);

void syncFromRoot(const boost::mpi::communicator& comm,
                  const std::vector<double>& vec);

void syncFromRoot(const boost::mpi::communicator& comm,
                  const std::vector<std::string>& strVec);

void syncFromRoot(const boost::mpi::communicator& comm,
                  const std::vector< std::pair<size_t, size_t> > &map);

void checkSync(MPI_Comm rawComm,
               const std::string& fileName,
               const std::string& functionName,
               const int line);

void printVector(const int myRank, const std::vector<double>& conc, const std::string& prefix = "");


template<typename Func>
void runFuncSerialized(const boost::mpi::communicator& comm, Func& func)
{
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    for(int i = 0; i < rank; ++i)
    {
        MPI_Barrier(comm);
    }

    if(rank == 0)
    {
        std::cout << "--------------snip--------------" << std::endl
                  << "rank: " << rank << std::endl;
    }
    else
    {
        std::cout << "rank: " << rank << std::endl;
    }

    func();
    std::cout << std::endl;

    for(int i = 0; i < (size - rank); ++i)
    {
        MPI_Barrier(comm);
    }

    if(rank == 0)
    {
        std::cout << "--------------snap--------------" << std::endl;
    }
}

}

} // end of namespace LBGeoChem

#ifdef ENABLE_MPI_CHECK_SYNC
#define MPI_CHECK_SYNC IPP::MpiTools::checkSync(MPI_COMM_WORLD, __FILE__, __func__, __LINE__);
#else
#define MPI_CHECK_SYNC ;
#endif

#endif // MPITOOLS_H

