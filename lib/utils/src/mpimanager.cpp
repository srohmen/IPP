#include "mpimanager.h"

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

namespace IPP
{

MPIManager MPIManager::m_instance;

MPIManager::MPIManager()
    : m_mpiEnv()
    , m_size(-1)
    , m_rank(-1)
    , m_isMainProc(false)
{

}

MPIManager::~MPIManager()
{

}

MPIManager& MPIManager::getInstance()
{
    return m_instance;
}

void MPIManager::init(int& argc, char**& argv, const boost::mpi::communicator& comm)
{
    m_comm = comm;
    boost::mpi::environment* env = new boost::mpi::environment(argc, argv);
    m_mpiEnv.reset(env);

    m_size = m_comm.size();
    m_rank = m_comm.rank();

    m_isMainProc = m_rank == 0;
}

boost::mpi::communicator MPIManager::getCommunicator() const
{
    return m_comm;
}

int MPIManager::getNProcs() const
{
    return m_size;
}

int MPIManager::getRank() const
{
    return m_rank;
}

bool MPIManager::isMainProc() const
{
    return m_isMainProc;
}

}
