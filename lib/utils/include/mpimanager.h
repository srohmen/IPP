#ifndef MPIMANAGER_H
#define MPIMANAGER_H

#include <memory>

#include <boost/mpi/communicator.hpp>

namespace boost {
namespace mpi {
class environment;
}
}

namespace IPP
{

class MPIManager
{
public:
    static MPIManager& getInstance();

    void init(int& argc, char** &argv, const boost::mpi::communicator& comm);

    boost::mpi::communicator getCommunicator() const;
    int getNProcs() const;
    int getRank() const;
    bool isMainProc() const;

private:
    MPIManager();
    ~MPIManager();

    static MPIManager m_instance;

    boost::mpi::communicator m_comm;
    std::unique_ptr<boost::mpi::environment> m_mpiEnv;

    int m_size;
    int m_rank;

    bool m_isMainProc;


};

}

#endif // MPIMANAGER_H
