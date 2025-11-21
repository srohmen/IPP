#ifndef SYNCPERIODICDOMAIN_H
#define SYNCPERIODICDOMAIN_H

#include <boost/mpi/communicator.hpp>

namespace plb
{
struct DotList2D;
struct DotList3D;

struct Box2D;
struct Box3D;
}

namespace IPP
{

class SyncPeriodicDomain
{
public:
    SyncPeriodicDomain(boost::mpi::communicator comm);

    typedef std::vector<int> Periodicity;
    void sync(const plb::Box2D& boundingBox,
              const Periodicity& periodic, plb::DotList2D& dots);
    void sync(const plb::Box3D& boundingBox,
              const Periodicity& periodic, plb::DotList3D& dots);

private:
    template<size_t dim, typename Surface, typename DotList>
    void sync(const Surface& surf, const Periodicity& periodic, DotList& dots);

    template<size_t dim, typename DotList>
    void syncFromOthers(const size_t nElements,
                        const DotList& periodicRelevant,
                        DotList& others);

    boost::mpi::communicator m_comm;
};

}

#endif // SYNCPERIODICDOMAIN_H

