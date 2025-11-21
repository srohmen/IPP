#include "syncperiodicdomain.h"


#include <palabos/core/geometry2D.h>
#include <palabos/core/geometry3D.h>


#include "plbtypededuction.h"
#include "ippexception.h"

namespace IPP
{

template<typename DotList>
DotList append(const DotList& a, const DotList& b)
{
    DotList merged = a;
    for(int i = 0; i < b.getN(); ++i )
    {
        merged.addDot(b.getDot(i));
    }

    return merged;
}

void findCorners(const plb::BlockSurface2D& surf, std::vector<plb::Box2D>& boxes)
{
    boxes.push_back(surf.cornerNN());
    boxes.push_back(surf.cornerPN());
    boxes.push_back(surf.cornerNP());
    boxes.push_back(surf.cornerPP());
}

std::vector<plb::Box2D> findBoundaryBoxes(const plb::BlockSurface2D& surf, const int periodicDim)
{
    std::vector<plb::Box2D> boxes;

    if(periodicDim == 0)
    {
        boxes.push_back(surf.edge0N());
        boxes.push_back(surf.edge0P());
        findCorners(surf, boxes);
    }
    else if(periodicDim == 1)
    {
        boxes.push_back(surf.edge1N());
        boxes.push_back(surf.edge1P());
        findCorners(surf, boxes);
    }
    else
    {
        throw std::runtime_error("invalid dimension: " + std::to_string(periodicDim));
    }

    return boxes;
}


void findCorners(const plb::BlockSurface3D& surf, std::vector<plb::Box3D>& boxes)
{
    boxes.push_back(surf.cornerNNN());
    boxes.push_back(surf.cornerNNP());
    boxes.push_back(surf.cornerNPN());
    boxes.push_back(surf.cornerNPP());
    boxes.push_back(surf.cornerPNN());
    boxes.push_back(surf.cornerPNP());
    boxes.push_back(surf.cornerPPN());
    boxes.push_back(surf.cornerPPP());
}

std::vector<plb::Box3D> findBoundaryBoxes(const plb::BlockSurface3D& surf, const int periodicDim)
{
    std::vector<plb::Box3D> boxes;

    if(periodicDim == 0)
    {
        boxes.push_back(surf.surface0N());
        boxes.push_back(surf.surface0P());

        boxes.push_back(surf.edge0NN());
        boxes.push_back(surf.edge0NP());
        boxes.push_back(surf.edge0PN());
        boxes.push_back(surf.edge0PP());

        findCorners(surf, boxes);
    }
    else if(periodicDim == 1)
    {
        boxes.push_back(surf.surface1N());
        boxes.push_back(surf.surface1P());

        boxes.push_back(surf.edge1NN());
        boxes.push_back(surf.edge1NP());
        boxes.push_back(surf.edge1PN());
        boxes.push_back(surf.edge1PP());

        findCorners(surf, boxes);
    }
    else if(periodicDim == 2)
    {
        boxes.push_back(surf.surface2N());
        boxes.push_back(surf.surface2P());

        boxes.push_back(surf.edge2NN());
        boxes.push_back(surf.edge2NP());
        boxes.push_back(surf.edge2PN());
        boxes.push_back(surf.edge2PP());

        findCorners(surf, boxes);
    }
    else
    {
        throw std::runtime_error("invalid dimension: " + std::to_string(periodicDim));
    }

    return boxes;
}


template<size_t dim, typename Surface, typename DotList>
void findOthers(const DotList& own, const Surface& surf, const int periodicDim, DotList& periodicRelevant)
{
    const auto boundaryBoxes = findBoundaryBoxes(surf, periodicDim);

    for(const auto& box : boundaryBoxes)
    {
        DotList otherBoundaryDots;
        plb::intersect(box, own, otherBoundaryDots);
        periodicRelevant = append(periodicRelevant, otherBoundaryDots);
    }
}

typedef std::vector<plb::plint>::iterator VectorPlIntIterator;
void setValue(const VectorPlIntIterator begin, const VectorPlIntIterator end, const plb::Dot2D& dot)
{
    assert(end - begin == 2);
    *(begin + 0) = dot.x;
    *(begin + 1) = dot.y;
}

void setValue(const VectorPlIntIterator begin, const VectorPlIntIterator end, const plb::Dot3D& dot)
{
    assert(end - begin == 3);
    *(begin + 0) = dot.x;
    *(begin + 1) = dot.y;
    *(begin + 2) = dot.z;
}

template<size_t dim>
struct MakeDot;

template<>
struct MakeDot<2>
{
    static plb::Dot2D make(const VectorPlIntIterator begin, const VectorPlIntIterator end)
    {
        assert(end - begin == 2);
        return plb::Dot2D(*begin, *(begin + 1));
    }
};

template<>
struct MakeDot<3>
{
    static plb::Dot3D make(const VectorPlIntIterator begin, const VectorPlIntIterator end)
    {
        assert(end - begin == 3);
        return plb::Dot3D(*begin, *(begin + 1), *(begin + 2));
    }
};


plb::BlockSurface2D makeSurface(const plb::Box2D &domain)
{
    return plb::BlockSurface2D(domain, 1);
}

plb::BlockSurface3D makeSurface(const plb::Box3D &domain)
{
    return plb::BlockSurface3D(domain, 1);
}


SyncPeriodicDomain::SyncPeriodicDomain(boost::mpi::communicator comm)
    : m_comm(comm)
{

}


void SyncPeriodicDomain::sync(const plb::Box2D& boundingBox, const Periodicity& periodic, plb::DotList2D& dots)
{
    sync<2>(boundingBox, periodic, dots);
}

void SyncPeriodicDomain::sync(const plb::Box3D& boundingBox, const Periodicity& periodic, plb::DotList3D& dots)
{
    sync<3>(boundingBox, periodic, dots);
}


template<size_t dim, typename Box, typename DotList>
void SyncPeriodicDomain::sync(const Box& boundingBox, const Periodicity& periodic, DotList& dots)
{
    const auto surf = makeSurface(boundingBox);

    DotList periodicRelevant;

    for(size_t i = 0; i < periodic.size(); ++i)
    {
        int periodicDim = periodic[i];
        findOthers<dim>(dots, surf, periodicDim, periodicRelevant);
    }

    const size_t maxElement = boundingBox.getMaxWidth();

    DotList others;
    syncFromOthers<dim>(maxElement, periodicRelevant, others);

    dots = append(dots, others);

    std::sort(dots.dots.begin(), dots.dots.end());
    dots.dots.erase(std::unique(dots.dots.begin(), dots.dots.end()), dots.dots.end());
}

template<size_t dim, typename DotList>
void SyncPeriodicDomain::syncFromOthers(const size_t nDots,
                                        const DotList& periodicRelevant,
                                        DotList& others)
{
    static constexpr plb::plint invalidCoord = std::numeric_limits<plb::plint>::lowest();

    const size_t nElements = nDots * dim * 2;
    std::vector<plb::plint> sendBuf(nElements, invalidCoord);
    for(int iDot = 0 ; iDot < periodicRelevant.getN(); ++iDot)
    {
        const auto& dot = periodicRelevant.getDot(iDot);

        const VectorPlIntIterator begin = sendBuf.begin() + iDot * dim;
        const VectorPlIntIterator end = begin + dim;
        setValue(begin, end, dot);
    }

    int nProcs;
    MPI_Comm_size(m_comm, &nProcs);
    MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<plb::plint>();


    std::vector<plb::plint> recvBuf(sendBuf.size() * nProcs, invalidCoord);
    MPI_Allgather(sendBuf.data(), nElements, mpiType,
                  recvBuf.data(), nElements, mpiType,
                  m_comm);

    for(VectorPlIntIterator it = recvBuf.begin(); it != recvBuf.end(); it += dim)
    {
        const auto& dot = MakeDot<dim>::make(it, it + dim);
        others.addDot(dot);
    }


    std::vector<plb::plint> dummy(dim, invalidCoord);
    const auto invalidDot = MakeDot<dim>::make(dummy.begin(), dummy.end());
    others.dots.erase(std::remove(others.dots.begin(), others.dots.end(), invalidDot),
                      others.dots.end());


}

}
