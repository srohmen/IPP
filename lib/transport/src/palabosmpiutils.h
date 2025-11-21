#ifndef PALABOSMPIUTILS_H
#define PALABOSMPIUTILS_H

#include <mpi.h>
#include <boost/mpi.hpp>

#include <palabos/core/globalDefs.h>
#include <palabos/core/geometry2D.h>
#include <palabos/core/geometry3D.h>

#include "mpitools.h"

namespace IPP
{

namespace PalabosMPIUtils
{

typedef std::pair<bool, size_t> ProcInfo;

void gatherAllIndices(const ProcInfo& procInfo,
                      const plb::Dot2D& offset,
                      const plb::Box2D& domain,
                      std::vector<plb::plint>& allGlobalIndices);

void gatherAllIndices(const ProcInfo& procInfo,
                      const plb::Dot3D& offset,
                      const plb::Box3D& domain,
                      std::vector<plb::plint>& allGlobalIndices);

void allGatherAllIndices(const PalabosMPIUtils::ProcInfo& procInfo,
                         const plb::Dot2D& offset, const plb::Box2D& domain,
                         std::vector<plb::plint>& allGlobalIndices);

void allGatherAllIndices(const PalabosMPIUtils::ProcInfo& procInfo,
                         const plb::Dot3D& offset, const plb::Box3D& domain,
                         std::vector<plb::plint>& allGlobalIndices);


template<typename SizeCalc>
void transferStrides(const std::vector<plb::plint>& allGlobalIndices,
                     const size_t nProcesses,
                     std::vector<int>& recvcounts,
                     std::vector<int>& displs);

}

}

#endif // PALABOSMPIUTILS_H
