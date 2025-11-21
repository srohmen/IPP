#include "palabosmpiutils.h"

#include "dimsizecalc.h"

namespace IPP
{

namespace PalabosMPIUtils
{

template<size_t nComp>
void gatherAllIndices(const ProcInfo& procInfo,
                      plb::plint globalIndices[nComp],
                      std::vector<plb::plint>& allGlobalIndices)
{
    if(procInfo.first)
    {
        allGlobalIndices.resize(procInfo.second * nComp);
    }

    const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<plb::plint>();
    MPI_Gather(globalIndices, nComp, mpiType,
               allGlobalIndices.data(), nComp, mpiType,
               0, plb::global::mpi().getGlobalCommunicator());

}

void gatherAllIndices(const PalabosMPIUtils::ProcInfo& procInfo,
                      const plb::Dot2D& offset, const plb::Box2D& domain,
                      std::vector<plb::plint>& allGlobalIndices)
{
    constexpr size_t nComp = 4;
    plb::plint globalIndices[nComp] = { offset.x + domain.x0,
                                        offset.x + domain.x1,
                                        offset.y + domain.y0,
                                        offset.y + domain.y1
                                      };

    gatherAllIndices<nComp>(procInfo, globalIndices, allGlobalIndices);
}

void gatherAllIndices(const PalabosMPIUtils::ProcInfo& procInfo,
                      const plb::Dot3D& offset, const plb::Box3D& domain,
                      std::vector<plb::plint>& allGlobalIndices)
{
    constexpr size_t nComp = 6;
    plb::plint globalIndices[nComp] = { offset.x + domain.x0,
                                        offset.x + domain.x1,
                                        offset.y + domain.y0,
                                        offset.y + domain.y1,
                                        offset.z + domain.z0,
                                        offset.z + domain.z1
                                      };

    gatherAllIndices<nComp>(procInfo, globalIndices, allGlobalIndices);
}


template<size_t nComp>
void allGatherAllIndices(const ProcInfo& procInfo,
                         plb::plint globalIndices[nComp],
                         std::vector<plb::plint>& allGlobalIndices)
{
    allGlobalIndices.resize(procInfo.second * nComp);

    const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<plb::plint>();
    MPI_Allgather(globalIndices, nComp, mpiType,
                  allGlobalIndices.data(), nComp, mpiType,
                  plb::global::mpi().getGlobalCommunicator());

}

void allGatherAllIndices(const PalabosMPIUtils::ProcInfo& procInfo,
                         const plb::Dot2D& offset, const plb::Box2D& domain,
                         std::vector<plb::plint>& allGlobalIndices)
{
    constexpr size_t nComp = 4;
    plb::plint globalIndices[nComp] = { offset.x + domain.x0,
                                        offset.x + domain.x1,
                                        offset.y + domain.y0,
                                        offset.y + domain.y1
                                      };

    allGatherAllIndices<nComp>(procInfo, globalIndices, allGlobalIndices);
}

void allGatherAllIndices(const PalabosMPIUtils::ProcInfo& procInfo,
                         const plb::Dot3D& offset, const plb::Box3D& domain,
                         std::vector<plb::plint>& allGlobalIndices)
{
    constexpr size_t nComp = 6;
    plb::plint globalIndices[nComp] = { offset.x + domain.x0,
                                        offset.x + domain.x1,
                                        offset.y + domain.y0,
                                        offset.y + domain.y1,
                                        offset.z + domain.z0,
                                        offset.z + domain.z1
                                      };

    allGatherAllIndices<nComp>(procInfo, globalIndices, allGlobalIndices);
}


template<typename SizeCalc>
void transferStrides(const std::vector<plb::plint>& allGlobalIndices,
                     const size_t nProcesses,
                     std::vector<int>& recvcounts,
                     std::vector<int>& displs)
{
    recvcounts.resize(nProcesses);
    displs.resize(nProcesses);

    size_t prevStride = 0;
    for (size_t iProc = 0; iProc < nProcesses; ++iProc)
    {
        const size_t startIndex = iProc * SizeCalc::nComp;
        const size_t size = SizeCalc::calc(allGlobalIndices, startIndex);

        recvcounts[iProc] = size;
        displs[iProc] = prevStride;
        prevStride += size;
    }
}

template
void transferStrides<DimSizeCalc2D>(
const std::vector<plb::plint>& allGlobalIndices,
const size_t nProcesses,
std::vector<int>& recvcounts,
std::vector<int>& displs);

template
void transferStrides<DimSizeCalc3D>(
const std::vector<plb::plint>& allGlobalIndices,
const size_t nProcesses,
std::vector<int>& recvcounts,
std::vector<int>& displs);



}


}
