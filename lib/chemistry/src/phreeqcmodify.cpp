#include "phreeqcmodify.h"

#include <boost/mpi.hpp>
#include <PhreeqcRM.h>

#include "ippexception.h"
#include "cellid_si_func.h"
#include "cellid_dissolveonly_func.h"
#include "mpitools.h"
#include "bench_tools.h"


namespace IPP
{

template<typename Func>
void collectAllData(const Func& func,
                    std::vector<typename Func::return_type>& newData)
{
    typedef typename Func::return_type return_type;

    const size_t nCells = func.nCells();


    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {
        std::vector<return_type> newValues;
        func.evaluate(iCell, newValues.begin(), newValues.end());

        for(size_t iPhase = 0; iPhase < newValues.size(); ++iPhase)
        {
            const return_type newSI = newValues[iPhase];
            newData.push_back(newSI);
        }
    }
}

template<typename Func>
void collectData(const Func& func,
                 const size_t nPhases,
                 const PhreeqcRM& phreeqc,
                 std::vector<typename Func::return_type>& newData)
{
    const size_t myTask = phreeqc.GetMpiMyself();
    assert(myTask >= 0);

    typedef typename Func::return_type return_type;
    std::vector<return_type> globaleNewData;
    if(myTask == 0)
    {
        BENCH(collectAllData(func, globaleNewData));
    }

    // distribute data according to cell IDs to which the current phreeqc process is responsible

    const int nTasks = phreeqc.GetMpiTasks();

    std::vector<int> displs(nTasks);
    std::vector<int> scounts(nTasks);

    for(int iTask = 0; iTask < nTasks; ++iTask)
    {
        const int startCell = phreeqc.GetStartCell().at(iTask);
        const int endCell = phreeqc.GetEndCell().at(iTask) + 1;
        const int diff = endCell - startCell;

        displs[iTask] = startCell * nPhases;
        scounts[iTask] = diff * nPhases;

    }

    const int startCell = phreeqc.GetStartCell().at(myTask);
    const int endCell = phreeqc.GetEndCell().at(myTask) + 1;
    const int diff = endCell - startCell;

    const MPI_Datatype mpiType = boost::mpi::get_mpi_datatype<return_type>();

    const int rcounts = diff * nPhases;
    newData.resize(rcounts);
    MPI_Scatterv(globaleNewData.data(), scounts.data(), displs.data(), mpiType,
                 newData.data(), rcounts, mpiType,
                 0, MPI_COMM_WORLD);


}

void PhreeqcModify::updateTargetSaturationIndices(const CellID_SI_Func& siFunc,
                                                  const std::vector<std::string>& phaseNames,
                                                  const size_t nCells,
                                                  PhreeqcRM& phreeqc,
                                                  std::vector<double>& newSImap)
{
    const size_t nPhases = phaseNames.size();

    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {
        const std::vector<double>::iterator begin = newSImap.begin() + iCell * nPhases;
        const std::vector<double>::iterator end = begin + nPhases;
        siFunc.evaluate(iCell, begin, end);
    }

    phreeqc.setTargetSaturationIndices(phaseNames, newSImap);
}

void PhreeqcModify::updateDissolveOnlyFlags(const CellID_DissolveOnly_Func &doFunc,
                                            const std::vector<std::string> &phaseNames,
                                            const size_t nCells,
                                            PhreeqcRM &phreeqc,
                                            std::vector<unsigned char>& currDissolveOnlyCells)
{
    const size_t nPhases = phaseNames.size();

    using Flags = std::vector<unsigned char>;
    Flags dissolveOnlyFlags(nCells * nPhases);

    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {
        const Flags::iterator begin = dissolveOnlyFlags.begin() + iCell * nPhases;
        const Flags::iterator end = begin + nPhases;
        const bool generalDissolveOnly = doFunc.evaluate(iCell, begin, end);
        currDissolveOnlyCells[iCell] = generalDissolveOnly;
    }

    phreeqc.setDissolveOnlyFlags(phaseNames, dissolveOnlyFlags);
}

} // end of namespace IPP
