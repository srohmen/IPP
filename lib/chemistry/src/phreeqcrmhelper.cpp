#include "phreeqcrmhelper.h"

#include <fstream>
#include <numeric>

#include <PhreeqcRM.h>
#include <IPhreeqcPhast.h>

#include "ippexception.h"
#include "reactionmoduleconfig.h"
#include "ippconfig.h"
#include "phreeqcglobalinitdata.h"
#include "cellindexcorrection.h"
#include "mpitools.h"
#include "mpimanager.h"
#include "ippstream.h"

namespace IPP
{


void PhreeqcRMHelper::muteConsoleOutput(PhreeqcRM& phreeqcRM, const boost::filesystem::path& outDir)
{
    // silence console output due to performance reasons
    phreeqcRM.SetScreenOn(false);

    // do not print warning messages to std::cerr
    const int myRank = phreeqcRM.GetMpiMyself();
    const std::vector<IPhreeqcPhast*>& workers = phreeqcRM.GetWorkers();

    for(size_t i = 0; i < workers.size(); ++i)
    {
        // mem leak?
        std::ofstream* errStrm = new std::ofstream;
        const boost::filesystem::path outFile = outDir / ("phreeqc_err_" + std::to_string(myRank)
                                                          + "_" + std::to_string(i));
        const std::string outFileStr = outFile.string();
        errStrm->open(outFileStr);
        IPPCheck::assertCheck(errStrm->is_open());

        IPhreeqcPhast* iphreeqc = workers[i];
        iphreeqc->Set_error_ostream(errStrm);
    }
}

void PhreeqcRMHelper::findFeatures(const size_t nCells, const CellIdToPhreeqcFeatureMask& featureMaskPerDomain,
                                   const std::vector<int>& cells, const size_t iDomain, std::vector<int>& ic)
{
    MPI_CHECK_SYNC;

    const unsigned int featureMask = featureMaskPerDomain.at(iDomain);

    for(size_t iFeature = 0; iFeature < s_allFeatures.size(); ++iFeature)
    {
        const PhreeqcFeatures feature = s_allFeatures[iFeature];
        if(featureMask & feature)
        {
            const size_t offset = iFeature * nCells;

            for(int iCell : cells)
            {
                // std::cout << "FT:\t" << iDomain << " / " << iCell << "\t->\t" << iFeature << std::endl;
                ic.at(offset + iCell) = iDomain;
            }
        }
    }

    MPI_CHECK_SYNC
}

void PhreeqcRMHelper::printNegativePorosWarnings(const std::vector<double>& porosities)
{
    const double tolerance = 0.001;
    std::vector<size_t> majorNegPorosCells;
    std::vector<size_t> minorNegPorosCells;
    for(size_t i = 0; i < porosities.size(); ++i)
    {
        const double& porosity = porosities[i];
        if(porosity < 0.0)
        {
            if(std::abs(porosity) > tolerance)
            {
                majorNegPorosCells.push_back(i);
            }
            else
            {
                minorNegPorosCells.push_back(i);
            }
        }
    }

    if(majorNegPorosCells.empty() == false)
    {
        std::stringstream ss;
        ss << "got major negative porosity value (tolerance: " << tolerance << ") in cells: ";
        for(size_t i = 0; i < majorNegPorosCells.size(); ++i)
        {
            const size_t index = majorNegPorosCells[i];
            ss << index << ", ";
        }

        throw std::runtime_error(ss.str());
    }

    if(minorNegPorosCells.empty() == false)
    {
        std::cout << "got minor negative porosity value (tolerance: " << tolerance << ") in cells: ";
        for(size_t i = 0; i < minorNegPorosCells.size(); ++i)
        {
            const size_t index = minorNegPorosCells[i];
            std::cout << index << ", ";
        }
        std::cout << std::endl;
    }
}

void PhreeqcRMHelper::findMappingAndFeatures(const ReactionModuleConfig &conf,
                                             const PhreeqcGlobalInitData &initData,
                                             const size_t nxyz, std::vector<int> &grid2chem,
                                             std::vector<int> &ic, std::vector<size_t> &cellToDomain)
{
    MPI_CHECK_SYNC;

    const CellIndexCorrection& indexCorr = conf.indexCorr;

    // setup mapping due to MPI decomposition
    grid2chem.resize(nxyz, -1);
    for(size_t i = 0; i < grid2chem.size(); ++i)
    {
        const size_t corrIndex = indexCorr.correctForDecomp(i);
        grid2chem[i] = corrIndex;
    }
    assert(std::find(grid2chem.begin(), grid2chem.end(), -1) == grid2chem.end());



    std::vector<int> missingCells(nxyz);
    std::iota(missingCells.begin(), missingCells.end(), 0);

    static const int definedCellTag = -1;
    static const int notInitializeTag = -1;

    ic.resize(nxyz * s_allFeatures.size(), notInitializeTag);
    cellToDomain.resize(nxyz, -1);

    const std::vector<Domain>& domains = conf.ippConfig.domains;

    size_t sumCellsDefined = 0;

    for(size_t iDomain = 0 ; iDomain < domains.size(); ++iDomain)
    {
        const Domain& domain = domains[iDomain];
        const std::vector<int>& cells = domain.cells;

        sumCellsDefined += cells.size();
        if(sumCellsDefined > nxyz)
        {
            throw std::runtime_error("more domain cells are defined than space in system");
        }

        pcout << "init domain x nCells: " << iDomain
              << " x " << cells.size() << std::endl;

        std::vector<int> correctedCells(cells.size(), -1);
        for(size_t i = 0; i < cells.size(); ++i)
        {
            const int oldIndex = cells[i];
            const size_t correctedIndex = indexCorr.correctForBCs(oldIndex);
            correctedCells[i] = correctedIndex;

            assert(missingCells.size() > (size_t)correctedIndex);
            missingCells[correctedIndex] = definedCellTag;

            assert(cellToDomain.size() > correctedIndex);
            assert(cellToDomain[correctedIndex] == (size_t)-1);
            cellToDomain[correctedIndex] = iDomain;
        }

        PhreeqcRMHelper::findFeatures(nxyz, initData.featureMask,
                                      correctedCells, iDomain, ic);


    }

    missingCells.erase(std::remove(missingCells.begin(), missingCells.end(),
                                   definedCellTag), missingCells.end());


    // init missing cells to generic cell ID
    if(missingCells.empty() == false)
    {
        // FIXME: gives false positive warnings for boundary conditions
        pcout << "WARNING: number of undefined cells (BC cells included): "
                  << missingCells.size() << std::endl;
    }

    PhreeqcRMHelper::findFeatures(nxyz, initData.featureMask,
                                  missingCells, initData.genericCellID, ic);

    MPI_CHECK_SYNC
}



} // end of namespace IPP
