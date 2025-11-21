#include "runphreeqccells.h"

#include "localaccessphreeqcrm.h"

#include "ippexception.h"
#include "phreeqcdump.h"
#include "phreeqcinputfilegenerator.h"
#include "phreeqcdefaultconvergencevalues.h"
#include "phreeqcsetfileprefix.h"
#include "setupphreeqcverbosity.h"
#include "mpitools.h"
#include "mpidelegate.h"

namespace IPP
{

RunPhreeqcCells::RunPhreeqcCells(LocalAccessPhreeqcRM &phreeqc,
                                 const PhreeqcSetFilePrefix &setPrefix,
                                 const SetupPhreeqcVerbosity &setupVerbosity)
    : m_phreeqc(phreeqc)
    , m_setPrefix(setPrefix)
    , m_setupVerbosity(setupVerbosity)
{

}

void RunPhreeqcCells::run() const
{
    MPI_CHECK_SYNC;

    m_setPrefix.runLocal();
    IRM_RESULT status = m_phreeqc.runCellsIfEnabled();
    this->checkConvergence(status);

    MPI_CHECK_SYNC;
}

void RunPhreeqcCells::checkConvergence(IRM_RESULT status) const
{
    MPI_CHECK_SYNC;

    if(status != IRM_OK)
    {
        std::cout << "WARNING: RunCells gone wrong. Try to recover with changing convergence parameters" << std::endl;

        // enable debug output
        m_phreeqc.SetScreenOn(true);
        std::vector<int> enablePrint(m_phreeqc.getNumberCells(), 1);
        m_phreeqc.setPrintChemistryMaskLocal(enablePrint);
        m_phreeqc.setPrintChemistryOnLocal(true, true, true);
        m_phreeqc.runCellsLocal();

        size_t iterations = PhreeqcDefaultConvergenceValues::iterations;
        double tolerance = PhreeqcDefaultConvergenceValues::tolerance;
        double convergenceTolerance = PhreeqcDefaultConvergenceValues::convergence_tolerance;
        bool enableDiagonalScale = PhreeqcDefaultConvergenceValues::diagonalScaling;

        if(enableDiagonalScale == false)
        {
            enableDiagonalScale = true;

            std::cout << "WARNING: Try to recover with enabled diagonal scale" << std::endl;
            std::stringstream enableSS;
            PhreeqcInputFileGenerator::generateConvergenceHacks(iterations, tolerance,
                                                                convergenceTolerance,
                                                                enableDiagonalScale,
                                                                enableSS);
            m_phreeqc.runStringLocal(true, false, false, enableSS);
            status = m_phreeqc.runCellsLocal();
        }

        size_t nTry = 0;
        bool isConverged = false;

        while (status != IRM_OK) {
            if (nTry >= 20) {
                isConverged = false;
                break;
            }

            convergenceTolerance *= 1.5;
            iterations *= 2.0;
            tolerance *= 0.5;

            std::ios::fmtflags f(std::cout.flags());

            std::cout << "WARNING: RunCells gone wrong (again...)\n\t-> Decreasing tolerance to: "
                      << std::scientific << tolerance
                      << "\n\t-> Increasing convergenceTolerance to: " << convergenceTolerance
                      << "\n\t-> Increasing iterations to: " << iterations << std::endl;

            std::cout.flags(f);

            std::stringstream enableSS;
            PhreeqcInputFileGenerator::generateConvergenceHacks(iterations,
                                                                tolerance,
                                                                convergenceTolerance,
                                                                enableDiagonalScale,
                                                                enableSS);
            m_phreeqc.runStringLocal(true, false, false, enableSS);
            status = m_phreeqc.runCellsLocal();

            if (status != IRM_OK) {
                isConverged = false;
            } else {
                isConverged = true;
            }

            ++nTry;
        }

        if(status == IRM_OK)
	{
            std::cout << "Finally found convergence. Exec final RunCells" << std::endl;
	}
	else	
	{
        std::cout << "Could not find convergence. Nevertheless, exec final RunCells" << std::endl;
    }

    if (status == IRM_OK || isConverged == false) {
        m_setupVerbosity.runLocal();
        status = m_phreeqc.runCellsLocal();
        std::cout << "Final RunCells status: " << status << std::endl;
        m_phreeqc.DecodeError(status);

        std::cout << "Resetting to default settings" << std::endl;
        // resetting to default
        std::stringstream disableSS;
        PhreeqcInputFileGenerator::generateConvergenceHacks(
            PhreeqcDefaultConvergenceValues::iterations,
            PhreeqcDefaultConvergenceValues::tolerance,
            PhreeqcDefaultConvergenceValues::convergence_tolerance,
            PhreeqcDefaultConvergenceValues::diagonalScaling,
            disableSS);
        m_phreeqc.runStringLocal(true, false, false, disableSS);
    } else {
        // throw std::runtime_error("Could not find phreeqc convergence");
    }
    }

    // IPPCheck::assertCheck(status == IRM_OK, "RunCells gone wrong");

    MPI_CHECK_SYNC;
}


}
