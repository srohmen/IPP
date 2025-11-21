#include "setupphreeqcverbosity.h"

#include <vector>

#include "localaccessphreeqcrm.h"

#include "ippexception.h"
#include "mpitools.h"
#include "mpidelegate.h"

namespace IPP
{

SetupPhreeqcVerbosity::SetupPhreeqcVerbosity(LocalAccessPhreeqcRM &phreeqc,
                                             const bool isPrintChemistryOn)
    : m_phreeqc(phreeqc)
    , m_isPrintChemistryOn(isPrintChemistryOn)
{

}

void SetupPhreeqcVerbosity::runLocal() const
{
    // MPI_CHECK_SYNC;


    // Set cells to print chemistry when print chemistry is turned on
    const size_t nCells = m_phreeqc.getNumberCells();

    if(m_isPrintChemistryOn)
    {
        m_phreeqc.SetScreenOn(true);

        std::vector<int> print_chemistry_mask(nCells, 1);
        m_phreeqc.setPrintChemistryMaskLocal(print_chemistry_mask);

        m_phreeqc.setPrintChemistryOnLocal(true, true, true); // workers, initial_phreeqc, utility
    }
    else
    {
        m_phreeqc.SetScreenOn(false);

        std::vector<int> print_chemistry_mask(nCells, 0);
        m_phreeqc.setPrintChemistryMaskLocal(print_chemistry_mask);

        m_phreeqc.setPrintChemistryOnLocal(false, false, false); // workers, initial_phreeqc, utility
    }


}

void SetupPhreeqcVerbosity::run() const
{
    assert(m_phreeqc.GetMpiMyself() == 0);

    // Set cells to print chemistry when print chemistry is turned on

    IRM_RESULT status = IRM_FAIL;

    if(m_isPrintChemistryOn)
    {
        status = m_phreeqc.SetScreenOn(true);
        IPPCheck::assertCheck(status == IRM_OK);

        std::vector<int> print_chemistry_mask(m_phreeqc.GetChemistryCellCount(), 1);
        status = m_phreeqc.SetPrintChemistryMask(print_chemistry_mask);
        IPPCheck::assertCheck(status == IRM_OK);

        status = m_phreeqc.SetPrintChemistryOn(true, true, true); // workers, initial_phreeqc, utility
        IPPCheck::assertCheck(status == IRM_OK);
    }
    else
    {
        status = m_phreeqc.SetScreenOn(false);
        IPPCheck::assertCheck(status == IRM_OK);

        std::vector<int> print_chemistry_mask(m_phreeqc.GetChemistryCellCount(), 0);
        status = m_phreeqc.SetPrintChemistryMask(print_chemistry_mask);
        IPPCheck::assertCheck(status == IRM_OK);

        status = m_phreeqc.SetPrintChemistryOn(false, false, false); // workers, initial_phreeqc, utility
        IPPCheck::assertCheck(status == IRM_OK);
    }

}

}
