#include "disablephreeqccells.h"

#include "localaccessphreeqcrm.h"

#include "indexhelper.h"

namespace IPP
{

DisablePhreeqcCells::DisablePhreeqcCells(LocalAccessPhreeqcRM& phreeqc,
                                         std::vector<char>& enabledCells)
    : m_phreeqc(phreeqc)
    , m_enabledCells(enabledCells)
{

}

DisablePhreeqcCells::~DisablePhreeqcCells()
{

}


void DisablePhreeqcCells::disable(const std::vector<size_t>& toDisable)
{
    for(size_t i = 0; i < toDisable.size(); ++i)
    {
        const size_t iCell = toDisable[i];
        this->disable(iCell);
    }

}

void DisablePhreeqcCells::disable(const size_t toDisable)
{
    this->disableCell(toDisable);
}

void DisablePhreeqcCells::disableCell(const size_t iCellLocal)
{
    assert(m_enabledCells.at(iCellLocal));
    m_enabledCells[iCellLocal] = false;
    m_phreeqc.disableCell(iCellLocal);
}


}
