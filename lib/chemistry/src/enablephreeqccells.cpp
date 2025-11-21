#include "enablephreeqccells.h"

#include <cassert>

#include "localaccessphreeqcrm.h"

namespace IPP
{

EnablePhreeqcCells::EnablePhreeqcCells(LocalAccessPhreeqcRM& phreeqc,
                                       std::vector<char>& enabledCells)
    : m_phreeqc(phreeqc)
    , m_enabledCells(enabledCells)
{

}

EnablePhreeqcCells::~EnablePhreeqcCells()
{

}


void EnablePhreeqcCells::enable(const std::vector<size_t>& toEnable)
{
    for(size_t i = 0; i < toEnable.size(); ++i)
    {
        const size_t iCell = toEnable[i];
        this->enableCell(iCell);
    }
}

void EnablePhreeqcCells::enable(const size_t toEnable)
{
    this->enableCell(toEnable);
}

bool EnablePhreeqcCells::needTotals() const
{
    return false;
}


void EnablePhreeqcCells::enableCell(const size_t iCellLocal)
{
    assert(m_enabledCells.at(iCellLocal) == false);
    m_enabledCells[iCellLocal] = true;
    m_phreeqc.enableCell(iCellLocal);

}

}
