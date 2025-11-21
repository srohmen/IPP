#ifndef ENABLEPHREEQCCELLS_H
#define ENABLEPHREEQCCELLS_H

#include <vector>
#include <cstddef>

namespace IPP
{

class LocalAccessPhreeqcRM;

class EnablePhreeqcCells
{
public:
    EnablePhreeqcCells(LocalAccessPhreeqcRM& phreeqc,
                       std::vector<char>& enabledCells);
    virtual ~EnablePhreeqcCells();

    virtual void enable(const std::vector<size_t>& toEnable);

    virtual void enable(const size_t toEnable);

    virtual bool needTotals() const;

protected:
    void enableCell(const size_t iCellLocal);

private:

    LocalAccessPhreeqcRM& m_phreeqc;
    std::vector<char>& m_enabledCells;

};

}

#endif // ENABLEPHREEQCCELLS_H
