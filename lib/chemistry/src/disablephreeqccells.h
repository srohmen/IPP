#ifndef DISABLEPHREEQCCELLS_H
#define DISABLEPHREEQCCELLS_H

#include <vector>
#include <cstddef>

namespace IPP
{

class LocalAccessPhreeqcRM;

class DisablePhreeqcCells
{
public:
    DisablePhreeqcCells(LocalAccessPhreeqcRM& phreeqc,
                        std::vector<char>& enabledCells);

    virtual ~DisablePhreeqcCells();

    void disable(const std::vector<size_t>& toDisable);

    void disable(const size_t toDisable);


private:
    void disableCell(const size_t iCellLocal);

    LocalAccessPhreeqcRM& m_phreeqc;
    std::vector<char>& m_enabledCells;

};

}

#endif // DISABLEPHREEQCCELLS_H
