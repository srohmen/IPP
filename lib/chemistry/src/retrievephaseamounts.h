#ifndef RETRIEVEPHASEAMOUNTS_H
#define RETRIEVEPHASEAMOUNTS_H

#include <vector>



namespace IPP
{
class LocalAccessPhreeqcRM;

class RetrievePhaseAmounts
{
public:
    RetrievePhaseAmounts(LocalAccessPhreeqcRM& phreeqc, const std::vector<char>& enabledCells);

    void retrieve(std::vector<double> &precipMoles) const;

private:
    LocalAccessPhreeqcRM& m_phreeqc;
    const std::vector<char>& m_enabledCells;
};

}
#endif // RETRIEVEPHASEAMOUNTS_H
