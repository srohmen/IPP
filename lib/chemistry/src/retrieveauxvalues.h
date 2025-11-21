#ifndef RETRIEVEAUXVALUES_H
#define RETRIEVEAUXVALUES_H

#include "auxdata.h"

namespace IPP
{
class LocalAccessPhreeqcRM;

class RetrieveAuxValues
{
public:
    RetrieveAuxValues(LocalAccessPhreeqcRM& phreeqc,
                      const std::vector<char>& enabledCells);

    void retrieve(AuxDataVec& auxDataVec) const;

private:
    LocalAccessPhreeqcRM& m_phreeqc;
    const std::vector<char>& m_enabledCells;
};


}

#endif // RETRIEVEAUXVALUES_H
