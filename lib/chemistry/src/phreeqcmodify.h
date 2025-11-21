#ifndef PHREEQCMODIFY_H
#define PHREEQCMODIFY_H

#include <string>
#include <vector>

class PhreeqcRM;

namespace IPP
{

class CellID_SI_Func;
class CellID_DissolveOnly_Func;

namespace PhreeqcModify
{

void updateTargetSaturationIndices(const CellID_SI_Func& siFunc,
                                   const std::vector<std::string>& phaseNames,
                                   const size_t nCells,
                                   PhreeqcRM& phreeqc,
                                   std::vector<double> &newSImap);

void updateDissolveOnlyFlags(const CellID_DissolveOnly_Func& doFunc,
                             const std::vector<std::string>& phaseNames,
                             const size_t nCells,
                             PhreeqcRM& phreeqc,
                             std::vector<unsigned char> &currDissolveOnlyCells);

}

} // end of namespace LBGeoChem

#endif // PHREEQCMODIFY_H
