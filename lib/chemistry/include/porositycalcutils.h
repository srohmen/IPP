#ifndef POROSITYCALCUTILS_H
#define POROSITYCALCUTILS_H

#include <vector>
#include <cstddef>
#include "array_view.h"

namespace IPP
{

class PorosityCalc
{
public:
    typedef av::array_view<const double, 2> ConstCellPhasesView;
    typedef av::array_view<double, 2> CellPhasesView;

    PorosityCalc(const size_t nDomains,
                 const std::vector<double>& phaseMolarVolumes, // L / mol
                 const std::vector<char>& isPhasePermeable,
                 const std::vector<double> &inertVolFrac);

    size_t getnDomains() const;
    size_t getnPhases() const;

    void calcPorosity(const ConstCellPhasesView& precipMol,
                      CellPhasesView& phasesVolumeFractions,
                      std::vector<double>& porosResult) const;

    void calcPorosity(const ConstCellPhasesView& precipMol,
                      std::vector<double>& porosResult) const;

private:
    void calcPhasesVolumeFractions(const ConstCellPhasesView& precipMol, CellPhasesView& phaseVolumes) const;

    void sumTotalSolidFractions(const ConstCellPhasesView& phaseFractions,
                                std::vector<double>& totalPhasesFractions) const;
    void sumTotalSolidFractions(const ConstCellPhasesView& phaseFractions,
                                std::vector<double>& totalPhasesFractions,
                                std::vector<double>& totalNonPermPhasesFractions) const;


    void calcPorosityTotalFromSolidFrac(const std::vector<double>& totalPhasesVolumes,
                                   std::vector<double>& porosResult) const;

    const size_t nDomains;

    const std::vector<double>& m_phaseMolarVolumes; // L/mol
    const std::vector<char>& m_isPhasePermeable;
    const std::vector<double>& m_inertVolFrac;

};

} // end of namespace LBGeoChem

#endif // POROSITYCALCUTILS_H

