#include "porositycalcutils.h"

#include <functional>
#include <algorithm>
#include <assert.h>
#include <iostream>

#include "ippexception.h"

namespace IPP
{


PorosityCalc::PorosityCalc(const size_t nDomains,
                           const std::vector<double>& phaseMolarVolumes,
                           const std::vector<char> &isPhasePermeable,
                           const std::vector<double> &inertVolFrac)
    : nDomains(nDomains),
      m_phaseMolarVolumes(phaseMolarVolumes),
      m_isPhasePermeable(isPhasePermeable),
      m_inertVolFrac(inertVolFrac)
{

}

size_t PorosityCalc::getnDomains() const
{
    return nDomains;
}

size_t PorosityCalc::getnPhases() const
{
    return m_phaseMolarVolumes.size();
}

void PorosityCalc::calcPorosity(const ConstCellPhasesView& precipMol,
                                CellPhasesView& phasesVolumeFractions,
                                std::vector<double>& porosResult) const
{
    this->calcPhasesVolumeFractions(precipMol, phasesVolumeFractions);

    std::vector<double> totalSolidFractions(nDomains, 0.0);
    std::vector<double> nonPermSolidFractions(nDomains, 0.0);
    this->sumTotalSolidFractions(phasesVolumeFractions,
                                 totalSolidFractions,
                                 nonPermSolidFractions);

    calcPorosityTotalFromSolidFrac(totalSolidFractions, porosResult);
}

void PorosityCalc::calcPorosity(const PorosityCalc::ConstCellPhasesView& precipMol,
                                std::vector<double>& porosResult) const
{
    const size_t nPhases = m_phaseMolarVolumes.size();

    av::bounds<2> bounds = { (ptrdiff_t)nDomains, (ptrdiff_t)nPhases };
    std::vector<double> phaseVolumesRaw(nPhases * nDomains, 0.0);
    CellPhasesView phasesVolumeFractions(phaseVolumesRaw, bounds);
    this->calcPhasesVolumeFractions(precipMol, phasesVolumeFractions);

    std::vector<double> totalSolidFractions(nDomains, 0.0);
    std::vector<double> nonPermSolidFractions(nDomains, 0.0);
    this->sumTotalSolidFractions(phasesVolumeFractions,
                                 totalSolidFractions,
                                 nonPermSolidFractions);

    this->calcPorosityTotalFromSolidFrac(totalSolidFractions, porosResult);

}

void PorosityCalc::calcPorosityTotalFromSolidFrac(const std::vector<double>& totalSolidFrac,
                                                  std::vector<double>& porosResult) const
{
    porosResult.resize(totalSolidFrac.size());
    for(size_t iDomain = 0; iDomain < totalSolidFrac.size(); ++iDomain)
    {
        const double solidFrac = totalSolidFrac[iDomain];
        const double inertFrac = m_inertVolFrac[iDomain];
        const double porosity = 1.0 - (solidFrac + inertFrac);
        porosResult[iDomain] = porosity;
    }
}

void PorosityCalc::sumTotalSolidFractions(const ConstCellPhasesView &phaseFractions,
                                          std::vector<double>& totalPhasesFractions) const
{
    const size_t nDomains = phaseFractions.bounds()[0];
    const size_t nPhases = phaseFractions.bounds()[1];

    assert(m_isPhasePermeable.size() == nPhases);

    totalPhasesFractions.resize(nDomains, 0.0);

    for(size_t iDomain = 0; iDomain < nDomains; ++iDomain)
    {
        double& totalFrac = totalPhasesFractions[iDomain];

        for(size_t iPhase = 0; iPhase < nPhases; ++iPhase)
        {
            const double& frac = phaseFractions[iDomain][iPhase];
            totalFrac += frac;
        }
    }
}

void PorosityCalc::sumTotalSolidFractions(const ConstCellPhasesView &phaseFractions,
                                          std::vector<double>& totalPhasesFractions,
                                          std::vector<double>& totalNonPermPhasesFractions) const
{
    const size_t nDomains = phaseFractions.bounds()[0];
    const size_t nPhases = phaseFractions.bounds()[1];

    assert(m_isPhasePermeable.size() == nPhases);

    totalPhasesFractions.resize(nDomains, 0.0);
    totalNonPermPhasesFractions.resize(nDomains, 0.0);

    for(size_t iDomain = 0; iDomain < nDomains; ++iDomain)
    {
        double totalFrac = 0.0;
        double totalFracNonPerm = 0.0;

        for(size_t iPhase = 0; iPhase < nPhases; ++iPhase)
        {
            const double& frac = phaseFractions[iDomain][iPhase];
            totalFrac += frac;

            if(m_isPhasePermeable[iPhase] == false)
            {
                totalFracNonPerm += frac;
            }
        }

        totalPhasesFractions[iDomain] = totalFrac;
        totalNonPermPhasesFractions[iDomain] = totalFracNonPerm;
    }
}

void PorosityCalc::calcPhasesVolumeFractions(const ConstCellPhasesView &precipMol,
                                     CellPhasesView &phaseVolumes) const
{
    for(size_t iDomain = 0; iDomain < nDomains; ++iDomain)
    {
        for(size_t iPhase = 0; iPhase < m_phaseMolarVolumes.size(); ++iPhase)
        {
            const double& molVol = m_phaseMolarVolumes[iPhase];

            const double& amount = precipMol[iDomain][iPhase];
            const double volume = amount * molVol; // L
            phaseVolumes[iDomain][iPhase] = volume; // V in liter is already the same as fraction
        }
    }

}


} // end of namespace IPP
