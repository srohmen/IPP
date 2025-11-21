#include "cellidpassingdissolveonlyfunc.h"

#include "dissolveonlyinfos.h"
#include <cassert>
#include <iostream>

namespace IPP
{

CellIDPassingDissolveOnlyFunc::CellIDPassingDissolveOnlyFunc(AbstractDissolveOnlyCalcPtr dissOnlyFunc,
                                                             const std::vector<double>& porosity,
                                                             const std::vector<CellNeighborInfo> &neighInfos,
                                                             const std::vector<double>& currentSI,
                                                             const std::vector<size_t>& phaseToNuclPhaseID,
                                                             const std::vector<size_t>& nuclPhaseToMonomerID,
                                                             const std::vector<double>& monomerConc,
                                                             const PhaseNameToInfos &phaseInfos)
    : m_dissOnlyFunc(dissOnlyFunc)
    , m_porosity(porosity)
    , m_neighInfos(neighInfos)
    , m_currentSI(currentSI)
    , m_phaseToNuclPhaseID(phaseToNuclPhaseID)
    , m_nuclPhaseToMonomerID(nuclPhaseToMonomerID)
    , m_monomerConc(monomerConc)
    , m_phaseInfos(phaseInfos)
{

}

void CellIDPassingDissolveOnlyFunc::init(const double& cellVol,
                                         const double& cellArea,
                                         const double& dt)
{
    m_dissOnlyFunc->init(m_phaseInfos, cellVol, cellArea, dt);
    m_nPhases = m_phaseInfos.size();

    const size_t nCells = m_porosity.size();
    assert(m_monomerConc.size() % nCells == 0);
    m_nMonomers = m_monomerConc.size() / nCells;
}

double CellIDPassingDissolveOnlyFunc::getLowerPorosityThresh() const
{
    return m_dissOnlyFunc->getLowerPorosityThresh();
}

bool CellIDPassingDissolveOnlyFunc::needsNeighInfos() const
{
    return m_dissOnlyFunc->needsNeighInfos();
}

bool CellIDPassingDissolveOnlyFunc::needsSaturationIndices() const
{
    return m_dissOnlyFunc->needsSaturationIndices();
}

bool CellIDPassingDissolveOnlyFunc::evaluate(const size_t iCell,
                                             const iterator& begin,
                                             const iterator& end) const
{
    // std::cout << "checking: " << iCell << std::endl;

    const size_t phaseIndex = iCell * m_nPhases;
    const size_t monomerIndex = iCell * m_nMonomers;

    const DissolveOnlyInfos infos = { m_porosity[iCell],
                                      m_neighInfos[iCell],
                                      m_currentSI.begin() + phaseIndex,
                                      m_phaseToNuclPhaseID,
                                      m_nuclPhaseToMonomerID,
                                      m_monomerConc.begin() + monomerIndex
                                    };

    return m_dissOnlyFunc->evaluate(infos, begin, end);
}


}
