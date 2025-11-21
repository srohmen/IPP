#include "runcellscalcporosity.h"

#include <PhreeqcRM.h>

#include "runphreeqccells.h"
#include "retrievephaseamounts.h"
#include "porositycalcutils.h"

#include "geometrytools.h"

namespace IPP
{

RunCellsCalcPorosity::RunCellsCalcPorosity(PhreeqcRM& phreeqc,
                                           const RunPhreeqcCells& runCells,
                                           const RetrievePhaseAmounts& retrievePhases,
                                           const PorosityCalc& porosCalc)
    : m_phreeqc(phreeqc),
      m_runCells(runCells),
      m_retrievePhases(retrievePhases),
      m_porosCalc(porosCalc)
{

}

RunCellsCalcPorosity::~RunCellsCalcPorosity()
{

}

void RunCellsCalcPorosity::start()
{
}

void RunCellsCalcPorosity::stop()
{
}

void RunCellsCalcPorosity::operator()(const std::vector<double>& pressure, std::vector<double>& porosities) const
{
    this->setPressureCalcPoros(pressure, porosities);
}

void RunCellsCalcPorosity::setPressureCalcPoros(const std::vector<double>& pressure, std::vector<double>& porosities) const
{
    m_phreeqc.SetPressure(pressure);
    m_runCells.run();

    calcCurrPoros(porosities);
}

void RunCellsCalcPorosity::calcCurrPoros(std::vector<double>& porosities) const
{
    std::vector<double> precipMol;
    m_retrievePhases.retrieve(precipMol);

    const size_t nDomains = m_phreeqc.GetGridCellCount();
    assert(precipMol.size() % nDomains == 0);
    const size_t nPhases = precipMol.size() / nDomains;
    GeometryTools::transpose(precipMol, nDomains, nPhases);

    const av::array_view<const double, 2> precipMolView(precipMol, { (ptrdiff_t)nDomains, (ptrdiff_t)nPhases });

    throw std::runtime_error("FIXME: this is unused / dead code. if still needed adapt to non permeable porosity feature");
    // m_porosCalc.calcPorosity(precipMolView, porosities);
}

void RunCellsCalcPorosity::setEnabledCells(const std::vector<size_t>& enabledCells) const
{

}


}
