#ifndef RUNCELLSCALCPOROSITY_H
#define RUNCELLSCALCPOROSITY_H

#include <vector>
#include <cstddef>

class PhreeqcRM;

namespace IPP
{
class RunPhreeqcCells;
class RetrievePhaseAmounts;
class PorosityCalc;

class RunCellsCalcPorosity
{
public:
    RunCellsCalcPorosity(PhreeqcRM& phreeqc,
                         const RunPhreeqcCells& runCells,
                         const RetrievePhaseAmounts& retrievePhases,
                         const PorosityCalc& porosCalc);

    ~RunCellsCalcPorosity();

    void start();
    void stop();

    void operator()(const std::vector<double>& pressure, std::vector<double>& porosities) const;
    void setPressureCalcPoros(const std::vector<double>& pressure, std::vector<double>& porosities) const;
    void calcCurrPoros(std::vector<double>& porosities) const;
    void setEnabledCells(const std::vector<size_t>& enabledCells) const;

private:
    PhreeqcRM& m_phreeqc;
    const RunPhreeqcCells& m_runCells;
    const RetrievePhaseAmounts& m_retrievePhases;
    const PorosityCalc& m_porosCalc;
    std::vector<double> m_initSaturation;
};

}
#endif // RUNCELLSCALCPOROSITY_H
