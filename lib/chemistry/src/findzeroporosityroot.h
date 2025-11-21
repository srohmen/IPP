#ifndef FINDZEROPOROSITYROOT_H
#define FINDZEROPOROSITYROOT_H

#include <vector>


namespace IPP
{

class LocalAccessPhreeqcRM;
class RetrievePhaseAmounts;
class PorosityCalc;

namespace FindZeroPorosityRoot
{


class SetTotalsAndRun
{
public:
    SetTotalsAndRun(LocalAccessPhreeqcRM& phreeqc,
                    const std::vector<double>& totalDiff,
                    const double& porosity,
                    const std::size_t iCell);

    void setAndRun(const double& t);

private:
    LocalAccessPhreeqcRM& m_phreeqc;
    const std::vector<double>& m_totalDiff;
    const double m_porosity;
    const std::size_t m_iCell;

    double m_currInterp;
};

struct Result
{
    double t;
    double porosity;
};

Result findRoot(SetTotalsAndRun& setAndRun,
                const RetrievePhaseAmounts& retrievePhases,
                const PorosityCalc& porosCalc,
                const double& oldPoros,
                const double& currPoros,
                const std::size_t iCell);

}

}

#endif // FINDZEROPOROSITYROOT_H
