#ifndef POROSITYHELPER_H
#define POROSITYHELPER_H

#include <vector>

namespace IPP
{

struct ReactionExchangeData;
class SimulationExchangeData;


namespace PorosityHelper
{

struct CalcPorosityInput
{
    const std::size_t nDomains;
    const std::size_t nPhases;
    const double invTotalCellVolume;

    const std::vector<double>& phaseMolarVolumes;
    const std::vector<char>& isPermeablePhase;
    const std::vector<double>& inertVolFrac;
    const std::vector<double>& precipField;
};

struct CalcPorosityOutput
{
    std::vector<double>& localPrecipMol;
    std::vector<double>& solidFractions;
    std::vector<double>& porosityTotal;
};

void calcPorosity(const std::vector<double>& localPrecipMol,
                  const std::vector<double>& solidFractions,
                  std::vector<double>& porosityPhysical,
                  std::vector<double>& porositiesNonPerm);

void calc(const ReactionExchangeData &reacData,
          CalcPorosityOutput &output, const bool verbose);




void calcPorosity(const CalcPorosityInput& input,
                  CalcPorosityOutput& output);

}

}
#endif // POROSITYHELPER_H
