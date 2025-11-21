#ifndef INITPOROSITY_H
#define INITPOROSITY_H

#include "phasenametoinfos.h"

namespace IPP
{

class AbstractPorosityCalc;
class AbstractMultiScaleDiffusionCalc;
class SimulationExchangeData;
class CellIndexCorrection;
struct Domain;

namespace InitPorosity
{

void extractPhasesMolarVolumes(const PhaseNameToInfos& phaseNameToInfos,
                               const std::vector<std::string>& phaseNames,
                               std::vector<double>& phaseMolarVolumesInLiter);

void run(const std::vector<Domain>& domains,
         const CellIndexCorrection &indexCorr,
         const size_t nCells, const PhaseNameToInfos &phaseInfos,
         const std::vector<double>& phaseMolarVolumes,
         const std::vector<char>& isPermeablePhase,
         const std::vector<std::string>& phaseNames,
         const AbstractPorosityCalc &porosCalc,
         const AbstractMultiScaleDiffusionCalc &diffCalc,
         std::vector<double>& porosities,
         std::vector<double> &capillaryPorosities,
         std::vector<double>& diffCoefs,
         std::vector<double>& inertVolFracs);

void syncData(const std::vector<double>& globalPorosities,
              const std::vector<double> &globalPorositiesCapillary,
              const std::vector<double>& globalDiffCoefs,
              const std::vector<double>& inertVolFrac,
              SimulationExchangeData& simData);

}

}

#endif // INITPOROSITY_H
