#ifndef PHREEQCCALCPOROSITY_H
#define PHREEQCCALCPOROSITY_H

#include "phasenametoinfos.h"
#include "simpleporosityinfosfwd.h"

namespace IPP
{

class AbstractPorosityCalc;


namespace PhreeqcCalcPorosity
{

void calc(const PhaseNameToInfos& phaseInfos,
          const std::vector<double>& solidFractions,
          const std::vector<double>& precipMol,
          const std::vector<double>& inertVolFrac,
          const AbstractPorosityCalc& porosCalc,
          std::vector<SimplePorosityInfosPtr>& porosInfos);

}


}

#endif // PHREEQCCALCPOROSITY_H
