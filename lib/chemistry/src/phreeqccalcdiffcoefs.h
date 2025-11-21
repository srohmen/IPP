#ifndef PHREEQCCALCDIFFCOEFS_H
#define PHREEQCCALCDIFFCOEFS_H

#include "phasenametoinfos.h"
#include "simpleporosityinfosfwd.h"

namespace IPP
{

class AbstractMultiScaleDiffusionCalc;


namespace PhreeqcCalcDiffCoefs
{

void calc(const AbstractMultiScaleDiffusionCalc& diffCalc,
          const std::vector<SimplePorosityInfosPtr>& porosInfos,
          std::vector<double>& diffCoefPerDomain);

}

}

#endif // PHREEQCCALCDIFFCOEFS_H
