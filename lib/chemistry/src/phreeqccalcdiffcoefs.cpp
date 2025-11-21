#include "phreeqccalcdiffcoefs.h"

#include <sstream>

#include "abstractmultiscalediffusioncalc.h"
#include "mpimanager.h"

namespace IPP
{

namespace PhreeqcCalcDiffCoefs
{


void calc(const AbstractMultiScaleDiffusionCalc& diffCalc,
          const std::vector<SimplePorosityInfosPtr>& porosInfos,
          std::vector<double>& diffCoefPerDomain)
{
    const size_t nDomains = porosInfos.size();
    diffCoefPerDomain.resize(nDomains);

    assert(porosInfos.size() == nDomains);

    for(size_t iDomain = 0; iDomain < nDomains; ++iDomain)
    {
        const SimplePorosityInfosPtr& porosInfo = porosInfos[iDomain];

        double& diffCoef = diffCoefPerDomain[iDomain];
        diffCoef = diffCalc.calc(*porosInfo);

        if(diffCoef < 0.0)
        {
            std::stringstream ss;
            ss << MPIManager::getInstance().getRank()
               << ": ERROR: negative diff coef at: " << iDomain
               << " diffCoef: " << diffCoef
               << std::endl;
            throw std::runtime_error(ss.str());
        }


    }
}


}

}
