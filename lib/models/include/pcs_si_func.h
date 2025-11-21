#ifndef PCS_SI_FUNC_H
#define PCS_SI_FUNC_H

#include "abstractsicalc.h"
#include <map>
#include "porositycontrolledsolubility.h"

namespace IPP
{

template<typename PCSFunc>
class PCS_SI_Func : public AbstractSICalc
{
public:
    PCS_SI_Func(const double &temp, const double& porosInit, const double& rInit,
                const double& surfTensFallbackConc,
                const std::map<std::string, double>& surfaceTensions);

    virtual void init(const std::vector<std::string> &phaseNames,
                      const PhaseNameToInfos &phaseInfos,
                      const std::vector<double>* poros) override;

    virtual void evaluate(const size_t iCell, const iterator& begin, const iterator& end) const override;



private:
    const double temp;
    const double porosInit;
    const double rInit;
    const double surfTensFallbackConc;
    const std::map<std::string, double> surfaceTensions;

    std::vector<PCSFunc> m_pcsFuncs;
    const std::vector<double>* m_poros;
};

using PCS_FuncCylinder = PCS_SI_Func<PorosityControlledSolubility::Cylinder>;

}

#endif // PCS_SI_FUNC_H
