#ifndef DISSOLVEONLYBASE_H
#define DISSOLVEONLYBASE_H

#include "abstractdissolveonlycalc.h"


namespace IPP
{

class DissolveOnlyBase : public AbstractDissolveOnlyCalc
{
public:
    DissolveOnlyBase();

    virtual void findNucleationPhases(std::vector<std::string>& /*nucleationPhases*/,
                                      std::vector<std::string>& /*monomers*/) const override
    {
        // do nothing
    }

    virtual void findNucleationMonomers(std::vector<std::string>& /*nucleationMonomers*/) const override
    {
        // do nothing
    }

    virtual void init(const PhaseNameToInfos &phaseInfos,
                      const double& cellVol,
                      const double& cellArea,
                      const double& dt) override;

    virtual bool evaluate(const DissolveOnlyInfos& infos,
                          const iterator& begin, const iterator& end) const override;

protected:
    enum PreventPrecipResult
    {
        PPR_PreventPrecip,
        PPR_AllowPrecip,
        PPR_AllPreventPrecip,
        PPR_AllAllowPrecip
    };

private:
    virtual PreventPrecipResult preventPrecip(const DissolveOnlyInfos &infos,
                                              const size_t index) const = 0;

};

}

#endif // DISSOLVEONLYBASE_H
