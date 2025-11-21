#ifndef DISSOLVEONLYONBELOWTHRESHOLD_H
#define DISSOLVEONLYONBELOWTHRESHOLD_H

#include "dissolveonlybase.h"


#include "dissprecipbehaviour.h"
#include "porositythresholds.h"

namespace IPP
{

class BelowThreshDissolveOnly : public DissolveOnlyBase
{
public:
    BelowThreshDissolveOnly(const double& porosLowerThresh);

    virtual double getLowerPorosityThresh() const override;

    virtual bool needsNeighInfos() const override;
    virtual bool needsSaturationIndices() const override;
    virtual bool isNucleation() const override;

    virtual void setPrecipDissolveOnlyBehaviour(const std::vector<DissPrecipBehaviour>& behav) override;

private:
    virtual DissolveOnlyBase::PreventPrecipResult preventPrecip(const DissolveOnlyInfos &infos,
                                                                const size_t index) const override;

    const double m_porosLowerThresh;
    std::vector<DissPrecipBehaviour> m_behav;
};

class BelowThreshNeighborAwareDissolveOnly : public DissolveOnlyBase
{
public:
    BelowThreshNeighborAwareDissolveOnly(const PorosityThresholds& thresholds);

    virtual double getLowerPorosityThresh() const override;

    virtual bool needsNeighInfos() const override;
    virtual bool needsSaturationIndices() const override;
    virtual bool isNucleation() const override;

    virtual void setPrecipDissolveOnlyBehaviour(const std::vector<DissPrecipBehaviour>& behav) override;

private:
    virtual PreventPrecipResult preventPrecip(const DissolveOnlyInfos &infos,
                                              const size_t index) const override;


    const PorosityThresholds m_thresholds;
    std::vector<DissPrecipBehaviour> m_behav;

};


}



#endif // DISSOLVEONLYONBELOWTHRESHOLD_H
