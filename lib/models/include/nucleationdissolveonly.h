#ifndef NUCLEATIONDISSOLVEONLY_H
#define NUCLEATIONDISSOLVEONLY_H

#include "dissolveonlybase.h"

#include <random>
#include <map>
#include <memory>

#include "porositythresholds.h"

namespace IPP
{

class NucleationProp_HON_HEN;

struct NucleationDissolveOnlyData
{
    std::string monomer;
    double V_molecule;
    double D;
    double T;
    double HON_IT_water;
    double HEN_IT_substrate;
    double HEN_IT_self;
    double HON_N0_water;
    double HEN_N0_substrate;
    double HEN_N0_self;
};


using NucleationData = std::map<std::string, NucleationDissolveOnlyData>;

template<typename _PCS_func>
class NucleationDissolveOnlyBase : public DissolveOnlyBase
{
public:
    using PCS_func = _PCS_func;

    NucleationDissolveOnlyBase(const PorosityThresholds& thresholds, const PCS_func& func);
    virtual ~NucleationDissolveOnlyBase();


    NucleationData& getNucleationData();
    PCS_func& getPCSfunc();

    virtual void init(const PhaseNameToInfos &phaseInfos,
                      const double& cellVol,
                      const double& cellArea,
                      const double& dt) override;

    virtual double getLowerPorosityThresh() const override;

    virtual bool needsNeighInfos() const override;
    virtual bool needsSaturationIndices() const override;
    virtual bool isNucleation() const override;

    virtual void findNucleationPhases(std::vector<std::string>& nucleationPhases,
                                      std::vector<std::string>& monomers) const override;
    virtual void findNucleationMonomers(std::vector<std::string>& nucleationMonomers) const override;

    virtual void setPrecipDissolveOnlyBehaviour(const std::vector<DissPrecipBehaviour>& behav) override;


    virtual PreventPrecipResult preventPrecip(const DissolveOnlyInfos &infos,
                                              const size_t phaseIndex) const override;
private:    
    void clear();

    bool isNucleating(const DissolveOnlyInfos &infos, const size_t iPhase) const;

    const PorosityThresholds m_thresholds;
    PCS_func m_pcsFunc;

    std::vector<DissPrecipBehaviour> m_behav;

    NucleationData m_rawData;

    std::vector<NucleationProp_HON_HEN*> m_nuclPropCalc;
    double m_cellVol;
    double m_cellArea;


    std::mt19937 m_generator;
    using Dist = std::uniform_real_distribution<double>;
    Dist m_distribution;
};

class NoPCS
{
public:
    static inline void init(const PhaseNameToInfos &phaseInfos, const NucleationData& nuclData);
    static inline bool isPCS_effective(const DissolveOnlyInfos &infos,
                                       const size_t phaseIndex,
                                       const double& porosLowerThresh);
};

class CylinderPCS
{
public:
    CylinderPCS(const std::string& masterPhase, const double &porosInit, const double& rInit);

    inline void init(const PhaseNameToInfos &phaseInfos, const NucleationData& nuclData);

    inline bool isPCS_effective(const DissolveOnlyInfos &infos,
                                const size_t phaseIndex,
                                const double& porosLowerThresh) const;

private:
    const std::string m_masterPhase;
    size_t m_index;
    const double m_cylLength;
    double k;
};

typedef NucleationDissolveOnlyBase<CylinderPCS> CylinderPCS_NucleationDissolveOnly;

class SpherePCS
{
    // TODO: this class is not implemented completely yet
public:
    SpherePCS(const std::string& masterPhase, const double &porosInit, const double& rInit);

    inline void init(const PhaseNameToInfos &phaseInfos, const NucleationData& nuclData);

    inline bool isPCS_effective(const DissolveOnlyInfos &infos,
                                const size_t phaseIndex,
                                const double& porosLowerThresh) const;

private:
    const std::string m_masterPhase;
};

class NucleationDissolveOnly : public NucleationDissolveOnlyBase<NoPCS>
{
public:
    NucleationDissolveOnly(const PorosityThresholds& thresholds);
};

typedef NucleationDissolveOnlyBase<SpherePCS> SpherePCS_NucleationDissolveOnly;


}

#endif // NUCLEATIONDISSOLVEONLY_H
