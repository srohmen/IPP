#include "nucleationdissolveonly.h"

#include <iostream>
#include <chrono>

#include "dissolveonlyutils.h"
#include "nucleationpropcalc.h"
#include "mpimanager.h"
#include "ippexception.h"

namespace IPP
{

class NucleationProp_HON_HEN
{
public:
    NucleationProp_HON_HEN(const NucleationDissolveOnlyData& nuclData,
                           const double& dt);

    double calcLnSurvivalProp(const double& SI,
                              const double& N1,
                              const double& cellVol,
                              const double& substrateSurf,
                              const double& selfSurf,
                              const bool homoOnly) const;

    double calcNuclProp(const double& SI,
                        const double& N1,
                        const double& cellVol,
                        const double& substrateSurf,
                        const double& selfSurf,
                        const bool homoOnly) const;

private:
    NucleationPropCalc m_HON;
    NucleationPropCalc m_HEN_substrate;
    NucleationPropCalc m_HEN_self;
};


NucleationProp_HON_HEN::NucleationProp_HON_HEN(const NucleationDissolveOnlyData& nuclData,
                                               const double& dt)
    : m_HON(dt, nuclData.T, nuclData.D, nuclData.V_molecule,
            nuclData.HON_N0_water, nuclData.HON_IT_water)
    , m_HEN_substrate(dt, nuclData.T, nuclData.D, nuclData.V_molecule,
                      nuclData.HEN_N0_substrate, nuclData.HEN_IT_substrate)
    , m_HEN_self(dt, nuclData.T, nuclData.D, nuclData.V_molecule,
                 nuclData.HEN_N0_self, nuclData.HEN_IT_self)
{

}

double NucleationProp_HON_HEN::calcLnSurvivalProp(const double &SI,
                                                  const double &N1,
                                                  const double &cellVol,
                                                  const double &substrateSurf,
                                                  const double &selfSurf,
                                                  const bool homoOnly) const
{
    const double SR = std::pow(10, SI);
    const double lnS = std::log(SR);

    const double propHON = m_HON.calcLnSurvivalProp(lnS, N1, cellVol);

    if(homoOnly)
    {
        return propHON;
    }
    else
    {
        const double propHEN_substrate = substrateSurf > 0.0 ?
                    m_HEN_substrate.calcLnSurvivalProp(lnS, N1, substrateSurf)
                  :
                    0.0;

        const double propHEN_self = selfSurf > 0.0 ?
                    m_HEN_self.calcLnSurvivalProp(lnS, N1, selfSurf)
                  :
                    0.0;

        const double sum = propHON + propHEN_substrate + propHEN_self;

        return sum;
    }


}

double NucleationProp_HON_HEN::calcNuclProp(const double& SI,
                                            const double& N1,
                                            const double& cellVol,
                                            const double& substrateSurf,
                                            const double& selfSurf,
                                            const bool homoOnly) const
{
    const double SR = std::pow(10, SI);
    const double lnS = std::log(SR);

    const double propHON = m_HON.calcNucleationProp(lnS, N1, cellVol);

    if(homoOnly)
    {
        return propHON;
    }
    else
    {

        const double propHEN_substrate = substrateSurf > 0.0 ?
                    m_HEN_substrate.calcNucleationProp(lnS, N1, substrateSurf)
                  :
                    0.0;

        const double propHEN_self = selfSurf > 0.0 ?
                    m_HEN_self.calcNucleationProp(lnS, N1, selfSurf)
                  :
                    0.0;

        // R = A+B+C - AB - AC - BC + ABC
        const double propCombined = propHON + propHEN_substrate + propHEN_self
                - (propHON * propHEN_substrate)
                - (propHON * propHEN_self)
                - (propHEN_substrate * propHEN_self)
                + (propHON * propHEN_substrate * propHEN_self);

        return propCombined;
    }
}

template<typename PCS_func>
NucleationDissolveOnlyBase<PCS_func>::NucleationDissolveOnlyBase(const PorosityThresholds &thresholds, const PCS_func &func)
    : m_thresholds(thresholds)
    , m_pcsFunc(func)
    , m_cellVol(-1)
    , m_cellArea(-1)
    , m_generator(std::chrono::system_clock::now().time_since_epoch().count())
    , m_distribution(0.0, 1.0)
{

}

template<typename PCS_func>
NucleationDissolveOnlyBase<PCS_func>::~NucleationDissolveOnlyBase()
{
    this->clear();
}

template<typename PCS_func>
NucleationData& NucleationDissolveOnlyBase<PCS_func>::getNucleationData()
{
    return m_rawData;
}

template<typename PCS_func>
PCS_func &NucleationDissolveOnlyBase<PCS_func>::getPCSfunc()
{
    return m_pcsFunc;
}

template<typename PCS_func>
void NucleationDissolveOnlyBase<PCS_func>::init(const PhaseNameToInfos &phaseInfos,
                                                const double& cellVol,
                                                const double& cellArea,
                                                const double& dt)
{
    this->clear();

    m_nuclPropCalc.resize(phaseInfos.size(), nullptr);
    m_cellVol = cellVol;
    m_cellArea = cellArea;

    size_t iPhase = 0;
    for(auto dataIt = phaseInfos.cbegin();
        dataIt != phaseInfos.cend();
        ++iPhase, ++dataIt)
    {
        const std::string& name = dataIt->first;

        auto it = m_rawData.find(name);
        if(it != m_rawData.end())
        {
            const NucleationDissolveOnlyData& data = it->second;
            m_nuclPropCalc[iPhase] = new NucleationProp_HON_HEN(data, dt);
        }
    }


    const size_t time = std::chrono::system_clock::now().time_since_epoch().count();
    const size_t rank = MPIManager::getInstance().getRank();
    const size_t seed = time + rank;
    m_generator.seed(seed);

    m_pcsFunc.init(phaseInfos, m_rawData);
}

template<typename PCS_func>
double NucleationDissolveOnlyBase<PCS_func>::getLowerPorosityThresh() const
{
    return m_thresholds.porosLower;
}

template<typename PCS_func>
bool NucleationDissolveOnlyBase<PCS_func>::needsNeighInfos() const
{
    return true;
}

template<typename PCS_func>
bool NucleationDissolveOnlyBase<PCS_func>::needsSaturationIndices() const
{
    return true;
}

template<typename PCS_func>
bool NucleationDissolveOnlyBase<PCS_func>::isNucleation() const
{
    return true;
}

template<typename PCS_func>
void NucleationDissolveOnlyBase<PCS_func>::findNucleationPhases(std::vector<std::string>& nucleationPhases,
                                                                std::vector<std::string>& monomers) const
{
    assert(nucleationPhases.empty());
    assert(monomers.empty());

    for(const std::pair<const std::string, NucleationDissolveOnlyData>& data : m_rawData)
    {
        assert(std::find(nucleationPhases.begin(), nucleationPhases.end(), data.first)
               == nucleationPhases.end());
        nucleationPhases.push_back(data.first);

        // monomers could be non unique
        monomers.push_back(data.second.monomer);
    }
}

template<typename PCS_func>
void NucleationDissolveOnlyBase<PCS_func>::findNucleationMonomers(std::vector<std::string>& nucleationMonomers) const
{
    assert(nucleationMonomers.empty());

    for(const std::pair<const std::string, NucleationDissolveOnlyData>& data : m_rawData)
    {
        nucleationMonomers.push_back(data.second.monomer);
    }

    std::sort(nucleationMonomers.begin(), nucleationMonomers.end());
    nucleationMonomers.erase(std::unique(nucleationMonomers.begin(), nucleationMonomers.end()),
                             nucleationMonomers.end());
}

template<typename PCS_func>
void NucleationDissolveOnlyBase<PCS_func>::setPrecipDissolveOnlyBehaviour(const std::vector<DissPrecipBehaviour> &behav)
{
    assert(m_behav.empty());
    m_behav = behav;
}

template<typename PCS_func>
DissolveOnlyBase::PreventPrecipResult
NucleationDissolveOnlyBase<PCS_func>::preventPrecip(const DissolveOnlyInfos &infos,
                                                    const size_t phaseIndex) const
{
    const bool isDissolveOnlyBehav = DissolveOnlyUtils::isDissolveOnlyBehav(m_behav, phaseIndex);

    if(isDissolveOnlyBehav)
    {
        return PPR_PreventPrecip;
    }
    else
    {
        const bool isBelowThresh = m_pcsFunc.isPCS_effective(infos, phaseIndex, m_thresholds.porosLower);

        if(isBelowThresh)
        {
            return PPR_AllPreventPrecip;
        }
        else
        {
            const bool isOwnPorosLowerMinumum = infos.porosity < m_thresholds.porosUpper;
            if(isOwnPorosLowerMinumum)
            {
                return PPR_AllowPrecip;
            }
            else
            {
                const bool isNucleating = this->isNucleating(infos, phaseIndex);

                if(isNucleating)
                {
                    return PPR_AllAllowPrecip;
                }
                else
                {
                    return PPR_PreventPrecip;
                }
            }
        }
    }

}

template<typename PCS_func>
void NucleationDissolveOnlyBase<PCS_func>::clear()
{
    for(size_t i = 0; i < m_nuclPropCalc.size(); ++i)
    {
        delete m_nuclPropCalc[i];
        m_nuclPropCalc[i] = nullptr;
    }

    m_nuclPropCalc.clear();
}

template<typename PCS_func>
bool NucleationDissolveOnlyBase<PCS_func>::isNucleating(const DissolveOnlyInfos &infos, const size_t iPhase) const
{    
    const NucleationProp_HON_HEN* nucProp = m_nuclPropCalc[iPhase];

    if(nucProp)
    {
        const DissolveOnlyInfos::ArrIterator satIndexIt = infos.SI_begin + iPhase;
        const double& SI = *satIndexIt;

        if(SI > 0.0)
        {
            assert(infos.phaseToNuclPhaseID.size() > iPhase);
            const size_t iNuclPhase = infos.phaseToNuclPhaseID[iPhase];

            assert(infos.neighInfo.phaseVolFrac.size() > iNuclPhase);
            const CellNeighborInfo::VolumeFractions& fracs = infos.neighInfo.phaseVolFrac[iNuclPhase];
            const double selfSurf = fracs.self * m_cellArea;
            const double substrateSurf = fracs.rest * m_cellArea;

            assert(infos.nuclPhaseToMonomerID.size() > iNuclPhase);
            const size_t iMonomer = infos.nuclPhaseToMonomerID[iNuclPhase];
            const DissolveOnlyInfos::ArrIterator monomerIt = infos.monomerConc_begin + iMonomer;
            const double N_A = 6.022E+23;
            const double& N1 = *monomerIt * N_A;

            const bool homoOnly = infos.neighInfo.location != CellNeighborInfo::L_Surface;
            const double lnSurvivalProp = nucProp->calcLnSurvivalProp(SI, N1, m_cellVol, substrateSurf,
                                                                      selfSurf, homoOnly);

            // save some runtime if survival prop is one
            // p = e^k -> k = ln(p)
            // p >= 1 -> k >= 0
            // nucleation only possible of k is negative
            if(lnSurvivalProp < 0.0)
            {
                // internal state of random number generator changes but this function is only used by one thread per process
                const double randNumTmp = const_cast<Dist&>(m_distribution)(const_cast<std::mt19937&>(m_generator));

                // the uniform distribution produces random numbers in the interval [0, 1) -> zero included and one excluded
                // but we need random numbers which excludes zero (because of later log(0) -> undefined) and including 1
                const double randNum = 1.0 - randNumTmp;
                const double lnRand = std::log(randNum);

                if(lnRand > lnSurvivalProp)
                {
                    // nucleation occured
                    std::cout << "nucleation event for phase " << iPhase << std::endl
                              << "\tHEN included: " << (int)(homoOnly == false) << std::endl
                              << "\tSI: " << SI << std::endl
                              << "\trandNum: " << randNum << std::endl
                              << "\tlnRand: " << lnRand << std::endl
                                 // << "\tnuclProp: " << nuclProp << std::endl;
                              << "\tlnSurvivalProp " << lnSurvivalProp << std::endl;

                    return true;
                }
            }
        }
    }

    return false;
}

NucleationDissolveOnly::NucleationDissolveOnly(const PorosityThresholds &thresholds)
    : NucleationDissolveOnlyBase<NoPCS>(thresholds, NoPCS())
{

}

void NoPCS::init(const PhaseNameToInfos &/*phaseInfos*/, const NucleationData& /*nuclData*/)
{
    // do nothing
}

bool NoPCS::isPCS_effective(const DissolveOnlyInfos &infos,
                            const size_t /*phaseIndex*/,
                            const double &porosLowerThresh)
{
    return DissolveOnlyUtils::isBelowThresh(porosLowerThresh, infos.porosity);
}

CylinderPCS::CylinderPCS(const std::string &masterPhase,
                         const double& porosInit,
                         const double &rInit)
    : m_masterPhase(masterPhase)
    , m_cylLength(porosInit / (M_PI * rInit * rInit ))
    , k(-1.0)
{
    IPPCheck::assertCheck(m_cylLength >= 0.0);
}

void CylinderPCS::init(const PhaseNameToInfos &phaseInfos, const NucleationData& nuclData)
{
    const auto it = phaseInfos.find(m_masterPhase);
    m_index = std::distance(std::begin(phaseInfos), it);
    IPPCheck::assertCheck(it != phaseInfos.end());
    const PhaseInfo& phaseInfo = it->second;

    const NucleationDissolveOnlyData& nuclInfo = nuclData.at(m_masterPhase);

    const double R = 8.314;
    const double RT = R *  nuclInfo.T;

    // convert from cm3 to to m3
    const double molarVol = phaseInfo.molarVolume * 1e-6;

    k = nuclInfo.HON_IT_water * molarVol / RT;
    IPPCheck::assertCheck(k >= 0.0);
}

bool CylinderPCS::isPCS_effective(const DissolveOnlyInfos &infos,
                                  const size_t phaseIndex,
                                  const double &porosLowerThresh) const
{
    assert(infos.porosity >= 0.0);

    // dense sphere packing
    const double critPoros = 0.74;

    if(infos.porosity < porosLowerThresh )
    {
        return true;
    }
    else if (infos.porosity >= critPoros)
    {
        return false;
    }

    if(m_index == phaseIndex)
    {
        const double SI = *(infos.SI_begin + phaseIndex);
        if(SI > 0.0)
        {
            const double r = std::sqrt(infos.porosity / (M_PI * m_cylLength));
            const double ln_f = k / r;
            const double f = std::exp(ln_f);

            const double slope = 1.0 / critPoros;
            const double linear = slope * infos.porosity;
            const double fAdapt = f * (1.0 - linear) + linear;

            // values lower 1 should be excluded by critPoros check above
            assert(fAdapt >= 1.0);

            const double SI_thresh = std::log10(fAdapt);

            if(SI <= SI_thresh)
            {
                std::cout << "PCS effect!!!!!!!  " << infos.porosity << "\t" << SI << std::endl;
                return true;
            }

        }
    }

    return false;
}


SpherePCS::SpherePCS(const std::string &masterPhase, const double &porosInit, const double &rInit)
{
    IPPCheck::assertCheck(false, "not implemented");
}

void SpherePCS::init(const PhaseNameToInfos &phaseInfos, const NucleationData &nuclData)
{
    IPPCheck::assertCheck(false, "not implemented");
}

bool SpherePCS::isPCS_effective(const DissolveOnlyInfos &infos, const size_t phaseIndex, const double &porosLowerThresh) const
{
    IPPCheck::assertCheck(false, "not implemented");
    return false;
}


template class NucleationDissolveOnlyBase<NoPCS>;
template class NucleationDissolveOnlyBase<CylinderPCS>;
template class NucleationDissolveOnlyBase<SpherePCS>;


}
