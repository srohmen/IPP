#include "tennisjenningscshporositycalc.h"

#include <cassert>
#include <algorithm>

#include "ippexception.h"
#include "cshporosityinfos.h"

namespace IPP
{

const double TennisJenningsCSHPorosityCalc::intrGelPorLD = 0.373;
const double TennisJenningsCSHPorosityCalc::intrGelPorHD = 0.237;
const double TennisJenningsCSHPorosityCalc::referenceVol = 1000;


TennisJenningsCSHPorosityCalc::TennisJenningsCSHPorosityCalc(const double &waterCementRatio,
                                                             const double &degreeHydration)
    : m_compNames(nullptr)
    , m_allPhases(nullptr)
    , m_phaseInfos(nullptr)
    , indexCa(-1)
    , indexSi(-1)
{
    // see paper Tennis + Jennings (2000)
    const double f1 = 3.017;
    const double slope = 1.347;
    const double intercept = 0.538;
    m_massFracCSH_LD = f1 * waterCementRatio * degreeHydration
                       - slope * degreeHydration
                       + intercept;

    IPPCheck::assertCheck(m_massFracCSH_LD >= 0.0 && m_massFracCSH_LD <= 1.0,
                          "degree of hydration is too high for this w/c: led to negative CSH-LD fraction: "
                          + std::to_string(m_massFracCSH_LD));

    // std::cout << "mass fraction CSH-LD: " << m_massFracCSH_LD << std::endl;
}


void TennisJenningsCSHPorosityCalc::setCompNames(const std::vector<std::string>* names)
{
    assert(names);
    assert(m_compNames == nullptr);
    m_compNames = names;

    const std::vector<std::string>::const_iterator itCa =
            std::find(m_compNames->begin(), m_compNames->end(), "Ca");
    IPPCheck::assertCheck(itCa != m_compNames->end(),
                          "could not find Calcium in components but try to use CSH diffusivity model");
    indexCa = itCa - m_compNames->begin();


    const std::vector<std::string>::const_iterator itSi =
            std::find(m_compNames->begin(), m_compNames->end(), "Si");
    IPPCheck::assertCheck(itSi != m_compNames->end(),
                          "could not find Silicon in components but try to use CSH diffusivity model");
    indexSi = itSi - m_compNames->begin();
}

void TennisJenningsCSHPorosityCalc::setPhaseNames(const std::vector<std::string> *names)
{
    assert(names);
    assert(m_allPhases == nullptr);
    m_allPhases = names;
}

void TennisJenningsCSHPorosityCalc::setPhaseInfos(const PhaseNameToInfos *phaseInfos)
{
    assert(phaseInfos);
    assert(m_phaseInfos == nullptr);
    m_phaseInfos = phaseInfos;

    for(const std::pair<const std::string, PhaseInfo>& infoPair : *m_phaseInfos)
    {
        const std::string& phaseName = infoPair.first;
        const PhaseInfo& info = infoPair.second;
        IPPCheck::assertCheck(info.gfw > 0.0, "molar mass is zero for: " + phaseName);
    }
}

const double& TennisJenningsCSHPorosityCalc::getMassFractionCSH_LD() const
{
    return m_massFracCSH_LD;
}

void TennisJenningsCSHPorosityCalc::add_CSH_Phase(const std::string &name)
{
    m_CSHphaseNames.insert(name);
}


AbstractPorosityCalc::SimplePorosityInfosPtr TennisJenningsCSHPorosityCalc::calc(const PorosityCalcInput& input) const
{
    std::unique_ptr<CSHPorosityInfos> output(new CSHPorosityInfos);

    const BeginEndIt& range = input.phaseAmounts;

    double dryVolCSH_LD = 0.0;
    double dryVolCSH_HD = 0.0;

    double& volNonPerm = output->volNonPerm;
    volNonPerm = 0.0;

    double sumCa = 0.0;
    double sumSi = 0.0;

    size_t counter = 0;
    for(ValueIt it = range.begin; it != range.end; ++it)
    {
        const std::string& phaseName = m_allPhases->at(counter);
        const PhaseInfo& info = m_phaseInfos->at(phaseName);

        const double& phaseAmount = *it;
        const double phaseVol = phaseAmount * info.molarVolume;

        if(m_CSHphaseNames.find(phaseName) != m_CSHphaseNames.end())
        {
            const double massCSH = phaseAmount * info.gfw;

            const double massCSH_LD = m_massFracCSH_LD * massCSH;
            const double massCSH_HD = (1.0 - m_massFracCSH_LD) * massCSH;

            // IMPORTANT: the molar mass and volume must
            // be the one for dry CSH i.e. unexpanded by gel porosity
            const double nCSH_LD = massCSH_LD / info.gfw;
            const double nCSH_HD = massCSH_HD / info.gfw;

            // TODO: check if molar volume is valid for LD vs. HD
            // maybe using seperation factor similar to mass frac?
            const double volCSH_LD = nCSH_LD * info.molarVolume;
            const double volCSH_HD = nCSH_HD * info.molarVolume;

            dryVolCSH_LD += volCSH_LD;
            dryVolCSH_HD += volCSH_HD;


            const double nCSH = nCSH_LD + nCSH_HD;

            assert(info.stoich.size() > indexCa);
            const double& stoichCa = info.stoich[indexCa];
            sumCa += stoichCa * nCSH;

            assert(info.stoich.size() > indexSi);
            const double& stoichSi = info.stoich[indexSi];
            sumSi += stoichSi * nCSH;
        }
        else
        {
            volNonPerm += phaseVol;
            assert(volNonPerm + input.inertFraction * referenceVol <= referenceVol);
        }

        ++counter;
    }


    volNonPerm += input.inertFraction * referenceVol;

    if(volNonPerm > referenceVol)
    {
        const std::string errMsg = "non diffusive volume is bigger than reference volume: "
                                   + std::to_string(volNonPerm) + " > " + std::to_string(referenceVol);
        throw std::runtime_error(errMsg);
    }


    const double volTotalDrySolid = volNonPerm + dryVolCSH_LD + dryVolCSH_HD;
    const double volFracTotalDrySolid = volTotalDrySolid / referenceVol;
    output->porosityTotal = 1.0 - volFracTotalDrySolid;

    if(output->porosityTotal < 0.0 || output->porosityTotal > 1.0)
    {
        const std::string errMsg =  "total porosity not in valid range: "
                                    + std::to_string(output->porosityTotal);
        throw std::runtime_error(errMsg);
    }


    if(dryVolCSH_LD > 0.0 || dryVolCSH_HD > 0.0)
    {
        const double satVolCSH_LD = dryVolCSH_LD / (1.0 - intrGelPorLD);
        const double satVolCSH_HD = dryVolCSH_HD / (1.0 - intrGelPorHD);

        IPPCheck::assertCheck(satVolCSH_LD >= 0.0);
        IPPCheck::assertCheck(satVolCSH_HD >= 0.0);

        double& volSatCSH = output->volSatCSH;
        volSatCSH = satVolCSH_LD + satVolCSH_HD;

        if(volSatCSH > referenceVol)
        {
            throw IPPException("saturated CSH exceeds reference vol: " + std::to_string(volSatCSH));
        }

        const double volFracLD = satVolCSH_LD / volSatCSH;
        output->volFracCSH_HD = 1.0 - volFracLD;

        const double volForCSH = referenceVol - volNonPerm;
        const double volCapillary = volForCSH - volSatCSH;
        const double porosityCapillary = volCapillary / referenceVol;

        if(porosityCapillary < 0.0)
        {
            std::cout << ("WARNING: Got negative (capillary) porosity( " + std::to_string(porosityCapillary)
                          + " ) when subtracting gel pores") << std::endl;

            std::cout << "phase\tn[mol]\tvol\n";

            size_t iPhase = 0;
            for(auto it = range.begin; it != range.end; ++it)
            {
                const std::string& name = m_allPhases->at(iPhase);
                const PhaseInfo& info = m_phaseInfos->at(name);
                std::cerr << name << "\t" << *it << "\t" << *it * info.molarVolume << std::endl;

                ++iPhase;
            }
            std::cout << "Inert\t0\t" << input.inertFraction << std::endl;

        }

        output->porosityCapillary = porosityCapillary;

        output->ratioCaSi = sumCa / sumSi;
        assert(output->ratioCaSi > 0.0);
    }
    else
    {
        const double volFracNonPerm = volNonPerm / referenceVol;
        output->porosityCapillary = 1.0 - volFracNonPerm;

        output->volSatCSH = 0.0;
        output->ratioCaSi = 0.0;
    }



    return output;

}

}
