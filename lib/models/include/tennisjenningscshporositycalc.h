#ifndef TENNISJENNINGSCSHPOROSITYCALC_H
#define TENNISJENNINGSCSHPOROSITYCALC_H

#include "abstractporositycalc.h"

#include <string>
#include <set>


namespace IPP
{

struct AbstractPorosityInfos;

class TennisJenningsCSHPorosityCalc : public AbstractPorosityCalc
{
public:
    TennisJenningsCSHPorosityCalc(const double& waterCementRatio,
                                  const double& degreeHydration);


    void setPhaseInfos(const PhaseNameToInfos *phaseInfos) override;
    void setPhaseNames(const std::vector<std::string> *names) override;
    void setCompNames(const std::vector<std::string> *names) override;

    virtual SimplePorosityInfosPtr calc(const PorosityCalcInput& input) const override;


    const double& getMassFractionCSH_LD() const;
    void add_CSH_Phase(const std::string& phaseName);

private:
    static const double intrGelPorLD;
    static const double intrGelPorHD;
    static const double referenceVol; // cm3

    double m_massFracCSH_LD;

    const std::vector<std::string>* m_compNames;
    const std::vector<std::string>* m_allPhases;
    std::set<std::string> m_CSHphaseNames;
    const PhaseNameToInfos* m_phaseInfos;

    size_t indexCa;
    size_t indexSi;
};

}

#endif // TENNISJENNINGSCSHPOROSITYCALC_H
