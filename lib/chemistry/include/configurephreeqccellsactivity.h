#ifndef CONFIGUREPHREEQCCELLSACTIVITY_H
#define CONFIGUREPHREEQCCELLSACTIVITY_H

#include <iostream>
#include <vector>
#include <cstddef>
#include <memory>
#include <cassert>



namespace IPP
{

class LocalAccessPhreeqcRM;
class IsInsideVolume;
class SignificantChangePrediction;
class PhaseIdPerElement;
class HandleInterfaceSwitch;


namespace TestChains
{
class TestChain;
}

class ConfigurePhreeqcCellsActivity
{
public:
    ConfigurePhreeqcCellsActivity(const double &tolerance,
                                  const size_t forceCalcIt,
                                  const SignificantChangePrediction* changePredict);

    ~ConfigurePhreeqcCellsActivity();

    void init(const size_t nCells, const size_t nComps,
              const size_t nPrimaryComps, const size_t nPhases,
              const std::vector<double>& currConc,
              const PhaseIdPerElement& phaseIdPerElement,
              std::vector<char>& isNonPermNonInterfaceCell,
              std::vector<char>& isInterfaceCell);

    HandleInterfaceSwitch& getInterfaceSwitch();
    void setStaticDisabledCells(const std::vector<size_t> &cellIDs);
    void setGloballyEnabledState(const bool isEnabled);
    void setGloballyDisabledState(const bool isDisabled);

    void run(const size_t it,
             const std::vector<char>& oldEnabled,
             std::vector<size_t> &newEnabled,
             std::vector<size_t> &newDisabled);


    std::vector<double>& getLastPhreeqcEqConc();
    std::vector<double>& getLastPhreeqcPhaseAmount();
    std::vector<double>& getLastPhreeqcTotalAmount();
    std::vector<double>& getConcDiffCummulative();

private:
    bool mustRecalc(const size_t it, const size_t iCellLocal);
    bool hasSignificantChange(const size_t iCell);


    // prechecks
    std::vector<char> m_staticDisabledCells;

    std::unique_ptr<HandleInterfaceSwitch> m_interfaceSwitch;

    bool m_isGloballyEnabled;
    bool m_isGloballyDisabled;

    std::vector<double> m_lastPhreeqcEqConc;
    const double m_tolerance;

    const size_t m_forceCalcIt;


    const SignificantChangePrediction* m_changePredict;
    std::vector<double> m_lastPhreeqcPhaseAmount;
    std::vector<double> m_lastPhreeqcTotalAmount;
    std::vector<double> m_concDiffCummulative;


    bool m_hasSignificancyChecks;
    std::vector<TestChains::TestChain*> m_chain;


};

}

#endif // CONFIGUREPHREEQCCELLSACTIVITY_H
