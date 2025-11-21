#include "phreeqcinputfilegenerator.h"

#include "phreeqcconstants.h"
#include "ippconstants.h"
#include "phreeqcdefaultconvergencevalues.h"

#include "phasenametodissprecipbehaviour.h"
#include "abstractdissolveonlycalc.h"

namespace IPP
{

namespace
{

void generateTargetSIPhaseBody(const AbstractDissPrecipOnlyInfo& phaseBehav,
                               const std::map<std::string, double>& phaseNameToAmount,
                               std::stringstream& ss)
{
    for(const auto& phase : phaseNameToAmount)
    {
        const std::string& phaseName = phase.first;
        const DissPrecipBehaviour& dissPrec = phaseBehav.getBehaviour(phaseName);

        std::string extension = "";

        if(DPB_DissOnly == dissPrec)
        {
            extension = "dissolve_only";
        }
        else if(DPB_PrecOnly == dissPrec)
        {
            extension = "precipitate_only";
        }
        else if(DPB_Normal == dissPrec)
        {
            extension = "";
        }
        else
        {
            throw std::runtime_error("Unknown dissolution/precipitation behaviour: "
                                     + std::to_string(dissPrec));
        }


        const double targetSI = 0.0;

        ss << "\t" << phaseName
           << "\t" << std::to_string(targetSI)
           << "\t" << phase.second // amount
           << "\t" << extension
           << std::endl;

    }
}



void generateHeader(std::stringstream& ss)
{
    ss << "TITLE no name\n";
}


static void generateEqPhaseBodyExcludeNonEq(std::map<std::string, double> phaseNameToAmount, // copy intended
                                            const Composition::PhaseDefVec& phases,
                                            const AbstractDissPrecipOnlyInfo& phaseBehav,
                                            const std::set<std::string>& nonEqPhases,
                                            std::stringstream& ss)
{
    for(const Composition::PhaseDefinition& phase : phases)
    {
        if(nonEqPhases.find(phase.name) == nonEqPhases.end())
        {
            assert(phaseNameToAmount.find(phase.name) != phaseNameToAmount.end());
            phaseNameToAmount[phase.name] = phase.amount;
        }
    }

    generateTargetSIPhaseBody(phaseBehav, phaseNameToAmount, ss);
}

void generateKineticPhaseBody(const Composition::PhaseDefVec& phases,
                              const std::set<std::string>& kineticPhases,
                              std::stringstream& ss)
{
    constexpr double tolerance = 1.0E-3;
    const std::string tolStr = std::to_string(tolerance);
    const std::string tolLine = "-tol " + tolStr +  "\n";

    for(const Composition::PhaseDefinition& phase : phases)
    {
        if(kineticPhases.find(phase.name) != kineticPhases.end())
        {
            const double amount = phase.amount;

            ss << phase.name << std::endl;
            ss << "-m\t" << amount << std::endl;
            ss << "-m0\t" << amount << std::endl;
            ss << tolLine;
        }
    }
}

bool isKineticsNeeded(const Composition::PhaseDefVec& phases,
                      const std::set<std::string>& kineticPhases)
{
    for(const Composition::PhaseDefinition& phase : phases)
    {
        if(kineticPhases.find(phase.name) != kineticPhases.end())
        {
            return true;
        }
    }

    return false;
}

void generateWaterLimitedRates(const std::set<std::string>& waterLimitedPhases, std::stringstream& ss)
{
    if(waterLimitedPhases.empty() == false)
    {
        ss << "RATES" << std::endl;
        for(const std::string& phaseName : waterLimitedPhases)
        {
            ss << phaseName << std::endl;
            ss << "-start\n";
            ss << "10 mole = 0\n";
            ss << "20 if(M <= 0) THEN GOTO 300\n";

//            ss << "100 cCa = TOT(\"Ca\") / SOLN_VOL \n";
            ss << "105 mass_water = TOT(\"Water\")\n";

            // ss << "110 if(cCa > 0.03) THEN GOTO 300\n";
            ss << "111 if(mass_water < 0.001) THEN GOTO 300\n";

            ss << "113 mole_water_poros = (mass_water - 0.9) * 10000 / GFW(\"H2O\")\n";
            ss << "120 rate = mole_water_poros\n";
            // ss << "125 rate = 100000\n";
            ss << "130 mole = rate * Time\n";
            ss << "140 if(mole > M) THEN mole = M\n";
            ss << "300 Save mole\n";
            ss << "-end\n";
        }
    }
}

void generateKineticConvergence(std::stringstream& ss)
{
    ss << "-steps 0\n";
    ss << "-cvode true\n";
    ss << "-cvode_order 5\n";

    // ss << "-step_divide 10\n";
    // -cvode_order 2
    // -cvode_steps 500
    // -bad_step_max 200
}

} // end of anonymous namespace

namespace PhreeqcInputFileGenerator
{

void generateEqPhaseBody(const Composition::PhaseDefVec& phases,
                         const AbstractDissPrecipOnlyInfo& phaseBehav,
                         std::stringstream& ss)
{
    std::map<std::string, double> phaseNameToAmount;

    for(const Composition::PhaseDefinition& phase : phases)
    {
        assert(phaseNameToAmount.find(phase.name) == phaseNameToAmount.end());
        phaseNameToAmount[phase.name] = phase.amount;
    }

    generateTargetSIPhaseBody(phaseBehav, phaseNameToAmount, ss);
}


void generateElementConcLine(const std::string& elemName,
                             const double& conc,
                             const std::string& specName,
                             std::stringstream& ss)
{
    ss << "\t" << elemName << "\t" << conc;
    if(specName.empty() == false)
    {
        ss << " as " << specName;
    }
    ss << std::endl;
}


void generateElementConcBody(const Composition& comp,
                             std::map<ElementSpecies, double> elemNameToAmount,
                             std::stringstream& ss)
{
    // copy of map is intended

    for(const Composition::ElementConcentration& compName : comp.elemConcVec)
    {
        elemNameToAmount[compName.elemSpecies] = compName.conc;
    }

    for(const auto& compName : elemNameToAmount)
    {
        generateElementConcLine(compName.first.elementName, compName.second, compName.first.speciesName, ss);
    }
}


void generateSolutionBody(const Composition& comp,
                          const std::map<ElementSpecies, double>& elemNameToAmount,
                          std::stringstream& ss)
{

    assert(comp.units.empty() == false);
    ss << "\t-units\t" << comp.units << std::endl;
    ss << "\t-temp\t" << comp.temp << std::endl;
    ss << "\t-pressure\t" << comp.pressure << std::endl;
    ss << "\t-pH\t" << comp.pH.value << " " << comp.pH.extension << std::endl;
    ss << "\t-pe\t" << comp.pe.value << " " << comp.pe.extension << std::endl;

    generateElementConcBody(comp, elemNameToAmount, ss);
}

void generateSolidSolutionBody(const std::vector<IPPConfig::SolidSolutionDefinition>& solidSolutionVec,
                               const Composition::PhaseDefVec& phases,
                               std::stringstream& ss)
{
    for(size_t iSS = 0; iSS < solidSolutionVec.size(); ++iSS)
    {
        const IPPConfig::SolidSolutionDefinition& def = solidSolutionVec[iSS];
        const std::string& name = def.name;
        ss << name << std::endl;

        for(const std::string& compName : def.phaseNames)
        {
            // find amount if defined
            double amount = 0;
            for(const Composition::PhaseDefinition& phase : phases)
            {
                const std::string& paseName = phase.name;
                if(compName == paseName)
                {
                    amount = phase.amount;
                    break;
                }
            }

            ss << "\t-comp\t" << compName << "\t" << amount << std::endl;
        }


    }

}

void generateDefinedComp(const IPPConfig& config,
                         const Composition& comp,
                         const std::map<std::string, double>& eqPhaseNameToAmount,
                         const AbstractDissPrecipOnlyInfo& dissPrecBehav,
                         const std::set<std::string>& nonEqPhases,
                         const std::map<ElementSpecies, double>& allElemConc,
                         const std::string& phreeqcID,
                         unsigned int& featureMask,
                         std::stringstream& ss)
{
    ss << "SOLUTION " << phreeqcID << "\t" << comp.name << std::endl;
    // append solution body
    generateSolutionBody(comp, allElemConc, ss);
    featureMask |= PF_SOLUTION;



    ss << "EQUILIBRIUM_PHASES " << phreeqcID << "\t" << comp.name << std::endl;
    generateEqPhaseBodyExcludeNonEq(eqPhaseNameToAmount, comp.phases,
                                    dissPrecBehav, nonEqPhases, ss);
    featureMask |= PF_EQUILIBRIUM_PHASES;



    ss << "KINETICS " << phreeqcID << "\t" << comp.name << std::endl;
    if(isKineticsNeeded(comp.phases, config.kineticPhases))
    {
        generateKineticPhaseBody(comp.phases,
                                 config.kineticPhases,
                                 ss);

        generateKineticConvergence(ss);

        featureMask |= PF_KINETICS;
    }

    ss << "SOLID_SOLUTIONS " << phreeqcID << "\t" << comp.name << std::endl;
    if(config.solidSolutionPhases.empty() == false)
    {
        generateSolidSolutionBody(config.solidSolutionVec, comp.phases, ss);
        featureMask |= PF_SOLID_SOLUTION;
    }
}

void generateUndefinedComp(const IPPConfig& config,
                           const std::map<std::string, double>& phaseNameToAmount,
                           const AbstractDissPrecipOnlyInfo& dissPrecBehav,
                           const std::map<ElementSpecies, double>& allElemConc,
                           const std::string& phreeqcID,
                           unsigned int& featureMask,
                           std::stringstream& ss)
{
    ss << "SOLUTION " << phreeqcID << "\tEMPTY\n";
    ss << "\t-units mol/L\n";
    ss << "\t-pH 7 charge\n";
    // append solution body
    for(const auto& compName : allElemConc)
    {
        generateElementConcLine(compName.first.elementName,
                                compName.second,
                                compName.first.speciesName, ss);
    }
    featureMask |= PF_SOLUTION;

    ss << "EQUILIBRIUM_PHASES " << phreeqcID << std::endl;
    generateTargetSIPhaseBody(dissPrecBehav, phaseNameToAmount, ss);
    featureMask |= PF_EQUILIBRIUM_PHASES;

    // dummys for filling cells correctly
    ss << "KINETICS " << phreeqcID << std::endl;
    // do not enable KINETICS in features

    ss << "SOLID_SOLUTIONS " << phreeqcID << std::endl;
    if(config.solidSolutionPhases.empty() == false)
    {
        const Composition::PhaseDefVec emptyPhases;
        generateSolidSolutionBody(config.solidSolutionVec, emptyPhases, ss);
        featureMask |= PF_SOLID_SOLUTION;
    }
}


void generateInertSoluteDatabaseExt(const std::vector<IPPConfig::Results::DiffEffInfo>& diffEffInfos,
                                    std::stringstream& ss)
{
    if(diffEffInfos.empty() == false)
    {
        ss << "SOLUTION_MASTER_SPECIES\n";
        for(size_t i = 0; i < diffEffInfos.size(); ++i)
        {
            const IPPConfig::Results::DiffEffInfo& info = diffEffInfos[i];
            ss << info.tracerName << "\t" << info.tracerName << "\t0\t1\t0" << std::endl;
        }

        ss << "SOLUTION_SPECIES\n";
        for(size_t i = 0; i < diffEffInfos.size(); ++i)
        {
            const IPPConfig::Results::DiffEffInfo& info = diffEffInfos[i];
            ss << info.tracerName << " = " << info.tracerName << std::endl;
        }
    }
}


void generateConvergenceHacksDefault(std::stringstream& ss)
{
//    ss << "SOLUTION_SPECIES\n";
//    ss << "\tH2O + 0.01e- = H2O-0.01\n";
//    ss << "\tlog_k   -9.0\n";

     generateConvergenceHacks(PhreeqcDefaultConvergenceValues::iterations,
                              PhreeqcDefaultConvergenceValues::tolerance,
                              PhreeqcDefaultConvergenceValues::convergence_tolerance,
                              PhreeqcDefaultConvergenceValues::diagonalScaling,
                              ss);
}

void generateConvergenceHacksInitRun(std::stringstream& ss)
{
//    ss << "SOLUTION_SPECIES\n";
//    ss << "\tH2O + 0.01e- = H2O-0.01\n";
//    ss << "\tlog_k   -9.0\n";

     generateConvergenceHacks(PhreeqcDefaultConvergenceValues::iterations * 2,
                              PhreeqcDefaultConvergenceValues::tolerance * 0.1,
                              PhreeqcDefaultConvergenceValues::convergence_tolerance,
                              PhreeqcDefaultConvergenceValues::diagonalScaling,
                              ss);

//    generateConvergenceHacksDefault(ss);
}


void generateConvergenceHacks(const size_t iterations,
                              const double& tolerance,
                              const double& convergence_tolerance,
                              const bool diagonalScaling,
                              std::stringstream& ss)
{
    ss << "KNOBS\n";
    ss << "\t-iterations " << iterations << std::endl;;
    ss << "\t-tolerance " << tolerance << std::endl;
    ss << "\t-convergence_tolerance " << convergence_tolerance << std::endl;


    if(diagonalScaling == true)
    {
        ss << "\t-diagonal_scale true\n";
    }
    else
    {
        ss << "\t-diagonal_scale false\n";
    }

}


static void generateSelectedOutput(const IPPConfig& config,
                                   const std::vector<ElementSpecies>& allElems,
                                   const std::string& allElemsStr,
                                   const std::string& allEqPhasesStr,
                                   const std::string& allSolidSolutionPhasesStr,
                                   const std::string& allKinPhasesStr,
                                   const std::vector<std::string>& nucleationMonomers,
                                   const bool enableTotals,
                                   const bool enableSaturationIndices,
                                   std::stringstream& ss)
{
    ss << "SELECTED_OUTPUT " << PhreeqcConstants::SO_PhasesAmount << std::endl;
    ss << "\t-reset false\n";
    ss << "\t-file out_" + std::to_string(PhreeqcConstants::SO_PhasesAmount) + "\n";
    ss << "\t-high_precision\n";
    ss << "\t-equilibrium_phases\t" << allEqPhasesStr << std::endl;
    ss << "\t-kinetic_reactants\t" << allKinPhasesStr << std::endl;
    ss << "\t-solid_solutions\t" << allSolidSolutionPhasesStr << std::endl;


    if(enableTotals)
    {
        ss << "SELECTED_OUTPUT " << PhreeqcConstants::SO_ComponentTotal << std::endl;
        ss << "\t-reset false\n";
        ss << "\t-file out_" + std::to_string(PhreeqcConstants::SO_ComponentTotal) + "\n";
        ss << "\t-high_precision\n";
        ss << "\t-user_punch true\n";

        ss << "USER_PUNCH " << PhreeqcConstants::SO_ComponentTotal << std::endl;
        ss << "-heading H20 H O cb ";
        ss << allElemsStr << std::endl;
        ss << "-start\n";


        ss << "10 m_H2O = TOTMOLE(\"water\")\n";
        ss << "11 PUNCH m_H2O\n";
        ss << "20 m_H = SYS(\"H\") - 2*m_H2O\n";
        ss << "21 PUNCH m_H\n";
        ss << "30 m_O = SYS(\"O\") - m_H2O\n";
        ss << "31 PUNCH m_O\n";
        ss << "40 PUNCH CHARGE_BALANCE\n";

        for(size_t i = 0; i < allElems.size(); ++i)
        {
            const size_t line = 10*i + 50;
            const std::string& elemName = allElems[i].elementName;
            const std::string elemVar = "m_" + elemName;
            ss << line << " " << elemVar << " = SYS(\"" << elemName << "\")\n";
            ss << line + 1 << " PUNCH " << elemVar << std::endl;
        }
        ss << "-end\n";
    }



    if(enableSaturationIndices)
    {
        const std::string allPhasesStr = allEqPhasesStr + " " + allSolidSolutionPhasesStr;
        ss << "SELECTED_OUTPUT " << PhreeqcConstants::SO_PhasesSaturationIndices << std::endl;
        ss << "\t-reset false\n";
        ss << "\t-file out_" + std::to_string(PhreeqcConstants::SO_PhasesSaturationIndices) + "\n";
        ss << "\t-high_precision\n";
        ss << "\t-saturation_indices\t" << allPhasesStr << std::endl;
    }

    const std::vector<IPPConfig::Results::NameCommand>& auxData = config.results.auxData;
    if(auxData.empty() == false)
    {
        ss << "SELECTED_OUTPUT " << PhreeqcConstants::SO_Aux << std::endl;
        ss << "\t-reset false\n";
        ss << "\t-file out_" + std::to_string(PhreeqcConstants::SO_Aux) + "\n";
        ss << "\t-high_precision\n";

        for(size_t i = 0; i < auxData.size(); ++i)
        {
            const IPPConfig::Results::NameCommand& nameCommand = auxData[i];
            ss << "\t-" + nameCommand.command + "\n";
        }
    }


    if(nucleationMonomers.empty() == false)
    {
        ss << "SELECTED_OUTPUT " << PhreeqcConstants::SO_NucleationMonomers << std::endl;
        ss << "\t-reset false\n";
        ss << "\t-file out_" + std::to_string(PhreeqcConstants::SO_NucleationMonomers) + "\n";
        ss << "\t-high_precision\n";

        ss << "\t-molalities";
        for(size_t i = 0; i < nucleationMonomers.size(); ++i)
        {
            const std::string& name = nucleationMonomers[i];
            ss << " " + name;

        }
        ss << std::endl;
    }
}

std::tuple<bool, AbstractDissPrecipOnlyInfo*>
findInitialDissPrecBehav(const IPPConfig& config,
                         const std::vector<std::string>& allPhases)
{
    AbstractDissPrecipOnlyInfo& persistentBehav = *config.dissPrecBehav;

    if(config.dissolveOnlyFunc->isNucleation())
    {
        PhaseNameToDissPrecipBehaviour* dissOnlyAll = new PhaseNameToDissPrecipBehaviour;

        for(const std::string& name : allPhases)
        {
            dissOnlyAll->addData(name, DPB_DissOnly);
        }

        return std::make_tuple(true, dissOnlyAll);
    }
    else
    {
        return std::make_tuple(false, &persistentBehav);
    }
}

void generateInitString(const IPPConfig& config,
                        const std::vector<std::string>& allPhases,
                        const bool enableTotals,
                        const bool enableSaturationIndex,
                        CellIdToPhreeqcFeatureMask& enabledFeatures,
                        size_t &genericCellID,
                        std::stringstream& compositionSS,
                        std::stringstream& outputSS)
{
    const std::vector<ElementSpecies>& allElems = config.getAllElements();
    std::string allElemsStr;
    std::map<ElementSpecies, double> allElemConc;
    for(size_t i = 0; i < allElems.size(); ++i)
    {
        const ElementSpecies& elemName = allElems[i];
        allElemsStr += elemName.elementName + " ";
        allElemConc[elemName] = 0.0;
    }


    std::string allEqPhasesStr, allKinPhasesStr, allSolidSolutionPhasesStr;
    std::map<std::string, double> eqPhaseNameToAmount;

    std::set<std::string> nonEqPhases;

    for(const std::string& phaseName : allPhases)
    {
        const bool isKineticPhase = config.kineticPhases.find(phaseName)
                != config.kineticPhases.end();
        if(isKineticPhase)
        {
            // kinetic phase
            allKinPhasesStr += phaseName + " ";
            nonEqPhases.insert(phaseName);
        }
        else
        {
            const bool isSolidSolutionPhase = config.solidSolutionPhases.find(phaseName)
                    != config.solidSolutionPhases.end();
            if(isSolidSolutionPhase)
            {
                // solid solution phase
                allSolidSolutionPhasesStr += phaseName + " ";
                nonEqPhases.insert(phaseName);
            }
            else
            {
                // eq phase
                eqPhaseNameToAmount[phaseName] = 0.0;
                allEqPhasesStr += phaseName + " ";
            }
        }

    }


    generateHeader(compositionSS);

    generateWaterLimitedRates(config.waterLimitedPhases, compositionSS);


    AbstractDissPrecipOnlyInfo* initialBehav = nullptr;
    bool mustDeleteBehav = false;
    std::tie(mustDeleteBehav, initialBehav) = findInitialDissPrecBehav(config, allPhases);


    {
        const std::vector<Domain>& domains = config.domains;

        for(size_t iDomain = 0; iDomain < domains.size(); ++iDomain)
        {
            const Domain& domain = domains[iDomain];

            assert(domain.composition);
            const Composition& comp = *(domain.composition);
            const std::string phreeqcID = std::to_string(iDomain);

            unsigned int& featureMask = enabledFeatures[iDomain];

            generateDefinedComp(config, comp, eqPhaseNameToAmount, *initialBehav,
                                nonEqPhases, allElemConc, phreeqcID,
                                featureMask, compositionSS);
        }
    }


    const size_t nxyz = config.nx * config.ny * config.nz;


    {
        // Dummy cell for tricking phreeqc to include the elements at all
        const std::string range = std::to_string(nxyz);
        std::map<ElementSpecies, double> dummyElemConc;
        for(const ElementSpecies& elemName : allElems)
        {
            dummyElemConc[elemName] = 1.0E-3;
        }

        compositionSS << "SOLUTION " << range << "\tDUMMY\n";
        compositionSS << "\t-units mol/L\n";
        compositionSS << "\t-pH 7 charge\n";
        // append solution body
        for(const auto& compName : dummyElemConc)
        {
            generateElementConcLine(compName.first.elementName,
                                    compName.second,
                                    compName.first.speciesName, compositionSS);
        }
    }

    {
        // append one empty (water or generic) cell at the very end
        // try to find the generic definition
        const std::vector<Composition>::const_iterator it = std::find_if(config.compositions.begin(),
                                                                         config.compositions.end(),
                                   [](const Composition& comp){ return comp.name == "generic";}
                );

        genericCellID = nxyz+1;
        const std::string genericCellIDStr = std::to_string(genericCellID);

        unsigned int& featureMask = enabledFeatures[genericCellID];
        if(it != config.compositions.end())
        {
            const Composition& comp = *it;
            generateDefinedComp(config, comp, eqPhaseNameToAmount, *initialBehav, nonEqPhases,
                                allElemConc, genericCellIDStr, featureMask, compositionSS);
        }
        else
        {
            generateUndefinedComp(config, eqPhaseNameToAmount, *initialBehav, allElemConc,
                                  genericCellIDStr, featureMask, compositionSS);
        }
    }

    if(mustDeleteBehav)
    {
        delete initialBehav;
    }


    // we need no more the inert solutes DB extension because this is handled in dummy reaction module
    //    const std::vector<IPPConfig::Results::DiffEffInfo>& diffEffInfos = config.results.diffEffInfos;
    //    generateInertSoluteDatabaseExt(diffEffInfos, ss);


    // define output data fields
    std::vector<std::string> nucleationMonomers;
    config.dissolveOnlyFunc->findNucleationMonomers(nucleationMonomers);

    generateSelectedOutput(config, allElems, allElemsStr,
                           allEqPhasesStr, allSolidSolutionPhasesStr,
                           allKinPhasesStr, nucleationMonomers,
                           enableTotals, enableSaturationIndex, outputSS);

    generateConvergenceHacksDefault(outputSS);

}

}

}
