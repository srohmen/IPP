#include "phreeqcprerunhandler.h"

#include <fstream>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include <PhreeqcRM.h>
#include <IPhreeqcPhast.h>

#include "phreeqcinputfilegenerator.h"
#include "ippconstants.h"
#include "phreeqcconstants.h"
#include "ippexception.h"
#include "isincontainer.h"
#include "mpitools.h"
#include "phaseidperelement.h"
#include "ippstream.h"

namespace IPP
{

namespace PhreeqcPrerunHandler
{

static const size_t s_selectedOutMolVol = 0;
static const size_t s_selectedOutPhaseFomulas = 1;
// static const size_t s_selectedOutMolWeight = 2;
static const size_t s_soFirstDataRow = 1;

static void generateMolVolString(std::stringstream& inputStrm)
{
    inputStrm << "USER_PUNCH " + std::to_string(s_selectedOutMolVol) << std::endl;
    inputStrm << "\t-headings count names vms \n";
    inputStrm << "\t-start\n";

    inputStrm << "10 sysReturn = SYS( \"phases\", count, names$, type$, phase_si )\n";

    inputStrm << "20 PUNCH count\n";

    inputStrm << "30 allNames$ = \"\"\n";
    inputStrm << "35 allVMs$ = \"\"\n";

    inputStrm << "40 FOR i = 1 TO count STEP 1\n";

    inputStrm << "50 name$ = names$(i)\n";
    inputStrm << "60 allnames$ = allnames$+name$+\"!\"\n";

    inputStrm << "70 vm$ = trim(STR$(PHASE_VM(name$)))\n";
    inputStrm << "80 allVMs$ = allVMs$+vm$+\"!\"\n";

    inputStrm << "90 NEXT i\n";

    inputStrm << "92 PUNCH allNames$\n";
    inputStrm << "95 PUNCH allVMs$\n";


    inputStrm << "\t-end\n";


    inputStrm << "SELECTED_OUTPUT " + std::to_string(s_selectedOutMolVol) << std::endl;
    inputStrm << "\t-reset false\n";
    inputStrm << "\t-file sel_out_molvol\n";
    inputStrm << "\t-high_precision\n";
    inputStrm << "\t-user_punch true\n";


}

static void generateFormulaeString(const std::string& phaseName, const size_t iPhase, std::stringstream& inputStrm)
{
    const size_t iSo = s_selectedOutPhaseFomulas + iPhase;
    inputStrm << "USER_PUNCH " + std::to_string(iSo) << std::endl;
    inputStrm << "\t-headings count names coefs gfw\n";
    inputStrm << "\t-start\n";

//    inputStrm << "10 PUNCH 1\n";
//    inputStrm << "20 PUNCH \"test1\"\n";
//    inputStrm << "30 PUNCH 2\n";

    inputStrm << "10 formReturn$ = PHASE_FORMULA$( \"" + phaseName + "\", count, elt$, coefs )\n";

    inputStrm << "20 PUNCH count\n";

    inputStrm << "30 allElems$ = \"\"\n";
    inputStrm << "35 allCoefs$ = \"\"\n";

    inputStrm << "40 FOR i = 1 TO count STEP 1\n";

    inputStrm << "50 name$ = elt$(i)\n";
    inputStrm << "60 allElems$ = allElems$+name$+\"!\"\n";
    inputStrm << "65 coef = coefs(i)\n";
    inputStrm << "70 coefStr$ = trim(STR$(coef))\n";
    inputStrm << "80 allCoefs$ = allCoefs$+coefStr$+\"!\"\n";

    inputStrm << "90 NEXT i\n";


    inputStrm << "91 weight = GFW(formReturn$)\n";

    inputStrm << "95 PUNCH allElems$\n";
    inputStrm << "96 PUNCH allCoefs$\n";
    inputStrm << "97 PUNCH weight\n";

    inputStrm << "\t-end\n";


    inputStrm << "SELECTED_OUTPUT " + std::to_string(iSo) << std::endl;
    inputStrm << "\t-reset false\n";
    inputStrm << "\t-file sel_out_" + std::to_string(iPhase) + "\n";
    inputStrm << "\t-high_precision\n";
    inputStrm << "\t-user_punch true\n";

}

static void generateFormulaeString(const std::vector<std::string>& phaseNames, std::stringstream& inputStrm)
{
    for(size_t iPhase = 0; iPhase < phaseNames.size(); ++iPhase)
    {
        const std::string& phaseName = phaseNames[iPhase];
        generateFormulaeString(phaseName, iPhase, inputStrm);
    }
}

void findElementPhaseNamesAndMolVols(PhreeqcRM& phreeqc,
                                     std::vector<std::string>& phaseNames,
                                     std::map<std::string, double>& molVols)
{
    const size_t countCol = 0;
    const size_t phasesCol = 1;
    const size_t volCol = 2;

    phreeqc.SetCurrentSelectedOutputUserNumber(s_selectedOutMolVol);

    const std::vector<IPhreeqcPhast*>& internalPhreeqcVec = phreeqc.GetWorkers();
    IPhreeqcPhast* internalPhreeqc = internalPhreeqcVec.front();

    const size_t nRows = internalPhreeqc->GetSelectedOutputRowCount();

    for(size_t iRow = s_soFirstDataRow; iRow < nRows; ++iRow)
    {

        VAR countVar;
        countVar.type = TT_EMPTY;
        internalPhreeqc->GetSelectedOutputValue(iRow, countCol, &countVar);

        assert(countVar.type == TT_DOUBLE);

        std::vector<std::string> phaseNamesRow;

        {
            VAR strVar;
            VarInit(&strVar);
            internalPhreeqc->GetSelectedOutputValue(iRow, phasesCol, &strVar);
            assert(strVar.type == TT_STRING);

            const std::string allPhasesStr(strVar.sVal);

            VarClear(&strVar);

            boost::split(phaseNamesRow, allPhasesStr, boost::is_any_of("!"));
            phaseNamesRow.erase(std::remove(phaseNamesRow.begin(), phaseNamesRow.end(), ""),
                                phaseNamesRow.end());

            assert(phaseNamesRow.size() == static_cast<size_t>(countVar.dVal + 0.5));

            phaseNames.reserve( phaseNames.size() + phaseNamesRow.size() );
            phaseNames.insert( phaseNames.end(), phaseNamesRow.begin(), phaseNamesRow.end() );
        }

        {
            VAR strVar;
            VarInit(&strVar);
            internalPhreeqc->GetSelectedOutputValue(iRow, volCol, &strVar);
            assert(strVar.type == TT_STRING);

            const std::string allVolStr(strVar.sVal);

            VarClear(&strVar);

            std::vector<std::string> volStrList;
            boost::split(volStrList, allVolStr, boost::is_any_of("!"));
            volStrList.erase(std::remove(volStrList.begin(), volStrList.end(), ""),
                             volStrList.end());

            assert(volStrList.size() == static_cast<size_t>(countVar.dVal + 0.5));

            for(size_t iPhase = 0; iPhase < phaseNamesRow.size(); ++iPhase)
            {
                const std::string& phaseName = phaseNamesRow[iPhase];
                const std::string& volStr = volStrList[iPhase];
                const double vol = std::stod(volStr);
                molVols[phaseName] = vol;
            }
        }
    }

    std::sort(phaseNames.begin(), phaseNames.end());
    phaseNames.erase(std::unique(phaseNames.begin(), phaseNames.end()), phaseNames.end());
}


static void generateComposition(const IPPConfig::ElemSpeciesVec& allElems,
                                const std::vector<std::string>& allInputPhases,
                                const size_t iComposition, const double& conc,
                                const double& phaseConc, const std::string& compName,
                                std::stringstream& ss)
{
    ss << "SOLUTION " << std::to_string(iComposition) << " " << compName << std::endl;
    ss << "pH 7 charge\n";
    ss << "\tunits mol/L\n";
    for(const ElementSpecies& elemSpec : allElems)
    {
        PhreeqcInputFileGenerator::generateElementConcLine(elemSpec.elementName, conc, elemSpec.speciesName, ss);
    }


    ss << "EQUILIBRIUM_PHASES " << std::to_string(iComposition) << " " << compName << std::endl;
    for(const std::string& phaseName : allInputPhases)
    {
        ss << "\t" << phaseName << " 0.0 " << phaseConc << std::endl;
    }
}


static void generateAllInComposition(const IPPConfig::ElemSpeciesVec& allElems,
                                const std::vector<std::string>& allInputPhases,
                                const size_t iComposition,
                                std::stringstream& ss)
{
    static const double dummyConc = 1.0E-9;
    static const double dummyPhaseMole = 1.0E-3;
    generateComposition(allElems, allInputPhases, iComposition, dummyConc, dummyPhaseMole, "allElements", ss);
}

// generate dummy config for retrieving aux infos
// solution -> purpose
// 0                => possible phases + molar volumes
// 1 -> (nBC+1)     => boundary conditions
// nBC+1 -> nCells  => water only for empty cells
void generateDummyInputString(const IPPConfig& conf, const size_t nCells,
                              std::stringstream& ss)
{
    ss << "TITLE no name\n";


    // IMPORTANT: all in cell must be the first!
    // generate all-in cell:
    const IPPConfig::ElemSpeciesVec& allElems = conf.getAllElements();
    const std::vector<std::string>& allInputPhases = conf.getInputAllPhases();
    generateAllInComposition(allElems, allInputPhases, 0, ss);



    // generate boundary condition cells:
    std::map<ElementSpecies, double> allElemConc;
    for(const ElementSpecies& elemName : allElems)
    {
        allElemConc[elemName] = 0.0;
    }

    const std::vector<std::string>& allPhases = conf.getInputAllPhases();
    std::map<std::string, double> phaseNameToAmount;
    for(const std::string& phaseName : allPhases)
    {
        phaseNameToAmount[phaseName] = 0.0;
    }

    const ConfigBoundaryConditions::DiffusiveBCVec& diffBC = conf.boundaryConditions->diffusiveBC;
    typedef ConfigBoundaryConditions::DiffusiveBoundaryConditionData DiffusiveBoundaryConditionData;
    const size_t nBCs = diffBC.size();
    for(size_t iComposition = 0; iComposition < nBCs; ++iComposition)
    {
        const size_t iPhreeqcComposition = iComposition+1;
        const DiffusiveBoundaryConditionData& bcData = diffBC[iComposition];
        const Composition* comp = bcData.composition;
        if(comp)
        {
            const std::string& compName = comp->name;

            ss << "SOLUTION " << std::to_string(iPhreeqcComposition) << " " << compName << std::endl;
            // append solution body
            PhreeqcInputFileGenerator::generateSolutionBody(*comp, allElemConc, ss);

            ss << "EQUILIBRIUM_PHASES " << std::to_string(iComposition+1) << " " << compName << std::endl;
            PhreeqcInputFileGenerator::generateEqPhaseBody(comp->phases, *conf.dissPrecBehav, ss);

        }
        else
        {
            generateAllInComposition(allElems, allInputPhases, iPhreeqcComposition, ss);
        }
    }


    // at least one more water cell
    const std::size_t firstWaterCell = nBCs+1;
    assert(firstWaterCell  < nCells);

    // generate empty water cells
    // fill rest of cells with empty cells
    for(size_t iCell = firstWaterCell; iCell < nCells; ++iCell)
    {
        generateComposition(allElems, allInputPhases, iCell, 0.0, 0.0, "empty", ss);
    }

    PhreeqcInputFileGenerator::generateInertSoluteDatabaseExt(conf.results.diffEffInfos, ss);

    generateMolVolString(ss);

    PhreeqcInputFileGenerator::generateConvergenceHacksInitRun(ss);

}

typedef BoundaryConditions::BCVec BCVec;
typedef BoundaryConditionDomain BCDomain;
void addBCData(const BCDomain& domain, const IPP::BoundaryConditionType& type, BCVec& resultBCVec)
{
    resultBCVec.push_back(AbstractBoundaryConditions::BoundaryConditionData());
    AbstractBoundaryConditions::BoundaryConditionData& outputBC = resultBCVec.back();

    outputBC.domain = domain;
    outputBC.type = type;
    outputBC.density = 0.0;
    outputBC.velocity = {{0.0, 0.0, 0.0}};
}

void setupBCData(const ConfigBoundaryConditions::DiffusiveBCVec& inputBCs,
                 const IPPConfig::Results::DiffEffInfos& passiveTracerInfos, const int nCells,
                 const std::vector<std::string>& compNames, const std::vector<double>& conc,
                 BoundaryConditions::ElementNameToBCVec& resultDiffusiveBCs)
{

    assert(resultDiffusiveBCs.empty());

    // insert empty BCs
    const std::size_t nSpecies = compNames.size();
    for(size_t iComp = 0; iComp < nSpecies; ++iComp)
    {
        const std::string& compName = compNames[iComp];

        // insert new BC
        IPPCheck::assertCheck(resultDiffusiveBCs.find(compName) == resultDiffusiveBCs.end());
        BCVec& resultBCVec = resultDiffusiveBCs[compName];


        // exclude passive tracers for eff diff calc
        if(std::find_if(passiveTracerInfos.begin(), passiveTracerInfos.end(),
                        [&](const IPPConfig::Results::DiffEffInfo& info) { return info.tracerName == compName;} )
                        != passiveTracerInfos.end())
        {
            continue;
        }

        for(size_t iBC = 0; iBC < inputBCs.size(); ++iBC)
        {
            const ConfigBoundaryConditions::DiffusiveBoundaryConditionData& diffInputBC = inputBCs[iBC];

            // special treatment for charge species in case of neumann type BC
            // charge species gives problems when having outflow flux BC
            if(diffInputBC.type == BCT_VelocityOutflow && compName == "Charge")
            {
                addBCData(diffInputBC.domain, BCT_DensityDirichlet, resultBCVec);
            }
            else
            {
                addBCData(diffInputBC.domain, diffInputBC.type, resultBCVec);
            }
        }

    }


    assert(resultDiffusiveBCs.size() == nSpecies);
    for(size_t iComp = 0; iComp < nSpecies; ++iComp)
    {
        const std::string& compName = compNames[iComp];

        // exclude passive tracers for eff diff calc
        if(std::find_if(passiveTracerInfos.begin(), passiveTracerInfos.end(),
                        [&](const IPPConfig::Results::DiffEffInfo& info) { return info.tracerName == compName;} )
                        != passiveTracerInfos.end())
        {
            continue;
        }

        IPPCheck::assertCheck(resultDiffusiveBCs.find(compName) != resultDiffusiveBCs.end());
        typedef BoundaryConditions::BCVec BCVec;
        BCVec& resultBCVec = resultDiffusiveBCs[compName];

        assert(resultBCVec.size() == inputBCs.size());
        const size_t cellOffset = iComp * nCells;

        for(size_t iBC = 0; iBC < inputBCs.size(); ++iBC)
        {
            const ConfigBoundaryConditions::DiffusiveBoundaryConditionData& inputData = inputBCs[iBC];

            if(inputData.composition)
            {
                AbstractBoundaryConditions::BoundaryConditionData& outputData = resultBCVec.at(iBC);
                const size_t resultIndex = (iBC+1) + cellOffset;
                outputData.density = conc.at(resultIndex);
            }

        }
    }
}

static void setupInertPhaseComposition(const size_t nCells, const size_t nSpecies,
                                  const std::vector<double>& conc,
                                  std::vector<double>& inertComposition)
{
    IPPCheck::assertCheck(nCells > 0);

    const size_t cellOffset = nCells - 1;

    inertComposition.resize(nSpecies);
    for(size_t iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    {
        const size_t srcIndex = iSpecies * nCells + cellOffset;
        inertComposition[iSpecies] = conc[srcIndex];
    }
}

void findPhaseCoefs(PhreeqcRM& phreeqc, const std::vector<std::string>& phaseNames,
                    const std::map<std::string, size_t>& elemNameToID,
                    PhaseNameToInfos& phaseToInfo)
{
    const size_t countCol = 0;
    const size_t elementsCol = 1;
    const size_t coefsCol = 2;
    const size_t gfwCol = 3;


    for(size_t iPhase = 0; iPhase < phaseNames.size(); ++iPhase)
    {
        const std::string& phaseName = phaseNames[iPhase];
        phreeqc.SetCurrentSelectedOutputUserNumber(s_selectedOutPhaseFomulas + iPhase);

        PhaseInfo& info = phaseToInfo.at(phaseName);
        PhaseInfo::ElemIdToStoich& stoich = info.stoich;
        stoich.resize(elemNameToID.size(), 0);

        const std::vector<IPhreeqcPhast*>& internalPhreeqcVec = phreeqc.GetWorkers();
        IPhreeqcPhast* internalPhreeqc = internalPhreeqcVec.front();

        VAR countVar;
        countVar.type = TT_EMPTY;
        internalPhreeqc->GetSelectedOutputValue(s_soFirstDataRow, countCol, &countVar);


        const VAR_TYPE type = countVar.type;
        IPPCheck::assertCheck(type == TT_DOUBLE);

        const size_t nElems = countVar.dVal;
        (void)nElems;

        std::vector<std::string> elemStrList;
        {
            VAR strVar;
            VarInit(&strVar);
            VRESULT status = internalPhreeqc->GetSelectedOutputValue(s_soFirstDataRow, elementsCol, &strVar);
            if(status != VR_OK)
            {
                std::cerr << status << std::endl;
                throw std::runtime_error(std::to_string(status));
            }

            assert(strVar.type == TT_STRING);

            const std::string allElemStr(strVar.sVal);

            VarClear(&strVar);


            boost::split(elemStrList, allElemStr, boost::is_any_of("!"));
            elemStrList.erase(std::remove(elemStrList.begin(), elemStrList.end(), ""),
                             elemStrList.end());

        }

        std::vector<std::string> coefsStrList;
        {
            VAR strVar;
            VarInit(&strVar);
            VRESULT status = internalPhreeqc->GetSelectedOutputValue(s_soFirstDataRow, coefsCol, &strVar);
            if(status != VR_OK)
            {
                std::cerr << status << std::endl;
                throw std::runtime_error(std::to_string(status));
            }

            assert(strVar.type == TT_STRING);

            const std::string allCoefsStr(strVar.sVal);

            VarClear(&strVar);

            boost::split(coefsStrList, allCoefsStr, boost::is_any_of("!"));
            coefsStrList.erase(std::remove(coefsStrList.begin(), coefsStrList.end(), ""),
                             coefsStrList.end());

        }

        IPPCheck::assertCheck(coefsStrList.size() == elemStrList.size());

        for(size_t iElem = 0; iElem < elemStrList.size(); ++iElem)
        {
            const std::string& elemName = elemStrList[iElem];
            const std::string& coefStr = coefsStrList[iElem];
            const double coef = std::stod(coefStr);

            const size_t elemId = elemNameToID.at(elemName);
            stoich[elemId] = coef;
        }

        {
            VAR floatVar;
            VarInit(&floatVar);

            VRESULT status = internalPhreeqc->GetSelectedOutputValue(s_soFirstDataRow, gfwCol, &floatVar);
            if(status != VR_OK)
            {
                std::cerr << status << std::endl;
                throw std::runtime_error(std::to_string(status));
            }

            assert(floatVar.type == TT_DOUBLE);

            const double val = floatVar.dVal;
            info.gfw = val;

            IPPCheck::assertCheck(info.gfw > 0.0, "molar mass is zero for: " + phaseName);

            VarClear(&floatVar);
        }
    }
}

void writePhreeqcInput(const std::string& inputStr, const IPPConfig& conf, const std::string& preRunFilePrefix)
{
    const std::string fileName = preRunFilePrefix + ".pqi";
    std::ofstream f(fileName);
    IPPCheck::assertCheck(f.is_open(), "Could not open file: " + fileName);
    f << "DATABASE " << conf.dataBaseName << std::endl;
    f << inputStr;
    f.close();
}

void getAuxDataFromPhreeqC(const IPPConfig& conf,
                           const boost::filesystem::path& outDir,
                           const bool calcH2OSeperately,
                           std::vector<std::string>& compNames,
                           std::vector<std::string>& phaseNames,
                           PhaseNameToInfos& phaseToInfos,
                           PhaseIdPerElement& phaseIDsPerElement,
                           BoundaryConditions::ElementNameToBCVec& resultBCs,
                           std::vector<double>& inertComposition)
{
    boost::mpi::communicator world;
    pcout << "starting prerun. this is rank: " << world.rank() << std::endl;

    const ConfigBoundaryConditions::DiffusiveBCVec& inputBCs = conf.boundaryConditions->diffusiveBC;
    const size_t nDiffBC = inputBCs.size();

    // we need at least allIn-Cell + nBC Cells + empty water cell -> nBCs + 2
    const int nCells = std::max((int)nDiffBC + 2, world.size());

    PhreeqcRM phreeqc(nCells, MPI_COMM_WORLD);


    std::vector<double> conc;


    if(world.rank() == 0)
    {
        IRM_RESULT status;

        phreeqc.SetErrorHandlerMode(1);

        const boost::filesystem::path outPath = outDir / "prerun";
        const std::string preRunFilePrefix = outPath.string();
        phreeqc.SetFilePrefix(preRunFilePrefix);
        phreeqc.OpenFiles();


        phreeqc.SetComponentH2O(calcH2OSeperately);
        phreeqc.SetSpeciesSaveOn(false);

        phreeqc.SetUnitsSolution(PhreeqcConstants::s_unitSolution);
        phreeqc.SetUnitsPPassemblage(PhreeqcConstants::s_unitPPassemblage);

        // Set printing of chemistry file
#ifdef DEBUG_OUTPUT
        status = phreeqc.SetPrintChemistryOn(true, true, true); // workers, initial_phreeqc, utility
#else
        status = phreeqc.SetPrintChemistryOn(false, false, false); // workers, initial_phreeqc, utility
#endif
        IPPCheck::assertCheck(status == IRM_OK);

        status = phreeqc.LoadDatabase(conf.getDataBaseRelativePath());
        IPPCheck::assertCheck(status == IRM_OK);

        // generate dummy config for retrieving aux infos
        // solution -> purpose
        // 0                => possible phases + molar volumes
        // 1 -> (nBC+1)     => boundary conditions
        // nBC+1 -> nCells  => water only for empty cells
        std::stringstream inputStrStrm;
        generateDummyInputString(conf, nCells, inputStrStrm);
        const std::string inputStr = inputStrStrm.str();


        // TODO: write in debug build only?
        writePhreeqcInput(inputStr, conf, preRunFilePrefix);

        // disable rebalancing because we are using internal phreeqc data
        phreeqc.SetRebalanceByCell(false);
        phreeqc.SetRebalanceFraction(0.0);

        bool workers = true;
        bool initial_phreeqc = true;
        bool utility = true;
        phreeqc.RunString(workers, initial_phreeqc, utility, inputStr);


        // Get component information
        const std::size_t nComp = phreeqc.FindComponents();
        // PhreeqC behaves strange... possible side effect above, but we ignore return value
        (void) nComp;

        compNames = phreeqc.GetComponents();


        assert(nComp == compNames.size());

        // retrieve phase names + molar volumes
        std::map<std::string, double> phaseNameToMolVol;
        findElementPhaseNamesAndMolVols(phreeqc, phaseNames, phaseNameToMolVol);

        typedef IsInContainer<std::set<std::string> > IsInStringSet;
        phaseNames.erase(std::remove_if(phaseNames.begin(), phaseNames.end(),
                                        IsInStringSet(conf.blackListedPhases)),
                         phaseNames.end());
        std::sort(phaseNames.begin(), phaseNames.end());

        const PhaseNameToInfos& inputPhaseInfos = conf.dataBaseExt.phaseNameToInfos;


        for(const std::string& phaseName : phaseNames)
        {
            PhaseInfo& output = phaseToInfos[phaseName];

            double molVol = -1.0;

            const auto it = inputPhaseInfos.find(phaseName);
            if(it != inputPhaseInfos.end())
            {
                if(it->second.molarVolume > 0.0)
                {
                    molVol = it->second.molarVolume;
                }
                else
                {
                    // invalid value found in DB extension
                    molVol = phaseNameToMolVol.at(phaseName);
                }

                output.isPermeable = it->second.isPermeable;
            }
            else
            {                
                // not found in DB extension
                molVol = phaseNameToMolVol.at(phaseName);
                output.isPermeable = false;
            }

            if(molVol <= 0.0)
            {
                const double minMolVol = 0.0;
                std::cout << "WARNING: No valid molar volume given for phase: " << phaseName
                          << "\tassuming " << std::to_string(minMolVol) << " cm3/mol" << std::endl;
                molVol = minMolVol;
            }

            output.molarVolume = molVol;
        }

        // retrieve BC concentrations
        phreeqc.GetConcentrations(conc);
        IPPCheck::assertCheck(conc.empty() == false);



        std::stringstream coefStrm;
        generateFormulaeString(phaseNames, coefStrm);
        const std::string coefStr = coefStrm.str();

        // TODO: write in debug build only?
        const std::string completeStr = inputStr + "\n" + coefStr;
        writePhreeqcInput(completeStr, conf, preRunFilePrefix + "_coef");

        status = phreeqc.RunString(workers, initial_phreeqc, utility, coefStr);
        IPPCheck::assertCheck(status == IRM_OK);

        status = phreeqc.RunCells();
        IPPCheck::assertCheck(status == IRM_OK);

        std::map<std::string, size_t> elemNameToID;
        for(size_t iElem = 0; iElem < compNames.size(); ++iElem)
        {
            const std::string& elemName = compNames[iElem];
            elemNameToID[elemName] = iElem;
        }
        findPhaseCoefs(phreeqc, phaseNames, elemNameToID, phaseToInfos);


        phreeqc.CloseFiles();

        status = phreeqc.MpiWorkerBreak();
        IPPCheck::assertCheck(status == IRM_OK);
    }
    else
    {
        phreeqc.MpiWorker();
    }


    // retrieve phase names for worker
    MpiTools::syncFromRoot(world, phaseNames);

    boost::mpi::broadcast(world, phaseToInfos, 0);

    // retrieve species names for worker
    MpiTools::syncFromRoot(world, compNames);

    // retrieve BC conc
    MpiTools::syncFromRoot(world, conc);


    phaseIDsPerElement.init(compNames.size());
    for(const auto& infoNamePair : phaseToInfos)
    {
        const std::string& phaseName = infoNamePair.first;

        // assuming phaseNames as sorted
        assert(std::is_sorted(phaseNames.begin(), phaseNames.end()));
        const std::vector<std::string>::const_iterator it =
                std::find(phaseNames.begin(), phaseNames.end(), phaseName);
        assert(it != phaseNames.end());

        const size_t phaseID = std::distance(phaseNames.cbegin(), it);

        const PhaseInfo& info = infoNamePair.second;
        const PhaseInfo::ElemIdToStoich& stoich = info.stoich;

        for(size_t iElement = 0; iElement < stoich.size(); ++iElement)
        {
            if(stoich[iElement] > 0)
            {
                phaseIDsPerElement.add(iElement, phaseID);
            }
        }


    }

    setupBCData(inputBCs, conf.results.diffEffInfos, nCells, compNames, conc, resultBCs);
    setupInertPhaseComposition(nCells, compNames.size(), conc, inertComposition);

    pcout << "finished prerun. this is rank: " << world.rank() << std::endl;

}

}

}
