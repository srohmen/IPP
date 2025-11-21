#include "initporosity.h"

#include "phreeqcconstants.h"
#include "abstractporositycalc.h"
#include "abstractmultiscalediffusioncalc.h"
#include "porositycalcutils.h"
#include "ippexception.h"

#include "ippconfig.h"
#include "simulationexchangedata.h"
#include "phreeqccalcporosity.h"
#include "phreeqccalcdiffcoefs.h"
#include "serialtomultidatasync.h"
#include "geometrytools.h"
#include "cellindexcorrection.h"
#include "simpleporosityinfos.h"
#include "geometrytools.h"

namespace IPP
{

namespace InitPorosity
{


void extractPhasesMolarVolumes(const PhaseNameToInfos& phaseNameToInfos,
                               const std::vector<std::string>& phaseNames,
                               std::vector<double>& phaseMolarVolumesInLiter)
{
    // REV of phreeqc cells is 1L = 1000 cm^3
    const double refCellVolFactor = 0.001;

    phaseMolarVolumesInLiter.resize(phaseNames.size());
    for(size_t iPhase = 0 ; iPhase < phaseNames.size(); ++iPhase)
    {
        const std::string& phaseName = phaseNames[iPhase];
        const PhaseInfo& info = phaseNameToInfos.at(phaseName);
        phaseMolarVolumesInLiter[iPhase] = info.molarVolume * refCellVolFactor;
    }
}


void run(const std::vector<Domain>& domains,
         const CellIndexCorrection& indexCorr,
         const size_t nCells,
         const PhaseNameToInfos& phaseInfos,
         const std::vector<double>& phaseMolarVolumes,
         const std::vector<char>& isPermeablePhase,
         const std::vector<std::string>& phaseNames,
         const AbstractPorosityCalc& porosCalc,
         const AbstractMultiScaleDiffusionCalc& diffCalc,
         std::vector<double>& porosities,
         std::vector<double>& capillaryPorosities,
         std::vector<double>& diffCoefs,
         std::vector<double>& inertVolFracs)
{

    // prepare initial phase amount table
    typedef std::map<std::string, size_t> PhaseNameToIdMap;
    PhaseNameToIdMap phaseNameToID;
    for(size_t i = 0; i < phaseNames.size(); ++i)
    {
        const std::string& name = phaseNames[i];
        phaseNameToID[name] = i;
    }

    const size_t nDomains = domains.size();
    const size_t nPhases = phaseNames.size();
    std::vector<double> precipMol(nDomains * nPhases);
    std::vector<double> inertVolFrac(nDomains, 0.0);

    for(size_t iDomain = 0 ; iDomain < nDomains; ++iDomain)
    {
        const Domain& domain = domains[iDomain];

        assert(domain.composition);
        const Composition& comp = *(domain.composition);
        const Composition::PhaseDefVec& phaseDefs = comp.phases;
        for(size_t iPhase = 0; iPhase < phaseDefs.size(); ++iPhase)
        {
            const Composition::PhaseDefinition& phaseDef = phaseDefs[iPhase];

            PhaseNameToIdMap::const_iterator it = phaseNameToID.find(phaseDef.name);
            if(it == phaseNameToID.end())
            {
                throw std::runtime_error("Could not convert phase to ID: " + phaseDef.name);
            }

            const size_t phaseIndex = it->second;

            const size_t index = iDomain + phaseIndex * nDomains;
            precipMol[index] = phaseDef.amount;

        }

        inertVolFrac[iDomain] = comp.inertVolFraction;
    }

    GeometryTools::transpose(precipMol, nDomains, nPhases);




    PorosityCalc totalPorosCalcFunc(nDomains, phaseMolarVolumes,
                                    isPermeablePhase, inertVolFrac);
    std::vector<double> porosPerDomain(nDomains);
    std::vector<double> solidFractions(nDomains * nPhases);

    const av::bounds<2> bounds = {(ptrdiff_t)nDomains, (ptrdiff_t)nPhases};
    PorosityCalc::ConstCellPhasesView precipMolView(precipMol, bounds);
    PorosityCalc::CellPhasesView solidFracView(solidFractions, bounds);
    totalPorosCalcFunc.calcPorosity(precipMolView, solidFracView, porosPerDomain);

    for(size_t iDomain = 0; iDomain < porosPerDomain.size(); ++iDomain)
    {
        const double& poros = porosPerDomain[iDomain];
        if(poros < 0.0 || poros > 1.0)
        {
            throw std::runtime_error("domain has invalid porosity: "
                                     + domains[iDomain].composition->name
                                     + " ( " + std::to_string(poros) + " )");
        }
    }


    std::vector<double> diffCoefPerDomain;
    std::vector<double> capillaryPorosPerDomain;

    std::vector<SimplePorosityInfosPtr> porosInfos;
    PhreeqcCalcPorosity::calc(phaseInfos, solidFractions, precipMol, inertVolFrac, porosCalc, porosInfos);

    capillaryPorosPerDomain.resize(porosInfos.size());
    for(size_t iDomain = 0; iDomain < porosInfos.size(); ++iDomain)
    {
        const SimplePorosityInfosPtr& info = porosInfos[iDomain];
        capillaryPorosPerDomain[iDomain] = info->porosityCapillary;
    }

    PhreeqcCalcDiffCoefs::calc(diffCalc, porosInfos, diffCoefPerDomain);


    for(size_t iDomain = 0 ; iDomain < nDomains; ++iDomain)
    {
        const double& diffCoef = diffCoefPerDomain[iDomain];
        IPPCheck::assertCheck(diffCoef >= 0.0, "got negative diffusion coeff in domain: "
                              + domains[iDomain].composition->name);
    }



    porosities.resize(nCells, 1.0);
    capillaryPorosities.resize(nCells, 1.0);
    diffCoefs.resize(nCells, 1.0E-9);
    inertVolFracs.resize(nCells, 0.0);

    for(size_t iDomain = 0 ; iDomain < nDomains; ++iDomain)
    {
        const Domain& domain = domains[iDomain];

        const double& poros = porosPerDomain[iDomain];
        const double& porosCapillary = capillaryPorosPerDomain[iDomain];
        const double& diffCoef = diffCoefPerDomain[iDomain];
        const double& currInertVolFrac = domain.composition->inertVolFraction;


        const std::vector<int>& cells = domain.cells;

        for(size_t iDomainCell = 0; iDomainCell < cells.size(); ++iDomainCell)
        {
            size_t cellID = cells[iDomainCell];
            cellID = indexCorr.correctForBCs(cellID);

            assert(cellID < nCells);

            porosities[cellID] = poros;
            capillaryPorosities[cellID] = porosCapillary;
            diffCoefs[cellID] = diffCoef;
            inertVolFracs[cellID] = currInertVolFrac;
        }

    }

}


static void pushScalarData(const std::vector<double>& globalData,
                           const FieldDecomposition& decomp,
                           std::vector<double>& localData)
{
    SerialToMultiDataSync<double> sync(localData, decomp, 1);
    sync.pushFromRoot(globalData);
}


void syncData(const std::vector<double>& globalPorosities,
              const std::vector<double>& globalPorositiesCapillary,
              const std::vector<double>& globalDiffCoefs,
              const std::vector<double>& inertVolFrac,
              SimulationExchangeData& simData)
{
    pushScalarData(globalPorosities, *simData.getDecomp(),
                   simData.getPorosity());

    pushScalarData(globalPorositiesCapillary, *simData.getDecomp(),
                   simData.getCapillaryPorosity());

    pushScalarData(globalDiffCoefs, *simData.getDecomp(),
                   simData.getDiffusionCoefs());

    std::vector<double>& localInertVolFrac = simData.getInertVolFrac();
    pushScalarData(inertVolFrac, *simData.getDecomp(),
                   localInertVolFrac);


    std::vector<char>& inertCells = simData.getInertSolidCells();
    inertCells.resize(localInertVolFrac.size(), 0);
    for(size_t i = 0; i < localInertVolFrac.size(); ++i)
    {
        if(localInertVolFrac[i] >= 1.0)
        {
            inertCells[i] = 1;
        }
        else
        {
            inertCells[i] = 0;
        }
    }

}


}

}

