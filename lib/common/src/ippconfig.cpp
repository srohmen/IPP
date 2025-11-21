#include "ippconfig.h"
#include <unordered_set>

#include "simpleporositycalc.h"
#include "thresholdmultiscalediffusioncalc.h"

namespace IPP
{

struct DisabledFlow : public AbstractFlowFunctor
{
    virtual bool isEnabled() const
    {
        return false;
    }

    virtual void getDefault(double& /*rho*/, IPPVector3D& /*u*/) const override
    {

    }

    virtual void operator()(int /*iX*/, int /*iY*/, double& /*rho*/,
                            IPPVector2D& /*u*/) const override
    {

    }

    virtual void operator()(int /*iX*/, int /*iY*/,
                            int /*iZ*/, double& /*rho*/,
                            IPPVector3D& /*u*/) const override
    {

    }

};


IPPConfig::IPPConfig(int argc, char* argv[])
    : argc(argc), argv(argv),
      globalComm(),
      verbose(false),
      maxIter(),
      maxTime(),
      outputFrequency(1),
      outletFluxFrequency(0),
      checkpointingTime(-1.0),
      checkpointingFreq(-1),
      diffusionCoefReference(1.0E-9),
      tauRef(1.0),
      porosRef(1.0),
      porosLow(1.0),
      latticeBoltzmannCollType(LBCT_PTRT),
      isChemistryEnabled(true),
      transportOnlyFlag(TOF_FullChemistry),
      isDegradedDiffCalc(false),
      chemistryTolerance(0.0),
      #ifdef DEBUG_OUTPUT
      write_chem(true),
      #else
      write_chem(false),
      #endif
      spatialResolution(-1.0),
      distTransPorosThresh(0.5),
      offset({{0.0, 0.0, 0.0}}),
      flowFunc(new DisabledFlow),
      boundaryConditions(nullptr),
      saturationIndexCalc(nullptr),
      dissolveOnlyFunc(nullptr),
      multiScaleDiffusionCalc(nullptr)
    , fluxConvergenceTolLow(0.0)
    , fluxConvergenceTolHigh(0.0)
    , fluxConvergenceTolFactor(2.0)
{

}

void IPPConfig::init()
{
    // check sanity of input

    if(tauRef <= 0.5)
    {
        throw std::runtime_error("choose relaxation parameter tauRef higher than 0.5");
    }

    if(isChemistryEnabled)
    {
        if(results.diffEffInfos.empty() == false)
        {
            throw std::runtime_error("chemistry enabled and inert diffusive tracer defined at the same time");
        }
    }
    else
    {
        if(results.diffEffInfos.empty())
        {
            throw std::runtime_error("neither chemistry enabled nor inert diffusive tracer defined");
        }

        if(boundaryConditions->diffusiveBC.empty() == false)
        {
            throw std::runtime_error("for chemistry disabled model chemical boundary conditions must not be defined but diffusive BC is defined");
        }

        if(boundaryConditions->advectiveBC.empty() == false)
        {
            throw std::runtime_error("for chemistry disabled model chemical boundary conditions must not be defined but advective BC is defined");
        }
    }


    // prepare all element set
    // collect all input phases
    std::unordered_set<ElementSpecies> elemSet;
    std::set<std::string> phasesSet = whiteListedPhases;

    for(const Composition& comp : compositions)
    {
        for(const Composition::ElementConcentration& elemConc : comp.elemConcVec)
        {
            elemSet.insert(elemConc.elemSpecies);
        }

        for(const Composition::PhaseDefinition& phase : comp.phases)
        {
            phasesSet.insert(phase.name);
        }
    }

    for(const Results::DiffEffInfo& info : results.diffEffInfos)
    {
        elemSet.insert(info.tracerName);
    }


    for(const std::string& phaseName : blackListedPhases)
    {
        phasesSet.erase(phaseName);
    }



    // copy to vectors
    for(const ElementSpecies& elem : elemSet)
    {
        m_allElements.push_back(elem);
    }

    std::sort(m_allElements.begin(), m_allElements.end());

    for(const std::string& phaseName : phasesSet)
    {
        m_allInputPhases.push_back(phaseName);
    }


    for(const SolidSolutionDefinition& def : solidSolutionVec)
    {
        for(const std::string& phaseName : def.phaseNames)
        {
            solidSolutionPhases.insert(phaseName);
        }
    }



    if(porosCalc == nullptr)
    {
        std::cout << "WARNING: No porosity model defined. Fallback to simple porosity calc model" << std::endl;
        porosCalc.reset(new SimplePorosityCalc);
    }

    if(multiScaleDiffusionCalc == nullptr)
    {
        const double fallBackDiff = 1.0E-9;
        std::cout << "WARNING: No diffusion calc defined. Fallback to free water diffusion coeffcient everywhere: "
                  << fallBackDiff << std::endl;
        multiScaleDiffusionCalc.reset(new ThresholdMultiScaleDiffusionCalc(0.0, fallBackDiff));
    }




    for(const Domain& domain : domains)
    {
        const double inertVolFrac = domain.composition->inertVolFraction;
        if(inertVolFrac < 1.0)
        {
            const double nonInertVolFrac = 1.0 - inertVolFrac;
            if(nonInertVolFrac < porosLow)
            {
                throw std::runtime_error("non inert volume fraction is lower than porosLow in composition: " + domain.composition->name);
            }
        }
    }

    for(const IPPBox3DInt& box : chemDisabledBoxes)
    {
        const bool isInters = isIntersecting(transportSyncEnabledBox, box);
        if(isInters)
        {
            std::stringstream ss;
            ss << "chemistry disabled box intersects with forced enabled box: " << box << std::endl;
            throw std::runtime_error(ss.str());
        }
    }
}

const std::string IPPConfig::getDataBaseRelativePath() const
{
    return dataBasePath + "/" + dataBaseName;
}

const IPPConfig::ElemSpeciesVec& IPPConfig::getAllElements() const
{
    return m_allElements;
}

const std::vector<std::string>& IPPConfig::getInputAllPhases() const
{
    return m_allInputPhases;
}

} // end of namespace LBGeoChem
