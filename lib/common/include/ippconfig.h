#ifndef IPPCONFIG_H
#define IPPCONFIG_H

#include <boost/mpi/communicator.hpp>

#include <array>
#include <boost/functional/hash.hpp>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include <phasenametoinfos.h>
#include <abstractflowfunctor.h>
#include <configboundaryconditions.h>
#include <abstractdisspreciponlyinfo.h>
#include "abstractsicalcfwd.h"
#include "abstractporositycalcfwd.h"
#include "abstractmultiscalediffusioncalcfwd.h"
#include "abstractdissolveonlycalcfwd.h"
#include "resultstowrite.h"
#include "latticeboltzmanncollisiontype.h"
#include "auxresultprocessing.h"

namespace IPP
{

struct ElementSpecies
{
    ElementSpecies()
        : elementName("N/A"),
          speciesName("N/A")
    {

    }

    ElementSpecies(const std::string& elementName)
        : elementName(elementName),
          speciesName()
    {

    }

    bool operator==(const ElementSpecies& other) const
    {
        return elementName == other.elementName &&
                speciesName == other.speciesName;
    }

    bool operator<(const ElementSpecies& other) const
    {
        const bool isElemEqual = elementName == other.elementName;
        if(isElemEqual == false)
        {
            return elementName < other.elementName;
        }
        else
        {
            return speciesName < other.speciesName;
        }
    }

    std::string elementName;
    std::string speciesName;
};

} // end of namespace



namespace std
{
template <>
struct hash<IPP::ElementSpecies>
{
    typedef IPP::ElementSpecies argument_type;
    typedef std::size_t result_type;

    result_type operator()(const IPP::ElementSpecies & t) const
    {
        std::size_t val { 0 };
        boost::hash_combine(val, t.elementName);
        boost::hash_combine(val, t.speciesName);
        return val;
    }
};
}


namespace IPP
{


class Composition
{
public:
    Composition()
        : temp(25.0),
          pressure(1),
          pH(7.0, "charge"),
          pe(4.0, ""),
          inertVolFraction(0.0)
    {

    }


    struct ElementConcentration
    {
        ElementConcentration()
            : elemSpecies(),
              conc(-1.0)
        {

        }

        ElementConcentration(const ElementSpecies& elemSpecies,
                             const double& conc)
            : elemSpecies(elemSpecies)
            , conc(conc)
        {

        }

        ElementSpecies elemSpecies;
        double conc;
    };

    struct ScalarAndString
    {
        ScalarAndString()
            : value(-1),
              extension()
        {

        }

        ScalarAndString(const double& value,
                        const std::string& extension)
            : value(value),
              extension(extension)
        {

        }

        double value;
        std::string extension;
    };

    struct PhaseDefinition
    {
        PhaseDefinition()
            : name("N/A"),
              amount(-1.0)
        {

        }

        PhaseDefinition(const std::string& name,
                        const double& amount)
            : name(name)
            , amount(amount)
        {

        }

        std::string name;
        double amount;
    };

    std::string name;
    double temp;
    double pressure;
    ScalarAndString pH;
    ScalarAndString pe;
    std::string units;
    std::vector<ElementConcentration> elemConcVec;

    typedef std::vector<PhaseDefinition> PhaseDefVec;
    PhaseDefVec phases;
    double inertVolFraction;
};


struct Domain
{
    explicit Domain(const Composition* composition)
        : composition(composition)
    {

    }

    const Composition* composition;
    std::vector<int> cells;
};

class IPPConfig
{
public:

    struct DataBaseExtension
    {
        PhaseNameToInfos phaseNameToInfos;
    };

    enum TransportOnlyFlag
    {
        TOF_FullChemistry,
        TOF_FakeChem,
        TOF_InitChemOnly
    };


    IPPConfig(int argc, char* argv[]);
    void init();

    const std::string getDataBaseRelativePath() const;

    typedef std::vector<ElementSpecies> ElemSpeciesVec;
    const ElemSpeciesVec& getAllElements() const;

    const std::vector<std::string>& getInputAllPhases() const;

    int argc;
    char** argv;
    const boost::mpi::communicator globalComm;

    bool verbose;

    size_t maxIter;
    double maxTime;
    size_t outputFrequency;
    size_t outletFluxFrequency;
    double checkpointingTime;
    size_t checkpointingFreq;
    double diffusionCoefReference;
    double tauRef;
    double porosRef;
    double porosLow;
    LatticeBoltzmannCollisionType latticeBoltzmannCollType;
    std::string outputDir;
    bool isChemistryEnabled;
    TransportOnlyFlag transportOnlyFlag;
    bool isDegradedDiffCalc;
    double chemistryTolerance;
    bool write_chem;
    std::string dataBasePath;
    std::string dataBaseName;

    size_t nx;
    size_t ny;
    size_t nz;

    bool is3DSimulation;

    double spatialResolution;
    double distTransPorosThresh;
    IPPVector3D offset;

    std::vector<IPPBox3DInt> chemDisabledBoxes;
    IPPBox3DInt transportSyncEnabledBox;
    std::vector<IPPBox3DInt> decomp;

    AbstractFlowFunctorPtr flowFunc;
    ConfigBoundaryConditionsPtr boundaryConditions;
    AbstractSICalcPtr saturationIndexCalc;
    AbstractDissolveOnlyCalcPtr dissolveOnlyFunc;

    AbstractPorosityCalcPtr porosCalc;
    AbstractMultiScaleDiffusionCalcPtr multiScaleDiffusionCalc;

    DataBaseExtension dataBaseExt;

    struct ChemistryPrediction
    {
        ChemistryPrediction()
            : enabled(false)
            , filename()
            , significantPorosDiff()
            , forceRecalcFreq()
        {

        }

        bool enabled;
        std::string filename;
        double significantPorosDiff;
        size_t forceRecalcFreq;
    };

    ChemistryPrediction chemPredict;

    struct SolidSolutionDefinition
    {
        std::string name;
        std::vector<std::string> phaseNames;
    };

    std::vector<SolidSolutionDefinition> solidSolutionVec;
    std::set<std::string> solidSolutionPhases;

    std::vector<Composition> compositions;
    std::vector<Domain> domains;

    std::set<std::string> whiteListedPhases;
    std::set<std::string> blackListedPhases;
    std::set<std::string> waterLimitedPhases;
    std::set<std::string> kineticPhases;


    std::unique_ptr<AbstractDissPrecipOnlyInfo> dissPrecBehav;


    double fluxConvergenceTolLow;
    double fluxConvergenceTolHigh;
    double fluxConvergenceTolFactor;


    struct Results
    {
        ResultsToWrite resultsToWrite;


        struct DiffEffInfo
        {
            DiffEffInfo(size_t dim, const std::string& tracerName)
                : dim(dim)
                , tracerName(tracerName)
            {

            }

            size_t dim;
            std::string tracerName;
        };

        typedef std::vector<DiffEffInfo> DiffEffInfos;
        DiffEffInfos diffEffInfos;

        struct NameCommand
        {
            NameCommand(const std::string& name, const std::string& command)
                : name(name),
                  command(command)
            {

            }

            const std::string name;
            const std::string command;
        };

        std::vector<NameCommand> auxData;
    };


    Results results;


    std::shared_ptr<AuxResultProcessing> resultProcessing;

private:
    ElemSpeciesVec m_allElements;
    std::vector<std::string> m_allInputPhases;

};

typedef std::shared_ptr<IPPConfig> IPPConfigPtr;

} // end of namespace IPP

#endif // IPPCONFIG_H

