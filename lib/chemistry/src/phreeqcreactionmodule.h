#ifndef PHREEQCREACTIONMODULE_H
#define PHREEQCREACTIONMODULE_H

#include <boost/mpi/communicator.hpp>
#include <boost/filesystem/path.hpp>

#include "reactionmodule.h"

#include "localaccessphreeqcrm.h"
#include "reactionexchangedata.h"
#include "arraydimensionconvert.h"
#include "phreeqcinitialdata.h"
#include "simpleporosityinfosfwd.h"
#include "nodeinfos.h"
#include "simpleporosityinfosfwd.h"

namespace IPP
{

class EnablePhreeqcCells;
class DisablePhreeqcCells;

class PhreeqcReactionModuleData;
class NonLocalOperations;
struct PhreeqcGlobalInitData;

class PhreeqcReactionModule : public ReactionModule
{
public:
    PhreeqcReactionModule(const size_t nCellsGlobal,
                          const boost::mpi::communicator &comm,
                          const ReactionModuleConfig& conf);
    virtual ~PhreeqcReactionModule() override;

    virtual ConfigurePhreeqcCellsActivity* getConfigCellsActivity() override;
    virtual AbstractGeometrySync* getGeomSync() override;

    virtual void saveCheckpoint() override;
    virtual void loadCheckpoint(const size_t iteration) override;

    virtual void errorDump() override;

    virtual void run() override;

    virtual void updateComponentConcentrations() override;

    virtual void updatePreCondition() override;
    virtual void updatePostCondition() override;

private:

    struct PrecipData
    {
        std::vector<double> localPrecipMol;
        std::vector<double> solidFractions;
    };


    // init functions
    void init(const size_t nCellsGlobal, const ReactionModuleConfig& conf);
    void preInit(const ReactionModuleConfig& conf);
    void initCells(const size_t nCellsGlobal,
                   const ReactionModuleConfig& conf,
                   const PhreeqcGlobalInitData& initData);
    void initConcentrations();
    void initRoot(const size_t nCellsGlobal,
                  const ReactionModuleConfig& conf,
                  PhreeqcGlobalInitData& initData);
    void initRootRun(const ReactionModuleConfig& conf, PhreeqcGlobalInitData& initData);
    void disableLowPorosCells(const std::vector<double>& porosVec,
                              const double& lowPorosThresh);
    void pushInitData(PhreeqcGlobalInitData& initData);
    void initPrecipMoles(const ReactionModuleConfig &conf,
                         const std::vector<size_t> &cellToDomain);
    void convergePorosity();
    void correctByPorosChange(std::vector<double>& concVec);
    void collectPorosityConvergenceData();


    // iterating functions
    void configCellsActivity();

    void runHomo();
    void runHetero();

    void runCells();
    void swapToOld();

    void collectDataAll();
    void collectData(PrecipData& tmpData);

    void calcCapillaryPorosity(const std::vector<double>& precipMol,
                               const std::vector<double>& solidFractions,
                               std::vector<SimplePorosityInfosPtr>& porosInfos);

    void getUpdatedConcentrations(std::vector<double>& postReacConc);
    void treatNegativePorosity();

    void retrieveConcentrations();
    //    void retrieveTotals();
    void retrieveSaturationIndices(std::vector<double>& satIndices);
    void retrieveMonomersConc();


    void initMapping(const std::vector<int>& grid2chem,
                     const std::vector<int>& ic);
    void setTime();

    // checkpoint functions
    void saveDump(const boost::filesystem::path& iterationPath);
    void loadDump(const boost::filesystem::path& iterationPath);

    void savePorosity(const boost::filesystem::path iterationPath) const;
    void loadPorosity(const boost::filesystem::path& iterationPath);

    void saveOldPorosity(const boost::filesystem::path iterationPath) const;
    void loadOldPorosity(const boost::filesystem::path &iterationPath);

    void saveNodeInfos(const boost::filesystem::path& iterationPath) const;
    void loadNodeInfos(const boost::filesystem::path& iterationPath);

    void saveAuxData(const boost::filesystem::path &iterationPath) const;
    void loadAuxData(const boost::filesystem::path &iterationPath);


    const boost::mpi::communicator &m_comm;

    bool m_isPrintChemistryOn;

    LocalAccessPhreeqcRM m_phreeqc;
    std::shared_ptr<SimulationExchangeData> m_simData;
    ReactionExchangeData m_reacData;

    boost::filesystem::path m_chemDir;
    boost::filesystem::path m_dumpDir;

    NonLocalOperations& m_nonLocalOperations;

    PhreeqcReactionModuleData& m_phreeqcData;


    std::vector<double> m_currTargetSImap;
    std::vector<unsigned char> m_currDissolveOnlyCells;


    NodeInfos m_nodeInfos;


    std::unique_ptr<ConfigurePhreeqcCellsActivity> m_configActivity;
    std::unique_ptr<EnablePhreeqcCells> m_enableCellsFunc;
    std::unique_ptr<DisablePhreeqcCells> m_disableCellsFunc;
    std::vector<char> m_enabledCells;


    std::vector<double> m_transportDiff;
    std::vector<double> m_oldPorosPhysical;
    std::vector<double> m_preReacConc;
    std::vector<double> m_oldPrecipAmount;

    PhreeqcInitialData m_initialData;

    const ArrayDimensionConvert m_indexConvLocal;

    std::unique_ptr<AbstractGeometrySync> m_geomSync;

    bool m_verbose;

    void collectSteadyState();
};

} // end of namespace IPP

#endif // PHREEQCREACTIONMODULE_H
