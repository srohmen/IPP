#ifndef DUMMYREACTIONMODULE_H
#define DUMMYREACTIONMODULE_H

#include <memory>

#include "reactionmodule.h"
#include "phasenametoinfos.h"

namespace IPP
{
class SimulationExchangeData;

class DummyReactionModule : public ReactionModule
{
public:
    DummyReactionModule(const size_t nCells, const ReactionModuleConfig& conf);

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
    void init(const ReactionModuleConfig& conf);

    size_t m_nxyz;
    size_t m_nTracer;
    std::shared_ptr<SimulationExchangeData> m_simData;

    std::vector<std::string> phaseNamesDummy;
    PhaseNameToInfos phaseNameToInfosDummy;
    std::vector<std::string> compNamesDummy;

};

}

#endif // DUMMYREACTIONMODULE_H
