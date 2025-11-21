#ifndef REACTIONMODULE_H
#define REACTIONMODULE_H

#include <vector>
#include <cstddef>

namespace IPP
{

class ReactionModuleConfig;
class ConfigurePhreeqcCellsActivity;
class AbstractGeometrySync;

class ReactionModule
{
public:
    ReactionModule()
    {

    }

    virtual ~ReactionModule()
    {

    }

    virtual ConfigurePhreeqcCellsActivity* getConfigCellsActivity() = 0;
    virtual AbstractGeometrySync* getGeomSync() = 0;

    virtual void saveCheckpoint() = 0;
    virtual void loadCheckpoint(const size_t iteration) = 0;

    virtual void errorDump() = 0;

    virtual void run() = 0;

    virtual void updateComponentConcentrations() = 0;

    virtual void updatePreCondition() = 0;
    virtual void updatePostCondition() = 0;

};

} // end of namespace LBGeoChem

#endif // REACTIONMODULE_H
