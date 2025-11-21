#include "phreeqcreactionmodulefactory.h"

#include "phreeqcreactionmodule.h"
#include "dummyreactionmodule.h"

namespace IPP
{

ReactionModule *PhreeqcReactionModuleFactory::create(const boost::mpi::communicator& comm,
                                                     const bool isChemistryEnabled,
                                                     const ReactionModuleConfig &reacConf,
                                                     const size_t nxyz)
{
    ReactionModule* reacModule = nullptr;

    if(isChemistryEnabled)
    {
        reacModule = new PhreeqcReactionModule(nxyz, comm, reacConf);
    }
    else
    {
        reacModule = new DummyReactionModule(nxyz, reacConf);
    }

    return reacModule;
}


}
