#ifndef PHREEQCREACTIONMODULEFACTORY_H
#define PHREEQCREACTIONMODULEFACTORY_H

#include <boost/mpi/communicator.hpp>


namespace IPP
{

class ReactionModule;
class ReactionModuleConfig;

namespace PhreeqcReactionModuleFactory
{
    ReactionModule *create(const boost::mpi::communicator& comm,
                           const bool isChemistryEnabled,
                           const ReactionModuleConfig& reacConf,
                           const size_t nxyz);
}

}

#endif // PHREEQCREACTIONMODULEFACTORY_H
