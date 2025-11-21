#ifndef TRANSPORTMODULEFACTORY_H
#define TRANSPORTMODULEFACTORY_H

#include "transportmoduleconfigfwd.h"

namespace IPP
{

class TransportModule;

namespace TransportModuleFactory
{
    TransportModule* create(TransportModuleConfigPtr& config);
}

}

#endif // TRANSPORTMODULEFACTORY_H
