#ifndef TRANSPORTMODULEFACTORY2D_H
#define TRANSPORTMODULEFACTORY2D_H

#include "transportmoduleconfigfwd.h"

namespace IPP
{

class TransportModule;

namespace TransportModuleFactory2D
{
    TransportModule* create(TransportModuleConfigPtr& config);
}


}

#endif // TRANSPORTMODULEFACTORY2D_H
