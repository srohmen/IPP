#ifndef TRANSPORTMODULEFACTORY3D_H
#define TRANSPORTMODULEFACTORY3D_H

#include "transportmoduleconfigfwd.h"

namespace IPP
{
class TransportModule;

namespace TransportModuleFactory3D
{
    TransportModule* create(TransportModuleConfigPtr& config);
}

}

#endif // TRANSPORTMODULEFACTORY3D_H
