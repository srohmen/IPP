#include "transportmodulefactory.h"

#include "transportmoduleconfig.h"
#include "transportmodulefactory2d.h"
#include "transportmodulefactory3d.h"
#include "transportmodule.h"

namespace IPP
{


TransportModule* TransportModuleFactory::create(TransportModuleConfigPtr& config)
{
    TransportModule* transModule;

    // TODO: make proper 3D check
    if(config->is3DSimulation)
    {
        transModule = TransportModuleFactory3D::create(config);
    }
    else
    {
        transModule = TransportModuleFactory2D::create(config);
    }

    return transModule;

}



}
