#include "transportmodulefactory3d.h"

#include "reactivepalabostransportmodule.h"
#include "inertpalabostransportmodule.h"

// #include "diffusionvelocitytransporttraits3d.h"
#include "porositytrttransporttraits3d.h"

namespace IPP
{

namespace TransportModuleFactory3D
{

TransportModule* create(TransportModuleConfigPtr& config)
{    
    if(config->latticeBoltzmannCollType == LBCT_DVSRT)
    {
        if(config->isChemistryEnabled)
        {
            // return new ReactivePalabosTransportModule<DiffusionVelocityTransportTraits3D>(config);
            return nullptr;
        }
        else
        {
            // return new InertPalabosTransportModule<DiffusionVelocityTransportTraits3D>(config);
            return nullptr;
        }
    }
    else if(config->latticeBoltzmannCollType == LBCT_PTRT)
    {
        if(config->isChemistryEnabled)
        {
            return new ReactivePalabosTransportModule<PorosityTRTTransportTraits3D>(config);
        }
        else
        {
            return new InertPalabosTransportModule<PorosityTRTTransportTraits3D>(config);
        }
    }
    else
    {
        throw std::runtime_error("unknown lattice boltzmann type: "
                                 + std::to_string(config->latticeBoltzmannCollType) );
    }
}

}

}
