#include "transportmodulefactory2d.h"

#include "reactivepalabostransportmodule.h"
#include "inertpalabostransportmodule.h"

// #include "diffusionvelocitytransporttraits2d.h"
#include "porositytrttransporttraits2d.h"

namespace IPP
{

namespace TransportModuleFactory2D
{

TransportModule* create(TransportModuleConfigPtr& config)
{
    if(config->latticeBoltzmannCollType == LBCT_DVSRT)
    {
        if(config->isChemistryEnabled)
        {
            // return new ReactivePalabosTransportModule<DiffusionVelocityTransportTraits2D>(config);
            return nullptr;
        }
        else
        {
            // return new InertPalabosTransportModule<DiffusionVelocityTransportTraits2D>(config);
            return nullptr;
        }
    }
    else if(config->latticeBoltzmannCollType == LBCT_PTRT)
    {
        if(config->isChemistryEnabled)
        {
            return new ReactivePalabosTransportModule<PorosityTRTTransportTraits2D>(config);
        }
        else
        {
            return new InertPalabosTransportModule<PorosityTRTTransportTraits2D>(config);
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
