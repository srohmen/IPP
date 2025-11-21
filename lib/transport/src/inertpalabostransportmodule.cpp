#ifdef EXPLICIT_INSTANTS

#include "inertpalabostransportmodule.hh"

#include "diffusionvelocitytransporttraits2d.h"
#include "diffusionvelocitytransporttraits3d.h"
#include "porositytrttransporttraits2d.h"
#include "porositytrttransporttraits3d.h"

namespace IPP
{

//template class InertPalabosTransportModule<DiffusionVelocityTransportTraits2D>;
//template class InertPalabosTransportModule<DiffusionVelocityTransportTraits3D>;

template class InertPalabosTransportModule<PorosityTRTTransportTraits2D>;
template class InertPalabosTransportModule<PorosityTRTTransportTraits3D>;

}

#endif
