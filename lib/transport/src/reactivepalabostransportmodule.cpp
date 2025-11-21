#ifdef EXPLICIT_INSTANTS

#include "reactivepalabostransportmodule.hh"

namespace IPP
{

//template class ReactivePalabosTransportModule<DiffusionVelocityTransportTraits2D>;
//template class ReactivePalabosTransportModule<DiffusionVelocityTransportTraits3D>;

template class ReactivePalabosTransportModule<PorosityTRTTransportTraits2D>;
template class ReactivePalabosTransportModule<PorosityTRTTransportTraits3D>;

}


#endif

