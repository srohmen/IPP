#ifndef EXTERNTEMPLATEHELPER_H
#define EXTERNTEMPLATEHELPER_H

#ifdef EXPLICIT_INSTANTS

#include "transporttypes.h"
#include "diffusionvelocitytransporttraits2d.h"
#include "diffusionvelocitytransporttraits3d.h"
#include "porositytrttransporttraits2d.h"
#include "porositytrttransporttraits3d.h"


#define EXTERN_TEMPLATE_TRANSPORT_TRAITS(CLASS_TYPE) \
    extern template CLASS_TYPE<IPP::PorosityTRTTransportTraits2D>; \
    extern template CLASS_TYPE<IPP::PorosityTRTTransportTraits3D>;
//    extern template CLASS_TYPE<IPP::DiffusionVelocityTransportTraits2D>;
//    extern template CLASS_TYPE<IPP::DiffusionVelocityTransportTraits3D>;


#define EXPLICIT_TEMPLATE_TRANSPORT_TRAITS(CLASS_TYPE) \
    template CLASS_TYPE<IPP::PorosityTRTTransportTraits2D>; \
    template CLASS_TYPE<IPP::PorosityTRTTransportTraits3D>;
//    template CLASS_TYPE<IPP::DiffusionVelocityTransportTraits2D>;
//    template CLASS_TYPE<IPP::DiffusionVelocityTransportTraits3D>;




#define EXTERN_TEMPLATE_SCALAR_DIM(CLASS_TYPE) \
    extern template CLASS_TYPE<IPP::TransportFloatingPointType, 2>; \
    extern template CLASS_TYPE<IPP::TransportFloatingPointType, 3>;

#define EXPLICIT_TEMPLATE_SCALAR_DIM_2D(CLASS_TYPE) \
    template CLASS_TYPE<IPP::TransportFloatingPointType, 2>;

#define EXPLICIT_TEMPLATE_SCALAR_DIM_3D(CLASS_TYPE) \
    template CLASS_TYPE<IPP::TransportFloatingPointType, 3>;

#define EXPLICIT_TEMPLATE_SCALAR_DIM(CLASS_TYPE) \
   EXPLICIT_TEMPLATE_SCALAR_DIM_2D(CLASS_TYPE) \
   EXPLICIT_TEMPLATE_SCALAR_DIM_3D(CLASS_TYPE)



#else

#define EXTERN_TEMPLATE_TRANSPORT_TRAITS(CLASS_TYPE)
#define EXPLICIT_TEMPLATE_TRANSPORT_TRAITS(CLASS_TYPE)

#define EXTERN_TEMPLATE_SCALAR_DIM(CLASS_TYPE)
#define EXPLICIT_TEMPLATE_SCALAR_DIM(CLASS_TYPE)
#define EXPLICIT_TEMPLATE_SCALAR_DIM_2D(CLASS_TYPE)
#define EXPLICIT_TEMPLATE_SCALAR_DIM_3D(CLASS_TYPE)


#endif

#endif // EXTERNTEMPLATEHELPER_H
