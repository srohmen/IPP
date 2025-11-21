#include "advectiondiffusionwithsourcetrtdynamics.h"

#include <palabos/core/dynamicsIdentifiers.h>

#include <palabos/latticeBoltzmann/dynamicsTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionDynamicsTemplates.h>

#include <palabos/latticeBoltzmann/momentTemplates.h>
#include <palabos/latticeBoltzmann/advectionDiffusionMomentTemplates.h>



namespace IPP
{

template<typename T, template<typename U> class Descriptor>
const T AdvectionDiffusionWithSourceTRTDynamics<T,Descriptor>::s_magic = 1.0/4.0;


template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionWithSourceTRTDynamics<T,Descriptor>::id =
        plb::meta::registerGeneralDynamics<T,Descriptor,AdvectionDiffusionWithSourceTRTDynamics<T,Descriptor> >("ADS_TRT");


template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionWithSourceTRTDynamics<T, Descriptor>::AdvectionDiffusionWithSourceTRTDynamics(const T &omega)
    : plb::TRTdynamics<T, Descriptor>(omega)
{

}

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionWithSourceTRTDynamics<T, Descriptor> *AdvectionDiffusionWithSourceTRTDynamics<T, Descriptor>::clone() const
{
    return new AdvectionDiffusionWithSourceTRTDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int AdvectionDiffusionWithSourceTRTDynamics<T, Descriptor>::getId() const
{
    return id;
}



}
