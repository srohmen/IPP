// must be included first because of descriptor dependency...
#include <palabos/latticeBoltzmann/advectionDiffusionLattices.h>

#include "advectiondiffusionporositywithsourcetrtdynamics.hh"

#include <palabos/core/dynamics.hh>
#include <palabos/basicDynamics/isoThermalDynamics.hh>
#include <palabos/complexDynamics/trtDynamics.hh>

#include "externtemplatedescriptor.h"


EXPLICIT_TEMPLATE_SCALAR_DIFFPOROSITYDESCRIPTOR_3D(class IPP::AdvectionDiffusionPorosityWithSourceTRTDynamics)


