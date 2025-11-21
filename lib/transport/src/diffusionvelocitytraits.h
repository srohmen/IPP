#ifndef DIFFUSIONVELOCITYTRAITS_H
#define DIFFUSIONVELOCITYTRAITS_H

#include "partialbouncebackbgkdynamics.h"
#include "advectiondiffusionvelocitywithsourcebgkdynamics.h"

#include "diffusionvelocitydiffusionomegacalc.h"
#include "diffusionvelocityporosityfunc.h"

namespace IPP
{

struct DiffusionVelocityTraits
{
    using Scalar = double;

    template<typename T, template<typename U> class Descriptor>
    using HydrodynamicDynamics = PartialBounceBackBGKDynamics<T, Descriptor>;

    template<typename T, template<typename U> class Descriptor>
    using DiffusionDynamics = AdvectionDiffusionVelocityWithSourceBGKdynamics<T, Descriptor>;
    // using DiffusionDynamic = PartialBounceBackAdvectionDiffusionWithSourceRLBdynamics<Scalar, DiffusionDescT>;
    // using DiffusionDynamic = AdvectionDiffusionVelocityWithSourceRLBdynamics<Scalar, DiffusionDescT>;


    struct HelperFunctions
    {        
        static constexpr bool isTRT = false;

        template<typename T>
        using DiffusionOmegaCalc = DiffusionVelocityDiffusionOmegaCalc<T>;


        using PorosityFunc = DiffusionVelocityPorosityFunc;

    };
};

}

#endif // DIFFUSIONVELOCITYTRAITS_H
