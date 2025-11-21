#ifndef POROSITYTRTTRAITS_H
#define POROSITYTRTTRAITS_H

#include "partialbouncebackbgkdynamics.h"
#include "advectiondiffusionporositywithsourcetrtdynamics.h"
#include "porositytrtdiffusionomegacalc.h"
#include "porositytrtporosityfunc.h"

namespace IPP
{

struct PorosityTRTTraits
{
    template<typename T, template<typename U> class Descriptor>
    using HydrodynamicDynamics = PartialBounceBackBGKDynamics<T, Descriptor>;

    template<typename T, template<typename U> class Descriptor>
    using DiffusionDynamics = AdvectionDiffusionPorosityWithSourceTRTDynamics<T, Descriptor>;

    struct HelperFunctions
    {        
        static constexpr bool isTRT = true;

        template<typename T>
        using DiffusionOmegaCalc = PorosityTRTDiffusionOmegaCalc<T>;

        using PorosityFunc = PorosityTRTPorosityFunc;
    };
};

}

#endif // POROSITYTRTTRAITS_H
