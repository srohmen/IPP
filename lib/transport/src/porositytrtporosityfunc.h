#ifndef POROSITYTRTPOROSITYFUNC_H
#define POROSITYTRTPOROSITYFUNC_H

#include <palabos/dataProcessors/dataAnalysisWrapper2D.hh>

#include "initequilibriumconcentration.h"
#include "trtequilibriumcorrectbyporosity.h"

#include "scopedfloatingpointexception.h"

namespace IPP
{

struct PorosityTRTPorosityFunc
{
    template<typename T,
             template<typename U> class Descriptor,
             template<typename U, template<typename V> class Desc> class Lattice,
             typename ScalarField>
    static void apply(Lattice<T,Descriptor>& diffLattice,
                      ScalarField& concentrations)
    {
        using SetDens = BoxInitEquilibriumConcentration_LS<T, Descriptor>;
        plb::applyProcessingFunctional<T,Descriptor,T>(
                    new SetDens, diffLattice.getBoundingBox(),
                    diffLattice, concentrations);

        // TRT rest velocity eq has to be corrected by porosity
        using CorrectByPoros = BoxTRTEquilibriumCorrectByPorosity<T, Descriptor>;
        plb::applyProcessingFunctional(new CorrectByPoros, diffLattice.getBoundingBox(), diffLattice);

    }

    template<typename T,
             template<typename U> class Descriptor,
             typename DotList,
             template<typename U, template<typename V> class Desc> class Lattice,
             typename ScalarField>
    static void apply(const DotList& dotList,
                      Lattice<T,Descriptor>& diffLattice,
                      ScalarField& concentrations)
    {
        using SetDens = DotInitEquilibriumConcentration_LS<T, Descriptor>;
        plb::applyProcessingFunctional(new SetDens, dotList, diffLattice, concentrations);

        // TRT rest velocity eq has to be corrected by porosity
        using CorrectByPoros = DotTRTEquilibriumCorrectByPorosity<T, Descriptor>;
        plb::applyProcessingFunctional(new CorrectByPoros, dotList, diffLattice);
    }

    template<typename T,
             template<typename U> class Descriptor,
             typename DotList,
             template<typename U, template<typename V> class Desc> class Lattice>
    static void apply(const DotList& dotList,
                      Lattice<T,Descriptor>& diffLattice,
                      const T& concentration)
    {
        using SetDens = DotInitEquilibriumConcentration_L<T, Descriptor>;
        plb::applyProcessingFunctional(new SetDens(concentration), dotList, diffLattice);

        // TRT rest velocity eq has to be corrected by porosity
        using CorrectByPoros = DotTRTEquilibriumCorrectByPorosity<T, Descriptor>;
        plb::applyProcessingFunctional(new CorrectByPoros, dotList, diffLattice);
    }

};

}

#endif // POROSITYTRTPOROSITYFUNC_H
