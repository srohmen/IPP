#ifndef DIFFUSIONVELOCITYPOROSITYFUNC_H
#define DIFFUSIONVELOCITYPOROSITYFUNC_H

#include <palabos/dataProcessors/dataAnalysisWrapper2D.h>
#include <palabos/dataProcessors/dataAnalysisWrapper3D.h>

#include "initequilibriumconcentration.h"

namespace IPP
{

class DiffusionVelocityPorosityFunc
{
public:
    DiffusionVelocityPorosityFunc() = default;

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
        plb::applyProcessingFunctional<T,Descriptor,T>(
                    new SetDens, dotList, diffLattice, concentrations);
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
        plb::applyProcessingFunctional<T,Descriptor,T>(
                    new SetDens(concentration), dotList, diffLattice);
    }
private:

};

}

#endif // DIFFUSIONVELOCITYPOROSITYFUNC_H
