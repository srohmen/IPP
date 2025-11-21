#ifndef ADVECTIONDIFFUSIONBOUNDARYCONDITION_H
#define ADVECTIONDIFFUSIONBOUNDARYCONDITION_H

#include <memory>

#include "palabosobjectfactoryutils.h"
#include "palabosboundaryconditionutils.h"

namespace IPP
{

namespace PalabosBoundaryConditions
{
template<size_t dim>
class AdvectionDiffusionBoundaryCondition
{
public:
    AdvectionDiffusionBoundaryCondition() = delete;

    template<typename T,
             template<typename U> class Descriptor,
             template<typename U, template<typename V> class Desc> class Lattice,
             typename Box
             >
    static void setDensityBoundaryCondition(Lattice<T, Descriptor>& lattice,
                                            const Box& domain,
                                            const std::vector<int>& periodicBC,
                                            const plb::boundary::BcType& bcType)
    {
        auto* concBC = PalabosObjectFactoryUtils<T,dim>::template
                createLocalAdvectionDiffusionBoundaryCondition<Descriptor>();

        using namespace PalabosBoundaryConditionUtils;
        setTemperatureConditionOnBlockBoundaries(*concBC, lattice, domain, periodicBC, bcType);

        delete concBC;
    }

    template<typename T,
             template<typename U> class Descriptor,
             template<typename U, template<typename V> class Desc> class Lattice,
             typename Box
             >
    static void setVelocityBoundaryCondition(Lattice<T, Descriptor>& lattice,
                                             const Box& domain,
                                             const std::vector<int>& periodicBC,
                                             const plb::boundary::BcType& bcType)
    {
        auto* velBC = PalabosObjectFactoryUtils<T,dim>::template
                createLocalBoundaryCondition<Descriptor>();

        namespace PlbBC = PalabosBoundaryConditionUtils;
        PlbBC::setVelocityConditionOnBlockBoundaries(*velBC, lattice, domain, periodicBC, bcType);

        delete velBC;
    }


};


}
}

#endif // ADVECTIONDIFFUSIONBOUNDARYCONDITION_H
