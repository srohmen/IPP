#ifndef HYDRODYNAMICBOUNDARYCONDITION_H
#define HYDRODYNAMICBOUNDARYCONDITION_H

#include <memory>

#include "palabosobjectfactoryutils.h"

namespace IPP
{

namespace PalabosBoundaryConditions
{

template<size_t dim>
struct HydrodynamicBoundaryCondition
{
public:
    HydrodynamicBoundaryCondition() = delete;

    template<typename T,
             template<typename U> class Descriptor,
             template<typename U, template<typename V> class Desc> class Lattice,
             typename Box
             >
    static void setDensityBoundaryCondition(Lattice<T, Descriptor>& lattice,
                                            const Box& domain,
                                            const std::vector<int>& /*periodicBC*/,
                                            const plb::boundary::BcType& bcType)
    {
        auto* velBC = PalabosObjectFactoryUtils<T, dim>::template createLocalBoundaryCondition<Descriptor>();
        velBC->setPressureConditionOnBlockBoundaries(lattice, domain, bcType);
        delete velBC;
    }

    template<typename T,
             template<typename U> class Descriptor,
             template<typename U, template<typename V> class Desc> class Lattice,
             typename Box
             >
    static void setVelocityBoundaryCondition(Lattice<T, Descriptor>& lattice,
                                             const Box& domain,
                                             const std::vector<int>& /*periodicBC*/,
                                             const plb::boundary::BcType& bcType)
    {
        auto* velBC = PalabosObjectFactoryUtils<T, dim>::template createLocalBoundaryCondition<Descriptor>();
        velBC->setVelocityConditionOnBlockBoundaries(lattice, domain, bcType);
        delete velBC;
    }

};


}

}


#endif // HYDRODYNAMICBOUNDARYCONDITION_H
