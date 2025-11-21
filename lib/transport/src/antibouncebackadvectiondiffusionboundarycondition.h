#ifndef ANTIBOUNCEBACKADVECTIONDIFFUSIONBOUNDARYCONDITION_H
#define ANTIBOUNCEBACKADVECTIONDIFFUSIONBOUNDARYCONDITION_H

#include <palabos/atomicBlock/dataProcessingFunctional2D.h>
#include <palabos/atomicBlock/dataProcessingFunctional3D.h>
#include <palabos/boundaryCondition/boundaryCondition.h>
#include <palabos/complexDynamics/advectionDiffusionBoundaries.h>
#include <palabos/complexDynamics/adiabaticBoundaryProcessor3D.h>
#include <palabos/boundaryCondition/boundaryDynamics.h>
#include <palabos/complexDynamics/advectionDiffusionBoundaryCondition2D.h>
#include <palabos/complexDynamics/advectionDiffusionBoundaryCondition3D.h>
#include <palabos/complexDynamics/advectionDiffusionBoundaryInstantiator2D.h>
#include <palabos/complexDynamics/advectionDiffusionBoundaryInstantiator3D.h>

namespace IPP
{

template<typename T,
         template<typename U> class Descriptor,
         int direction, int orientation>
class AntiBounceBackAdvectionDiffusionBoundaryDynamics : public plb::StoreDensityDynamics<T,Descriptor>
{
public:
    /// Constructor
    AntiBounceBackAdvectionDiffusionBoundaryDynamics(plb::Dynamics<T,Descriptor>* baseDynamics,
                                                     bool automaticPrepareCollision = true)
        : plb::StoreDensityDynamics<T,Descriptor>(baseDynamics, automaticPrepareCollision)
    {

    }

    AntiBounceBackAdvectionDiffusionBoundaryDynamics(plb::HierarchicUnserializer& unserializer)
        : plb::StoreDensityDynamics<T,Descriptor>(0, false)
    {
        this->unserialize(unserializer);
    }

    /// Return a unique ID for this class.
    virtual int getId() const
    {
        return id;
    }

    /// Clone the object, based on its dynamic type
    virtual AntiBounceBackAdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>* clone() const
    {
        return new AntiBounceBackAdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
    }

    /// Execute completion scheme before base collision
    virtual void completePopulations(plb::Cell<T,Descriptor>& cell) const
    {
        typedef Descriptor<T> D;

        const T rho    = this->computeDensity(cell);
        const T rhoBar = D::rhoBar(rho);

        // TODO: add the closure impementation
        assert(false);
    }

private:
    static int id;
};


template<typename T, template<typename U> class Descriptor>
class AntiBounceBackAdvectionDiffusionBoundaryManager2D {
public:
    template<int direction, int orientation>
    static plb::BoundaryCompositeDynamics<T,Descriptor>*
    getAdvectionDiffusionBoundaryDynamics(plb::Dynamics<T,Descriptor>* baseDynamics)
    {
        return new AntiBounceBackAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>(baseDynamics);
    }

    template<int direction, int orientation>
    static plb::BoxProcessingFunctional2D_L<T,Descriptor>* getAdvectionDiffusionBoundaryProcessor()
    {
        return nullptr;
    }

    template<int xNormal, int yNormal>
    static plb::BoundaryCompositeDynamics<T,Descriptor>*
    getAdvectionDiffusionCornerDynamics(plb::Dynamics<T,Descriptor>* baseDynamics)
    {
        return new AntiBounceBackAdvectionDiffusionBoundaryDynamics<T, Descriptor, xNormal, yNormal>(baseDynamics);
    }

    template<int xNormal, int yNormal>
    static plb::BoxProcessingFunctional2D_L<T,Descriptor>* getAdvectionDiffusionCornerProcessor()
    {
        return nullptr;
    }
};


template<typename T, template<typename U> class Descriptor>
static plb::OnLatticeAdvectionDiffusionBoundaryCondition2D<T,Descriptor>*
        createAntiBounceBackAdvectionDiffusionBoundaryCondition2D()
{
    return new plb::AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor,
                       AntiBounceBackAdvectionDiffusionBoundaryManager2D<T,Descriptor> > ();
}






template<typename T, template<typename U> class Descriptor>
class AntiBounceBackAdvectionDiffusionBoundaryManager3D {
public:
    template<int direction, int orientation>
    static plb::BoundaryCompositeDynamics<T,Descriptor>*
    getAdvectionDiffusionBoundaryDynamics(plb::Dynamics<T,Descriptor>* baseDynamics)
    {
        return new AntiBounceBackAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>(baseDynamics);
    }

    template<int direction, int orientation>
    static plb::BoxProcessingFunctional3D_L<T,Descriptor>* getAdvectionDiffusionBoundaryProcessor()
    {
        return nullptr;
    }

    template<int xNormal, int yNormal>
    static plb::BoundaryCompositeDynamics<T,Descriptor>*
    getAdvectionDiffusionCornerDynamics(plb::Dynamics<T,Descriptor>* baseDynamics)
    {
        return new AntiBounceBackAdvectionDiffusionBoundaryDynamics<T, Descriptor, xNormal, yNormal>(baseDynamics);
    }

    template<int xNormal, int yNormal>
    static plb::BoxProcessingFunctional3D_L<T,Descriptor>* getAdvectionDiffusionCornerProcessor()
    {
        return nullptr;
    }
};

template<typename T, template<typename U> class Descriptor>
static plb::OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>*
        createAntiBounceBackAdvectionDiffusionBoundaryCondition3D()
{
    return new plb::AdvectionDiffusionBoundaryConditionInstantiator3D<T, Descriptor,
                       AntiBounceBackAdvectionDiffusionBoundaryManager3D<T,Descriptor> > ();
}


}

#endif // ANTIBOUNCEBACKADVECTIONDIFFUSIONBOUNDARYCONDITION_H
