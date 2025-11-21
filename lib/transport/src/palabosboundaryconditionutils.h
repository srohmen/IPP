#ifndef PALABOSBOUNDARYCONDITIONUTILS_H
#define PALABOSBOUNDARYCONDITIONUTILS_H

#include <palabos/complexDynamics/advectionDiffusionBoundaryCondition2D.h>
#include <palabos/multiBlock/multiBlockLattice2D.h>
#include <palabos/boundaryCondition/boundaryCondition.h>
#include <palabos/complexDynamics/advectionDiffusionBoundaryCondition3D.h>
#include <palabos/boundaryCondition/boundaryCondition2D.h>
#include <palabos/boundaryCondition/boundaryCondition3D.h>

namespace IPP
{

namespace PalabosBoundaryConditionUtils
{

static inline plb::Box2D enlargeForPeriodicity(const plb::Box2D& block, const std::vector<int>& periodicBC)
{
    plb::Array<plb::plint,4> arr = block.to_plbArray();

    const bool isXperiodic = std::find(periodicBC.begin(), periodicBC.end(), 0) != periodicBC.end();
    if(isXperiodic)
    {
        --arr[0];
        ++arr[1];
    }

    const bool isYperiodic = std::find(periodicBC.begin(), periodicBC.end(), 1) != periodicBC.end();
    if(isYperiodic)
    {
        --arr[2];
        ++arr[3];
    }

    plb::Box2D extendedBlock;
    extendedBlock.from_plbArray(arr);
    return extendedBlock;
}

template<typename T, template<typename U> class Descriptor>
void setTemperatureConditionOnBlockBoundaries(plb::OnLatticeAdvectionDiffusionBoundaryCondition2D<T,Descriptor>& bc,
                                              plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                                              const plb::Box2D& block,
                                              const plb::Box2D& applicationDomain,
                                              const std::vector<int>& periodicBC = std::vector<int>(),
                                              const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    // only dirichlet type is supported for now
    assert(bcType == plb::boundary::dirichlet);

    const plb::plint boundaryWidth = 1;
    const plb::BlockSurface2D surf(block, boundaryWidth);
    plb::Box2D intersection;

    // boundary walls
    if (plb::intersect(surf.edge0N(), applicationDomain, intersection))
    {
        plb::Box2D wholeEdge = surf.edge0N();
        plb::merge(wholeEdge, surf.cornerNN());
        plb::merge(wholeEdge, surf.cornerNP());

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary0N(intersection, lattice);
    }

    if (plb::intersect(surf.edge0P(), applicationDomain, intersection))
    {
        plb::Box2D wholeEdge = surf.edge0P();
        plb::merge(wholeEdge, surf.cornerPN());
        plb::merge(wholeEdge, surf.cornerPP());

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary0P(intersection, lattice);
    }

    if (plb::intersect(surf.edge1N(), applicationDomain, intersection))
    {
        plb::Box2D wholeEdge = surf.edge1N();
        plb::merge(wholeEdge, surf.cornerNN());
        plb::merge(wholeEdge, surf.cornerPN());

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary1N(intersection, lattice);
    }

    if (plb::intersect(surf.edge1P(), applicationDomain, intersection))
    {
        plb::Box2D wholeEdge = surf.edge1P();
        plb::merge(wholeEdge, surf.cornerNP());
        plb::merge(wholeEdge, surf.cornerPP());

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary1P(intersection, lattice);
    }


    // boundary corners:
    // do not add corners since there implementations do not work with PTRT
    // corners are treated with edge dynamics above

    //    if (intersect(surf.cornerNN(), applicationDomain, intersection))
    //    {
    //        bc.addTemperatureCornerNN(intersection.x0, intersection.y0, lattice);
    //    }

    //    if (intersect(surf.cornerNP(), applicationDomain, intersection))
    //    {
    //        bc.addTemperatureCornerNP(intersection.x0, intersection.y0, lattice);
    //    }

    //    if (intersect(surf.cornerPN(), applicationDomain, intersection))
    //    {
    //        bc.addTemperatureCornerPN(intersection.x0, intersection.y0, lattice);
    //    }

    //    if (intersect(surf.cornerPP(), applicationDomain, intersection))
    //    {
    //        bc.addTemperatureCornerPP(intersection.x0, intersection.y0, lattice);
    //    }
}

template<typename T, template<typename U> class Descriptor>
void setTemperatureConditionOnBlockBoundaries(plb::OnLatticeAdvectionDiffusionBoundaryCondition2D<T,Descriptor>& bc,
                                              plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                                              const plb::Box2D& applicationDomain,
                                              const std::vector<int>& periodicBC = std::vector<int>(),
                                              const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    setTemperatureConditionOnBlockBoundaries(bc, lattice, lattice.getBoundingBox(),
                                             applicationDomain, periodicBC, bcType);
}


template<typename T, template<typename U> class Descriptor>
void setTemperatureConditionOnBlockBoundaries(plb::OnLatticeAdvectionDiffusionBoundaryCondition2D<T,Descriptor>& bc,
                                              plb::MultiBlockLattice2D<T, Descriptor>& lattice,
                                              const std::vector<int>& periodicBC = std::vector<int>(),
                                              const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    setTemperatureConditionOnBlockBoundaries(bc, lattice, lattice.getBoundingBox(),
                                             periodicBC, bcType);
}




////////////////////////////
// 3D part

static inline plb::Box3D enlargeForPeriodicity(const plb::Box3D& block, const std::vector<int>& periodicBC)
{
    plb::Array<plb::plint,6> arr = block.to_plbArray();

    const bool isXperiodic = std::find(periodicBC.begin(), periodicBC.end(), 0) != periodicBC.end();
    if(isXperiodic)
    {
        --arr[0];
        ++arr[1];
    }

    const bool isYperiodic = std::find(periodicBC.begin(), periodicBC.end(), 1) != periodicBC.end();
    if(isYperiodic)
    {
        --arr[2];
        ++arr[3];
    }

    const bool isZperiodic = std::find(periodicBC.begin(), periodicBC.end(), 2) != periodicBC.end();
    if(isZperiodic)
    {
        --arr[4];
        ++arr[5];
    }

    plb::Box3D extendedBlock;
    extendedBlock.from_plbArray(arr);
    return extendedBlock;
}


template<typename T, template<typename U> class Descriptor>
void setTemperatureConditionOnBlockBoundaries(plb::OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>& bc,
                                              plb::MultiBlockLattice3D<T,Descriptor>& lattice,
                                              const plb::Box3D& block,
                                              const plb::Box3D& applicationDomain,
                                              const std::vector<int>& periodicBC = std::vector<int>(),
                                              const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    // only dirichlet type is supported for now
    assert(bcType == plb::boundary::dirichlet);

    const plb::Box3D extendedBlock = enlargeForPeriodicity(block, periodicBC);
    const plb::plint boundaryWidth = 1;
    const plb::BlockSurface3D surf(extendedBlock, boundaryWidth);
    plb::Box3D intersection;

    if (plb::intersect(surf.surface0N(), applicationDomain, intersection))
    {
        plb::Box3D wholeEdge = surf.surface0N();

        plb::merge(wholeEdge, surf.edge1NN());
        plb::merge(wholeEdge, surf.edge1PN());

        plb::Box3D sideA = surf.edge2NN();
        plb::merge(sideA, surf.cornerNNN());
        plb::merge(sideA, surf.cornerNNP());
        plb::merge(wholeEdge, sideA);

        plb::Box3D sideB = surf.edge2NP();
        plb::merge(sideB, surf.cornerNPN());
        plb::merge(sideB, surf.cornerNPP());
        plb::merge(wholeEdge, sideB);

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary0N(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.surface0P(), applicationDomain, intersection))
    {
        plb::Box3D wholeEdge = surf.surface0P();

        plb::merge(wholeEdge, surf.edge1PP());
        plb::merge(wholeEdge, surf.edge1NP());

        plb::Box3D sideA = surf.edge2PP();
        plb::merge(sideA, surf.cornerPPP());
        plb::merge(sideA, surf.cornerPPN());
        plb::merge(wholeEdge, sideA);

        plb::Box3D sideB = surf.edge2PN();
        plb::merge(sideB, surf.cornerPNP());
        plb::merge(sideB, surf.cornerPNN());
        plb::merge(wholeEdge, sideB);

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary0P(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.surface1N(), applicationDomain, intersection))
    {
        plb::Box3D wholeEdge = surf.surface1N();

        plb::merge(wholeEdge, surf.edge0NN());
        plb::merge(wholeEdge, surf.edge0NP());

        plb::Box3D sideA = surf.edge2NN();
        plb::merge(sideA, surf.cornerNNN());
        plb::merge(sideA, surf.cornerNNP());
        plb::merge(wholeEdge, sideA);

        plb::Box3D sideB = surf.edge2PN();
        plb::merge(sideB, surf.cornerPNN());
        plb::merge(sideB, surf.cornerPNP());
        plb::merge(wholeEdge, sideB);

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary1N(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.surface1P(), applicationDomain, intersection))
    {
        plb::Box3D wholeEdge = surf.surface1P();

        plb::merge(wholeEdge, surf.edge0PP());
        plb::merge(wholeEdge, surf.edge0PN());

        plb::Box3D sideA = surf.edge2PP();
        plb::merge(sideA, surf.cornerPPP());
        plb::merge(sideA, surf.cornerPPN());
        plb::merge(wholeEdge, sideA);

        plb::Box3D sideB = surf.edge2NP();
        plb::merge(sideB, surf.cornerNPN());
        plb::merge(sideB, surf.cornerNPP());
        plb::merge(wholeEdge, sideB);

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary1P(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.surface2N(), applicationDomain, intersection))
    {
        plb::Box3D wholeEdge = surf.surface2N();

        plb::merge(wholeEdge, surf.edge0NN());
        plb::merge(wholeEdge, surf.edge0PN());

        plb::Box3D sideA = surf.edge1NN();
        plb::merge(sideA, surf.cornerNNN());
        plb::merge(sideA, surf.cornerNPN());
        plb::merge(wholeEdge, sideA);

        plb::Box3D sideB = surf.edge1NP();
        plb::merge(sideB, surf.cornerPNN());
        plb::merge(sideB, surf.cornerPPN());
        plb::merge(wholeEdge, sideB);

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary2N(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.surface2P(), applicationDomain, intersection))
    {
        plb::Box3D wholeEdge = surf.surface2P();

        plb::merge(wholeEdge, surf.edge0PP());
        plb::merge(wholeEdge, surf.edge0NP());

        plb::Box3D sideA = surf.edge1PN();
        plb::merge(sideA, surf.cornerNNP());
        plb::merge(sideA, surf.cornerNPP());
        plb::merge(wholeEdge, sideA);

        plb::Box3D sideB = surf.edge1PP();
        plb::merge(sideB, surf.cornerPNP());
        plb::merge(sideB, surf.cornerPPP());
        plb::merge(wholeEdge, sideB);

        plb::intersect(wholeEdge, applicationDomain, intersection);
        bc.addTemperatureBoundary2P(intersection, lattice, bcType);
    }



    // see comment in 2D implementation

    //    if (plb::intersect(surf.edge0NN(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge0NN(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge0NP(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge0NP(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge0PN(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge0PN(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge0PP(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge0PP(intersection, lattice, bcType);
    //    }

    //    if (plb::intersect(surf.edge1NN(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge1NN(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge1NP(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge1NP(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge1PN(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge1PN(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge1PP(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge1PP(intersection, lattice, bcType);
    //    }

    //    if (plb::intersect(surf.edge2NN(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge2NN(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge2NP(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge2NP(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge2PN(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge2PN(intersection, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.edge2PP(), applicationDomain, intersection)) {
    //        bc.addTemperatureEdge2PP(intersection, lattice, bcType);
    //    }

    //    if (plb::intersect(surf.cornerNNN(), applicationDomain, intersection)) {
    //        bc.addTemperatureCornerNNN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.cornerNNP(), applicationDomain, intersection)) {
    //        bc.addTemperatureCornerNNP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.cornerNPN(), applicationDomain, intersection)) {
    //        bc.addTemperatureCornerNPN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.cornerNPP(), applicationDomain, intersection)) {
    //        bc.addTemperatureCornerNPP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.cornerPNN(), applicationDomain, intersection)) {
    //        bc.addTemperatureCornerPNN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.cornerPNP(), applicationDomain, intersection)) {
    //        bc.addTemperatureCornerPNP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.cornerPPN(), applicationDomain, intersection)) {
    //        bc.addTemperatureCornerPPN(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    //    }
    //    if (plb::intersect(surf.cornerPPP(), applicationDomain, intersection)) {
    //        bc.addTemperatureCornerPPP(intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    //    }

}


template<typename T, template<typename U> class Descriptor>
void setTemperatureConditionOnBlockBoundaries(plb::OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>& bc,
                                              plb::MultiBlockLattice3D<T,Descriptor>& lattice,
                                              const plb::Box3D& applicationDomain,
                                              const std::vector<int>& periodicBC = std::vector<int>(),
                                              const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    setTemperatureConditionOnBlockBoundaries(bc, lattice, lattice.getBoundingBox(),
                                             applicationDomain, periodicBC, bcType);
}


template<typename T, template<typename U> class Descriptor>
void setTemperatureConditionOnBlockBoundaries(plb::OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>& bc,
                                              plb::MultiBlockLattice3D<T, Descriptor>& lattice,
                                              const std::vector<int>& periodicBC = std::vector<int>(),
                                              const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    setTemperatureConditionOnBlockBoundaries(bc, lattice, lattice.getBoundingBox(),
                                             periodicBC, bcType);
}



// velocity

template<typename T, template<typename U> class Descriptor>
void setVelocityConditionOnBlockBoundaries(plb::OnLatticeBoundaryCondition2D<T,Descriptor>& bc,
                                           plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                                           const plb::Box2D& block,
                                           const plb::Box2D& applicationDomain,
                                           const std::vector<int>& periodicBC = std::vector<int>(),
                                           const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    const plb::Box2D extendedBlock = enlargeForPeriodicity(block, periodicBC);
    const plb::plint boundaryWidth = 1;
    const plb::BlockSurface2D surf(extendedBlock, boundaryWidth);
    plb::Box2D intersection;


    if (plb::intersect(surf.edge0N(), applicationDomain, intersection)) {
        bc.addVelocityBoundary0N(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge0P(), applicationDomain, intersection)) {
        bc.addVelocityBoundary0P(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge1N(), applicationDomain, intersection)) {
        bc.addVelocityBoundary1N(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge1P(), applicationDomain, intersection)) {
        bc.addVelocityBoundary1P(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.cornerNN(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerNN(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerNP(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerNP(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerPN(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerPN(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerPP(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerPP(intersection.x0, intersection.y0, lattice, bcType);
    }
}

template<typename T, template<typename U> class Descriptor>
void setVelocityConditionOnBlockBoundaries(plb::OnLatticeBoundaryCondition2D<T,Descriptor>& bc,
                                           plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                                           const plb::Box2D& applicationDomain,
                                           const std::vector<int>& periodicBC = std::vector<int>(),
                                           const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    setVelocityConditionOnBlockBoundaries(bc, lattice, lattice.getBoundingBox(),
                                          applicationDomain, periodicBC, bcType);
}

template<typename T, template<typename U> class Descriptor>
void setVelocityConditionOnBlockBoundaries(plb::OnLatticeBoundaryCondition2D<T,Descriptor>& bc,
                                           plb::MultiBlockLattice2D<T, Descriptor>& lattice,
                                           const std::vector<int>& periodicBC = std::vector<int>(),
                                           const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    setVelocityConditionOnBlockBoundaries(bc, lattice, lattice.getBoundingBox(),
                                          periodicBC, bcType);
}




template<typename T, template<typename U> class Descriptor>
void setVelocityConditionOnBlockBoundaries(plb::OnLatticeBoundaryCondition3D<T,Descriptor>& bc,
                                           plb::MultiBlockLattice3D<T,Descriptor>& lattice,
                                           const plb::Box3D& block,
                                           const plb::Box3D& applicationDomain,
                                           const std::vector<int>& periodicBC = std::vector<int>(),
                                           const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    const plb::Box3D extendedBlock = enlargeForPeriodicity(block, periodicBC);
    const plb::plint boundaryWidth = 1;
    const plb::BlockSurface3D surf(extendedBlock, boundaryWidth);
    plb::Box3D intersection;


    if (plb::intersect(surf.surface0N(), applicationDomain, intersection)) {
        bc.addVelocityBoundary0N(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.surface0P(), applicationDomain, intersection)) {
        bc.addVelocityBoundary0P(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.surface1N(), applicationDomain, intersection)) {
        bc.addVelocityBoundary1N(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.surface1P(), applicationDomain, intersection)) {
        bc.addVelocityBoundary1P(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.surface2N(), applicationDomain, intersection)) {
        bc.addVelocityBoundary2N(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.surface2P(), applicationDomain, intersection)) {
        bc.addVelocityBoundary2P(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.edge0NN(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge0NN(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge0NP(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge0NP(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge0PN(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge0PN(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge0PP(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge0PP(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.edge1NN(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge1NN(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge1NP(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge1NP(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge1PN(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge1PN(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge1PP(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge1PP(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.edge2NN(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge2NN(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge2NP(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge2NP(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge2PN(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge2PN(intersection, lattice, bcType);
    }
    if (plb::intersect(surf.edge2PP(), applicationDomain, intersection)) {
        bc.addExternalVelocityEdge2PP(intersection, lattice, bcType);
    }

    if (plb::intersect(surf.cornerNNN(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerNNN (
                    intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerNNP(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerNNP (
                    intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerNPN(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerNPN (
                    intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerNPP(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerNPP (
                    intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerPNN(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerPNN (
                    intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerPNP(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerPNP (
                    intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerPPN(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerPPN (
                    intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (plb::intersect(surf.cornerPPP(), applicationDomain, intersection)) {
        bc.addExternalVelocityCornerPPP (
                    intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
}

template<typename T, template<typename U> class Descriptor>
void setVelocityConditionOnBlockBoundaries(plb::OnLatticeBoundaryCondition3D<T,Descriptor>& bc,
                                           plb::MultiBlockLattice3D<T,Descriptor>& lattice,
                                           const plb::Box3D& applicationDomain,
                                           const std::vector<int>& periodicBC = std::vector<int>(),
                                           const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    setVelocityConditionOnBlockBoundaries(bc, lattice, lattice.getBoundingBox(),
                                          applicationDomain, periodicBC, bcType);
}

template<typename T, template<typename U> class Descriptor>
void setVelocityConditionOnBlockBoundaries(plb::OnLatticeBoundaryCondition3D<T,Descriptor>& bc,
                                           plb::MultiBlockLattice3D<T, Descriptor>& lattice,
                                           const std::vector<int>& periodicBC = std::vector<int>(),
                                           const plb::boundary::BcType& bcType = plb::boundary::dirichlet)
{
    setVelocityConditionOnBlockBoundaries(bc, lattice, lattice.getBoundingBox(),
                                          periodicBC, bcType);
}


} // end of namespace PalabosBoundaryConditionUtils

} // end of namespace IPP

#endif // PALABOSBOUNDARYCONDITIONUTILS_H
