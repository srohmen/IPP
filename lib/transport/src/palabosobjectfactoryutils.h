#ifndef PALABOSOBJECTFACTORYUTILS_H
#define PALABOSOBJECTFACTORYUTILS_H


#include <type_traits>
#include <memory>

#include <palabos2D.h>
#include <palabos2D.hh>

#include <palabos3D.h>
#include <palabos3D.hh>

#include "ippexception.h"
#include "plbtypededuction.h"
#include "is_equal.h"

namespace IPP
{

// we are relying on RVO and move c++ move semantics here.
// in worst case of bad compiler the copy is just a map of pointers

template<typename Scalar, size_t dim>
struct PalabosObjectFactoryUtils;

template<typename Scalar>
struct PalabosObjectFactoryUtils<Scalar, 2>
{
    static plb::Dot2D makeDot(plb::plint x, plb::plint y, plb::plint /*z*/)
    {
        return plb::Dot2D(x, y);
    }



    template<template<typename U> class Descriptor,
             typename LatticeType,
             typename DynamicType>
    static std::auto_ptr<LatticeType>
    createMultiLattice(plb::SparseBlockStructure2D const& blockStructure, DynamicType* dyn)
    {
        typedef std::auto_ptr<LatticeType> LatticeTypePtr;
        LatticeTypePtr lattice(new LatticeType(
                                   plb::MultiBlockManagement2D(
                                       blockStructure,
                                       plb::defaultMultiBlockPolicy2D().getThreadAttribution(), 1 ),
                                   plb::defaultMultiBlockPolicy2D().getBlockCommunicator(),
                                   plb::defaultMultiBlockPolicy2D().getCombinedStatistics(),
                                   plb::defaultMultiBlockPolicy2D().getMultiCellAccess<Scalar, Descriptor>(),
                                   dyn)
                               );

        return lattice;
    }


    // multi proc types
    template<typename LatticeType, typename DynamicType>
    static std::auto_ptr<LatticeType>
    createMultiLattice(const size_t nx, const size_t ny, const size_t nz, DynamicType* dyn)
    {
        IPPCheck::assertCheck(nz == 1);

        typedef std::auto_ptr<LatticeType> LatticeTypePtr;
        LatticeTypePtr lattice(new LatticeType(nx, ny, dyn));
        return lattice;
    }



    typedef plb::MultiScalarField2D<Scalar> MultiScalarField;

    static MultiScalarField createMultiScalarField(plb::SparseBlockStructure2D const& blockStructure,
                                                   const Scalar& iniVal = Scalar())
    {
        return MultiScalarField(plb::MultiBlockManagement2D(
                                    blockStructure,
                                    plb::defaultMultiBlockPolicy2D().getThreadAttribution(), 1 ),
                                plb::defaultMultiBlockPolicy2D().getBlockCommunicator(),
                                plb::defaultMultiBlockPolicy2D().getCombinedStatistics(),
                                plb::defaultMultiBlockPolicy2D().getMultiScalarAccess<Scalar>(),
                                iniVal);

    }

    static MultiScalarField createMultiScalarField(const size_t nx,
                                                   const size_t ny,
                                                   const size_t nz,
                                                   const Scalar& iniVal = Scalar())
    {
        IPPCheck::assertCheck(nz == 1);
        MultiScalarField field(nx, ny, iniVal);
        return field;
    }


    static MultiScalarField createMultiScalarField(const plb::Box2D& box,
                                                   const Scalar& iniVal = Scalar())
    {
        return createMultiScalarField(box.getNx(), box.getNy(), 1, iniVal);
    }

    template<typename Vector>
    static MultiScalarField createMultiScalarField(const Vector& vec,
                                                   const Scalar& iniVal = Scalar())
    {
        return createMultiScalarField(vec[0], vec[1], 1, iniVal);
    }

    static MultiScalarField createMultiScalarField(const plb::Dot2D& vec,
                                                   const Scalar& iniVal = Scalar())
    {
        return createMultiScalarField(vec.x, vec.y, 1, iniVal);
    }




    static MultiScalarField* createMultiScalarFieldPtr(plb::SparseBlockStructure2D const& blockStructure,
                                                       const Scalar& iniVal = Scalar())
    {
        return new MultiScalarField(plb::MultiBlockManagement2D(
                                        blockStructure,
                                        plb::defaultMultiBlockPolicy2D().getThreadAttribution(), 1 ),
                                    plb::defaultMultiBlockPolicy2D().getBlockCommunicator(),
                                    plb::defaultMultiBlockPolicy2D().getCombinedStatistics(),
                                    plb::defaultMultiBlockPolicy2D().getMultiScalarAccess<Scalar>(),
                                    iniVal);
    }


    static MultiScalarField* createMultiScalarFieldPtr(const size_t nx,
                                                       const size_t ny,
                                                       const size_t nz,
                                                       const Scalar& iniVal = Scalar())
    {
        IPPCheck::assertCheck(nz == 1);
        MultiScalarField* field = new MultiScalarField(nx, ny, iniVal);
        return field;
    }

    static MultiScalarField* createMultiScalarFieldPtr(const plb::Box2D& box,
                                                       const Scalar& iniVal = Scalar())
    {
        return createMultiScalarFieldPtr(box.getNx(), box.getNy(), 1, iniVal);
    }

    template<typename Vector>
    static MultiScalarField* createMultiScalarFieldPtr(const Vector& vec,
                                                       const Scalar& iniVal = Scalar())
    {
        return createMultiScalarFieldPtr(vec[0], vec[1], 1, iniVal);
    }

    static MultiScalarField* createMultiScalarFieldPtr(const plb::Dot2D& vec,
                                                       const Scalar& iniVal = Scalar())
    {
        return createMultiScalarFieldPtr(vec.x, vec.y, 1, iniVal);
    }


    template<size_t TensorDim>
    using MultiTensorField = plb::MultiTensorField2D<Scalar, TensorDim>;

    template<size_t tensorDim>
    static MultiTensorField<tensorDim> createMultiTensorField(plb::SparseBlockStructure2D const& blockStructure)
    {
        return MultiTensorField<tensorDim>(plb::MultiBlockManagement2D(
                                               blockStructure,
                                               plb::defaultMultiBlockPolicy2D().getThreadAttribution(), 1 ),
                                           plb::defaultMultiBlockPolicy2D().getBlockCommunicator(),
                                           plb::defaultMultiBlockPolicy2D().getCombinedStatistics(),
                                           plb::defaultMultiBlockPolicy2D().getMultiTensorAccess<Scalar, tensorDim>());

    }


    template<size_t tensorDim>
    static MultiTensorField<tensorDim> createMultiTensorField(const size_t nx,
                                                              const size_t ny,
                                                              const size_t nz)
    {
        IPPCheck::assertCheck(nz == 1);
        MultiTensorField<tensorDim> field(nx, ny);
        return field;
    }

    template<size_t tensorDim>
    static MultiTensorField<tensorDim> createMultiTensorField(const plb::Box2D& box)
    {
        MultiTensorField<tensorDim> field(box.getNx(), box.getNy());
        return field;
    }


    template<size_t tensorDim>
    static MultiTensorField<tensorDim>* createMultiTensorFieldPtr(plb::SparseBlockStructure2D const& blockStructure)
    {
        return new MultiTensorField<tensorDim>(plb::MultiBlockManagement2D(
                                                   blockStructure,
                                                   plb::defaultMultiBlockPolicy2D().getThreadAttribution(), 1 ),
                                               plb::defaultMultiBlockPolicy2D().getBlockCommunicator(),
                                               plb::defaultMultiBlockPolicy2D().getCombinedStatistics(),
                                               plb::defaultMultiBlockPolicy2D().getMultiTensorAccess<Scalar, tensorDim>());
    }

    template<size_t tensorDim>
    static MultiTensorField<tensorDim>* createMultiTensorFieldPtr(const size_t nx,
                                                                  const size_t ny,
                                                                  const size_t nz)
    {
        IPPCheck::assertCheck(nz == 1);
        MultiTensorField<tensorDim>* field = new MultiTensorField<tensorDim>(nx, ny);
        return field;
    }



    using MultiNTensorField = plb::MultiNTensorField2D<Scalar>;

    static MultiNTensorField createMultiNTensorField(plb::SparseBlockStructure2D const& blockStructure,
                                                     const size_t nDim)
    {
        return MultiNTensorField(nDim,
                                 plb::MultiBlockManagement2D(
                                     blockStructure,
                                     plb::defaultMultiBlockPolicy2D().getThreadAttribution(), 1 ),
                                 plb::defaultMultiBlockPolicy2D().getBlockCommunicator(),
                                 plb::defaultMultiBlockPolicy2D().getCombinedStatistics(),
                                 plb::defaultMultiBlockPolicy2D().getMultiNTensorAccess<Scalar>());
    }

    static MultiNTensorField createMultiNTensorField(const size_t nx,
                                                     const size_t ny,
                                                     const size_t nz,
                                                     const size_t nDim)
    {
        IPPCheck::assertCheck(nz == 1);
        return MultiNTensorField(nx, ny, nDim);
    }

    template<typename Vector>
    static MultiNTensorField createMultiNTensorField(const Vector& vec,
                                                     const size_t nDim)
    {
        return MultiNTensorField(vec[0], vec[1], nDim);
    }

    static MultiNTensorField createMultiNTensorField(const size_t nx,
                                                     const size_t ny,
                                                     const size_t nz,
                                                     const size_t nDim,
                                                     const Scalar* init)
    {
        IPPCheck::assertCheck(nz == 1);
        return MultiNTensorField(nx, ny, nDim, init);
    }

    template<typename Vector>
    static MultiNTensorField createMultiNTensorField(const Vector& vec,
                                                     const size_t nDim,
                                                     const Scalar* init)
    {
        return MultiNTensorField(vec[0], vec[1], nDim, init);
    }


    static MultiNTensorField* createMultiNTensorFieldPtr(plb::SparseBlockStructure2D const& blockStructure,
                                                         const size_t nDim)
    {
        return new MultiNTensorField(nDim, plb::MultiBlockManagement2D(
                                         blockStructure,
                                         plb::defaultMultiBlockPolicy2D().getThreadAttribution(), 1 ),
                                     plb::defaultMultiBlockPolicy2D().getBlockCommunicator(),
                                     plb::defaultMultiBlockPolicy2D().getCombinedStatistics(),
                                     plb::defaultMultiBlockPolicy2D().getMultiNTensorAccess<Scalar>());
    }

    static MultiNTensorField* createMultiNTensorFieldPtr(plb::SparseBlockStructure2D const& blockStructure,
                                                         const size_t nDim,
                                                         const Scalar* ini)
    {
        return new MultiNTensorField(nDim, ini,
                                     plb::MultiBlockManagement2D(
                                         blockStructure,
                                         plb::defaultMultiBlockPolicy2D().getThreadAttribution(), 1 ),
                                     plb::defaultMultiBlockPolicy2D().getBlockCommunicator(),
                                     plb::defaultMultiBlockPolicy2D().getCombinedStatistics(),
                                     plb::defaultMultiBlockPolicy2D().getMultiNTensorAccess<Scalar>()
                                     );
    }



    static MultiNTensorField* createMultiNTensorFieldPtr(const size_t nx,
                                                         const size_t ny,
                                                         const size_t nz,
                                                         const size_t nDim)
    {
        IPPCheck::assertCheck(nz == 1);
        return new MultiNTensorField(nx, ny, nDim);
    }

    template<typename Vector>
    static MultiNTensorField* createMultiNTensorFieldPtr(const Vector& vec,
                                                         const size_t nDim)
    {
        return new MultiNTensorField(vec[0], vec[1], nDim);
    }

    static MultiNTensorField* createMultiNTensorFieldPtr(const size_t nx,
                                                         const size_t ny,
                                                         const size_t nz,
                                                         const size_t nDim,
                                                         const Scalar* init)
    {
        IPPCheck::assertCheck(nz == 1);
        return new MultiNTensorField(nx, ny, nDim, init);
    }

    template<typename Vector>
    static MultiNTensorField* createMultiNTensorFieldPtr(const Vector& vec,
                                                         const size_t nDim,
                                                         const Scalar* init)
    {
        return new MultiNTensorField(vec[0], vec[1], nDim, init);
    }


    //////////////////////////

    // serial types
    typedef plb::ScalarField2D<Scalar> ScalarField;

    static ScalarField createScalarField(const plb::Box2D& box,
                                         const Scalar& iniVal = Scalar())
    {
        ScalarField field(box.getNx(), box.getNy(), 1, iniVal);
        return field;
    }

    static ScalarField createScalarField(const size_t nx, const size_t ny, const size_t nz,
                                         const Scalar& iniVal = Scalar())
    {
        IPPCheck::assertCheck(nz == 1);
        ScalarField field(nx, ny, iniVal);
        return field;
    }
    //////////////////////////

    // boundary condition types
    template<template<class U> class Descriptor>
    static plb::OnLatticeBoundaryCondition2D<Scalar,Descriptor>*
    createLocalBoundaryCondition()
    {
        return plb::createLocalBoundaryCondition2D<Scalar,Descriptor>();
    }

    template<template<class U> class Descriptor>
    static plb::OnLatticeAdvectionDiffusionBoundaryCondition2D<Scalar,Descriptor>*
    createLocalAdvectionDiffusionBoundaryCondition()
    {
        return plb::createLocalAdvectionDiffusionBoundaryCondition2D<Scalar,Descriptor>();
    }
    //////////////////////////

    static plb::Box2D createBox(const plb::plint x0, const plb::plint x1,
                                const plb::plint y0, const plb::plint y1,
                                const plb::plint /*z0*/, const plb::plint /*z1*/)
    {
        return plb::Box2D(x0, x1, y0, y1);
    }
};

template<typename Scalar>
struct PalabosObjectFactoryUtils<Scalar, 3>
{

    static plb::Dot3D makeDot(plb::plint x, plb::plint y, plb::plint z)
    {
        return plb::Dot3D(x, y, z);
    }

    template<template<typename U> class Descriptor,
             typename LatticeType,
             typename DynamicType>
    static std::auto_ptr<LatticeType>
    createMultiLattice(plb::SparseBlockStructure3D const& blockStructure, DynamicType* dyn)
    {
        typedef std::auto_ptr<LatticeType> LatticeTypePtr;
        LatticeTypePtr lattice(new LatticeType(
                                   plb::MultiBlockManagement3D(
                                       blockStructure,
                                       plb::defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
                                   plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                   plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                   plb::defaultMultiBlockPolicy3D().getMultiCellAccess<Scalar,Descriptor>(),
                                   dyn)
                               );

        return lattice;
    }

    // multi proc types
    template<typename LatticeType, typename DynamicType>
    static std::auto_ptr<LatticeType>
    createMultiLattice(const size_t nx, const size_t ny, const size_t nz, DynamicType* dyn)
    {
        typedef std::auto_ptr<LatticeType> LatticeTypePtr;
        LatticeTypePtr lattice(new LatticeType(nx, ny, nz, dyn));
        return lattice;
    }


    typedef plb::MultiScalarField3D<Scalar> MultiScalarField;

    static MultiScalarField createMultiScalarField(plb::SparseBlockStructure3D const& blockStructure,
                                                   const Scalar& iniVal = Scalar())
    {
        MultiScalarField field(plb::MultiBlockManagement3D(
                                   blockStructure,
                                   plb::defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
                               plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
                               plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
                               plb::defaultMultiBlockPolicy3D().getMultiScalarAccess<Scalar>(),
                               iniVal);
        return field;
    }

    static MultiScalarField createMultiScalarField(const size_t nx, const size_t ny, const size_t nz,
                                                   const Scalar& iniVal = Scalar())
    {
        IPPCheck::assertCheck(nz > 1);
        MultiScalarField field(nx, ny, nz, iniVal);
        return field;
    }

    static MultiScalarField createMultiScalarField(const plb::Box3D& box,
                                                   const Scalar& iniVal = Scalar())
    {
        MultiScalarField field(box.getNx(), box.getNy(), box.getNz(), iniVal);
        return field;
    }

    template<typename Vector>
    static MultiScalarField createMultiScalarField(const Vector& vec,
                                                   const Scalar& iniVal = Scalar())
    {
        MultiScalarField field(vec[0], vec[1], vec[2], iniVal);
        return field;
    }

    static MultiScalarField createMultiScalarField(const plb::Dot3D& vec,
                                                   const Scalar& iniVal = Scalar())
    {
        return createMultiScalarField(vec.x, vec.y, vec.z, iniVal);
    }


    static MultiScalarField* createMultiScalarFieldPtr(plb::SparseBlockStructure3D const& blockStructure,
                                                       const Scalar& iniVal = Scalar())
    {
        return new MultiScalarField(plb::MultiBlockManagement3D(
                                        blockStructure,
                                        plb::defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
                                    plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                    plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                    plb::defaultMultiBlockPolicy3D().getMultiScalarAccess<Scalar>(),
                                    iniVal);
    }


    static MultiScalarField* createMultiScalarFieldPtr(const size_t nx, const size_t ny, const size_t nz,
                                                       const Scalar& iniVal = Scalar())
    {
        IPPCheck::assertCheck(nz > 1);
        MultiScalarField* field = new MultiScalarField(nx, ny, nz, iniVal);
        return field;
    }

    static MultiScalarField* createMultiScalarFieldPtr(const plb::Box3D& box,
                                                       const Scalar& iniVal = Scalar())
    {
        return createMultiScalarFieldPtr(box.getNx(), box.getNy(), box.getNz(), iniVal);
    }

    template<typename Vector>
    static MultiScalarField* createMultiScalarFieldPtr(const Vector& vec,
                                                       const Scalar& iniVal = Scalar())
    {
        return createMultiScalarFieldPtr(vec[0], vec[1], vec[2], iniVal);
    }

    static MultiScalarField* createMultiScalarFieldPtr(const plb::Dot3D& vec,
                                                       const Scalar& iniVal = Scalar())
    {
        return createMultiScalarFieldPtr(vec.x, vec.y, vec.z, iniVal);
    }


    template<size_t TensorDim>
    using MultiTensorField = plb::MultiTensorField3D<Scalar, TensorDim>;


    template<size_t tensorDim>
    static MultiTensorField<tensorDim> createMultiTensorField(plb::SparseBlockStructure3D const& blockStructure)
    {
        return MultiTensorField<tensorDim>(plb::MultiBlockManagement3D(
                                               blockStructure,
                                               plb::defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
                                           plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                           plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                           plb::defaultMultiBlockPolicy3D().getMultiTensorAccess<Scalar, tensorDim>());
    }


    template<size_t tensorDim>
    static MultiTensorField<tensorDim> createMultiTensorField(const size_t nx, const size_t ny, const size_t nz)
    {
        IPPCheck::assertCheck(nz > 1);
        MultiTensorField<tensorDim> field(nx, ny, nz);
        return field;
    }

    template<size_t tensorDim>
    static MultiTensorField<tensorDim> createMultiTensorField(const plb::Box3D& box)
    {
        IPPCheck::assertCheck(box.getNz() > 1);
        MultiTensorField<tensorDim> field(box.getNx(), box.getNy(), box.getNz());
        return field;
    }


    template<size_t tensorDim>
    static MultiTensorField<tensorDim>* createMultiTensorFieldPtr(plb::SparseBlockStructure3D const& blockStructure)
    {
        return new MultiTensorField<tensorDim>(plb::MultiBlockManagement3D(
                                                   blockStructure,
                                                   plb::defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
                                               plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                               plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                               plb::defaultMultiBlockPolicy3D().getMultiTensorAccess<Scalar, tensorDim>());
    }

    template<size_t tensorDim>
    static MultiTensorField<tensorDim>* createMultiTensorFieldPtr(const size_t nx, const size_t ny, const size_t nz)
    {
        IPPCheck::assertCheck(nz > 1);
        MultiTensorField<tensorDim>* field = new MultiTensorField<tensorDim>(nx, ny, nz);
        return field;
    }



    using MultiNTensorField = plb::MultiNTensorField3D<Scalar>;

    static MultiNTensorField createMultiNTensorField(plb::SparseBlockStructure3D const& blockStructure,
                                                     const size_t nDim)
    {
        return MultiNTensorField(nDim, plb::MultiBlockManagement3D(
                                     blockStructure,
                                     plb::defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
                                 plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                 plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                 plb::defaultMultiBlockPolicy3D().getMultiNTensorAccess<Scalar>());
    }

    static MultiNTensorField createMultiNTensorField(const size_t nx,
                                                     const size_t ny,
                                                     const size_t nz,
                                                     const size_t nDim)
    {
        IPPCheck::assertCheck(nz > 1);
        return MultiNTensorField(nx, ny, nz, nDim);
    }

    template<typename Vector>
    static MultiNTensorField createMultiNTensorField(const Vector& vec,
                                                     const size_t nDim)
    {
        return MultiNTensorField(vec[0], vec[1], vec[2], nDim);
    }

    static MultiNTensorField createMultiNTensorField(const size_t nx,
                                                     const size_t ny,
                                                     const size_t nz,
                                                     const size_t nDim,
                                                     const Scalar* init)
    {
        IPPCheck::assertCheck(nz > 1);
        return MultiNTensorField(nx, ny, nz, nDim, init);
    }

    template<typename Vector>
    static MultiNTensorField createMultiNTensorField(const Vector& vec,
                                                     const size_t nDim,
                                                     const Scalar* init)
    {
        return MultiNTensorField(vec[0], vec[1], vec[2], nDim, init);
    }


    static MultiNTensorField* createMultiNTensorFieldPtr(plb::SparseBlockStructure3D const& blockStructure,
                                                         const size_t nDim)
    {
        return new MultiNTensorField(nDim,
                                     plb::MultiBlockManagement3D(
                                         blockStructure,
                                         plb::defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
                                     plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                     plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                     plb::defaultMultiBlockPolicy3D().getMultiNTensorAccess<Scalar>());
    }


    static MultiNTensorField* createMultiNTensorFieldPtr(plb::SparseBlockStructure3D const& blockStructure,
                                                         const size_t nDim,
                                                         const Scalar* ini)
    {
        return new MultiNTensorField(nDim, ini,
                                     plb::MultiBlockManagement3D(
                                         blockStructure,
                                         plb::defaultMultiBlockPolicy3D().getThreadAttribution(), 1 ),
                                     plb::defaultMultiBlockPolicy3D().getBlockCommunicator(),
                                     plb::defaultMultiBlockPolicy3D().getCombinedStatistics(),
                                     plb::defaultMultiBlockPolicy3D().getMultiNTensorAccess<Scalar>());
    }

    static MultiNTensorField* createMultiNTensorFieldPtr(const size_t nx,
                                                         const size_t ny,
                                                         const size_t nz,
                                                         const size_t nDim)
    {
        IPPCheck::assertCheck(nz > 1);
        return new MultiNTensorField(nx, ny, nz, nDim);
    }

    template<typename Vector>
    static MultiNTensorField* createMultiNTensorFieldPtr(const Vector& vec,
                                                         const size_t nDim)
    {
        return new MultiNTensorField(vec[0], vec[1], vec[2], nDim);
    }

    static MultiNTensorField* createMultiNTensorFieldPtr(const size_t nx,
                                                         const size_t ny,
                                                         const size_t nz,
                                                         const size_t nDim,
                                                         const Scalar* init)
    {
        IPPCheck::assertCheck(nz > 1);
        return new MultiNTensorField(nx, ny, nz, nDim, init);
    }

    template<typename Vector>
    static MultiNTensorField* createMultiNTensorFieldPtr(const Vector& vec,
                                                         const size_t nDim,
                                                         const Scalar* init)
    {
        return new MultiNTensorField(vec[0], vec[1], vec[2], nDim, init);
    }



    //////////////////////////

    // serial types
    typedef plb::ScalarField3D<Scalar> ScalarField;
    static ScalarField createScalarField(const size_t nx, const size_t ny, const size_t nz,
                                         const Scalar& iniVal = Scalar())
    {
        IPPCheck::assertCheck(nz > 1);
        ScalarField field(nx, ny, nz, iniVal);
        return field;
    }
    //////////////////////////

    // boundary condition types
    template<template<class U> class Descriptor>
    static plb::OnLatticeBoundaryCondition3D<Scalar,Descriptor>*
    createLocalBoundaryCondition()
    {
        return plb::createLocalBoundaryCondition3D<Scalar,Descriptor>();
    }

    template<template<class U> class Descriptor>
    static plb::OnLatticeAdvectionDiffusionBoundaryCondition3D<Scalar,Descriptor>*
    createLocalAdvectionDiffusionBoundaryCondition()
    {
        return plb::createLocalAdvectionDiffusionBoundaryCondition3D<Scalar,Descriptor>();
    }
    //////////////////////////

    static plb::Box3D createBox(const plb::plint x0, const plb::plint x1,
                                const plb::plint y0, const plb::plint y1,
                                const plb::plint z0, const plb::plint z1)
    {
        return plb::Box3D(x0, x1, y0, y1, z0, z1);
    }
};

} // end of namespace IPP
#endif // PALABOSOBJECTFACTORYUTILS_H
