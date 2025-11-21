#ifndef PALABOSOBJECTFACTORY_H
#define PALABOSOBJECTFACTORY_H

#include "palabosobjectfactoryutils.h"
#include "plbtypededuction.h"
#include "ippbox.h"
#include "palabosconversiontools.h"
#include "setupblockperiodicity.h"

namespace IPP
{

template<size_t dim>
struct CreateSparseBlock;

template<>
struct CreateSparseBlock<2>
{
    static plb::SparseBlockStructure2D create(const plb::plint nx,
                                              const plb::plint ny,
                                              const plb::plint /*nz*/)
    {
        return plb::SparseBlockStructure2D(nx, ny);
    }
};


template<>
struct CreateSparseBlock<3>
{
    static plb::SparseBlockStructure3D create(const plb::plint nx,
                                              const plb::plint ny,
                                              const plb::plint nz)
    {
        return plb::SparseBlockStructure3D(nx, ny, nz);
    }
};


template<size_t dim>
class PalabosObjectFactory
{
public:
    PalabosObjectFactory(const plb::plint nx,
                         const plb::plint ny,
                         const plb::plint nz,
                         const std::vector<IPPBox3DInt>& decomp,
                         const std::vector<int>& periodicity)
        : nx(nx)
        , ny(ny)
        , nz(nz)
        , m_blockStructure(CreateSparseBlock<dim>::create(nx, ny, nz))
        , m_isBlocksDefined(false)
        , m_periodicity(periodicity)
    {
        for(size_t i = 0; i < decomp.size(); ++i)
        {
            const IPPBox3DInt& box = decomp[i];
            const auto plbBox = PalabosConversionTools::ToPlb<dim>::convertBox(box);
            m_blockStructure.addBlock(plbBox, m_blockStructure.nextIncrementalId());
        }

        if(decomp.empty() == false)
        {
            m_isBlocksDefined = true;
        }
    }



    template<typename T,
             template<typename U> class Descriptor,
             typename LatticeType,
             typename DynamicType>
    std::auto_ptr<LatticeType> createMultiLattice(DynamicType* dyn) const
    {
        if(m_isBlocksDefined)
        {
            auto field = PlbFactory<T>::template createMultiLattice<Descriptor, LatticeType>(m_blockStructure, dyn);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
        else
        {
            auto field = PlbFactory<T>::template createMultiLattice<LatticeType>(nx, ny, nz, dyn);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
    }

    template<typename T>
    auto createMultiScalarField(const T& iniVal = T()) const
    {
        if(m_isBlocksDefined)
        {
            auto field = PlbFactory<T>::createMultiScalarField(m_blockStructure, iniVal);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, field);
            return field;
        }
        else
        {
            auto field = PlbFactory<T>::createMultiScalarField(nx, ny, nz, iniVal);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, field);
            return field;
        }
    }

    template<typename T>
    auto createMultiScalarFieldPtr(const T& iniVal = T()) const
    {
        if(m_isBlocksDefined)
        {
            auto field = PlbFactory<T>::createMultiScalarFieldPtr(m_blockStructure, iniVal);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
        else
        {
            auto field = PlbFactory<T>::createMultiScalarFieldPtr(nx, ny, nz, iniVal);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
    }

    template<typename T, size_t tensorDim>
    auto createMultiTensorField() const
    {
        if(m_isBlocksDefined)
        {
            auto field = PlbFactory<T>::template createMultiTensorField<tensorDim>(m_blockStructure);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, field);
            return field;
        }
        else
        {
            auto field = PlbFactory<T>::template createMultiTensorField<tensorDim>(nx, ny, nz);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, field);
            return field;
        }
    }


    template<typename T, size_t tensorDim>
    auto createMultiTensorFieldPtr() const
    {
        if(m_isBlocksDefined)
        {
            auto field = PlbFactory<T>::template createMultiTensorFieldPtr<tensorDim>(m_blockStructure);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
        else
        {
            auto field = PlbFactory<T>::template createMultiTensorFieldPtr<tensorDim>(nx, ny, nz);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
    }


    template<typename T>
    auto createMultiNTensorField(const size_t tensorDim) const
    {
        if(m_isBlocksDefined)
        {
            auto field = PlbFactory<T>::createMultiNTensorField(m_blockStructure, tensorDim);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, field);
            return field;
        }
        else
        {
            auto field = PlbFactory<T>::createMultiNTensorField(nx, ny, nz, tensorDim);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, field);
            return field;
        }
    }

    template<typename T>
    auto createMultiNTensorFieldPtr(const size_t tensorDim) const
    {
        if(m_isBlocksDefined)
        {
            auto field = PlbFactory<T>::createMultiNTensorFieldPtr(m_blockStructure, tensorDim);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
        else
        {
            auto field = PlbFactory<T>::createMultiNTensorFieldPtr(nx, ny, nz, tensorDim);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
    }

    template<typename T>
    auto createMultiNTensorFieldPtr(const size_t tensorDim, const T* ini) const
    {
        if(m_isBlocksDefined)
        {
            auto field = PlbFactory<T>::createMultiNTensorFieldPtr(m_blockStructure, tensorDim, ini);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
        else
        {
            auto field = PlbFactory<T>::createMultiNTensorFieldPtr(nx, ny, nz, tensorDim, ini);
            SetupBlockPeriodicity::definePeriodicity(m_periodicity, *field);
            return field;
        }
    }

private:
    template<typename T>
    using PlbFactory = PalabosObjectFactoryUtils<T, dim>;

    using SparseBlockStructure = typename PlbTypeDeduction::GetSparseBlockStructure_XD<dim>::value;


    plb::plint nx, ny, nz;
    SparseBlockStructure m_blockStructure;
    bool m_isBlocksDefined;

    const std::vector<int>& m_periodicity;


};

}


#endif // PALABOSOBJECTFACTORY_H
