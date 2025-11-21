#ifndef PUSHBOUNCEBACKPOPULATIONS_H
#define PUSHBOUNCEBACKPOPULATIONS_H

#include "plbtypededuction.h"

#include <palabos/latticeBoltzmann/indexTemplates.h>

#include "palabosdataaccess.h"
#include "domainiterator.h"
#include "permeabiltyflag.h"

namespace IPP
{

template<size_t dim>
struct ToDot;

template<>
struct ToDot<2>
{
    static plb::Dot2D convert(const int arr[2])
    {
        return plb::Dot2D(arr[0], arr[1]);
    }
};

template<>
struct ToDot<3>
{
    static plb::Dot3D convert(const int arr[3])
    {
        return plb::Dot3D(arr[0], arr[1], arr[2]);
    }
};


template<typename T, template<typename U> class Descriptor, typename MaskType>
class PushBounceBackPopulations : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_LS<T,Descriptor, MaskType>::value
{
public:
    PushBounceBackPopulations()
    {

    }

    virtual PushBounceBackPopulations<T, Descriptor, MaskType>* clone() const override
    {
        return new PushBounceBackPopulations<T, Descriptor, MaskType>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulk;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_LS<T,Descriptor, MaskType>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;
    using MaskField = typename Traits::template argument<3>::type;

public:
    virtual void process(Box domain, Lattice& lattice, MaskField& maskField) override
    {
        const auto offsetMask = plb::computeRelativeDisplacement(lattice, maskField);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto cellPos = *it;

            const MaskType& maskVal = DataAccess::get(maskField, cellPos + offsetMask);

            if(maskVal & PF_isPermeable)
            {
                // pull into perm neighbors only

                plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, cellPos);

#if !defined(NDEBUG)
                plb::BounceBack<T,Descriptor>* toCheck =
                        dynamic_cast<plb::BounceBack<T,Descriptor>*>(&cell.getDynamics());
                assert(toCheck == nullptr);
#endif

                for (plb::plint iPop = 1; iPop < Descriptor<T>::q; ++iPop)
                {
                    T& ownVal = cell[iPop];

                    const plb::plint iOpp = plb::indexTemplates::opposite<Descriptor<T>>(iPop);

                    const auto& dir = Descriptor<T>::c[iOpp];
                    const auto dirVec = ToDot<Descriptor<T>::d>::convert(dir);
                    const auto nextCellPos = cellPos + dirVec;

                    plb::Cell<T,Descriptor>& nextCell = DataAccess::get(lattice, nextCellPos);
                    const MaskType& nextCellMaskVal = DataAccess::get(maskField, nextCellPos + offsetMask);

                    // shifting signed integer is undefined behavior
                    const unsigned int nPermNeighbors = (unsigned int)nextCellMaskVal >> 8;

                    // only pull from solid neighbors
                    if((nextCellMaskVal & PF_isPermeable) == false && nPermNeighbors > 0)
                    {
                        ownVal += nextCell[iOpp];
                        nextCell[iOpp] = 0.0;
                    }

                }
            }

        }

        setToZero(domain, lattice, maskField);
    }

private:

    static void setToZero(Box domain, Lattice& lattice, MaskField& maskField)
    {
        const auto offsetMask = plb::computeRelativeDisplacement(lattice, maskField);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto cellPos = *it;

            const MaskType& maskVal = DataAccess::get(maskField, cellPos + offsetMask);

            if((maskVal & PF_isPermeable) == false)
            {
                plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, cellPos);

#if !defined(NDEBUG)
                plb::BounceBack<T,Descriptor>* toCheck =
                        dynamic_cast<plb::BounceBack<T,Descriptor>*>(&cell.getDynamics());
                assert(toCheck);
#endif

                // residual population must be zero already
                assert(cell[0] == 0);

                cell.getRawPopulations().resetToZero();
            }
        }

    }
};


}

#endif // PUSHBOUNCEBACKPOPULATIONS_H
