#ifndef DISTRIBUTEBOUNCEBACKPOPULATIONS_H
#define DISTRIBUTEBOUNCEBACKPOPULATIONS_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"
#include "dotfindneighbors.h"
#include "permeabiltyflag.h"
#include "domainiterator.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor, typename MaskType>
class DistributeBounceBackPopulations : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_LS<T, Descriptor, MaskType>::value
{
public:
    DistributeBounceBackPopulations()
    {

    }

    virtual DistributeBounceBackPopulations<T, Descriptor, MaskType>* clone() const
    {
        return new DistributeBounceBackPopulations<T, Descriptor, MaskType>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        assert(modified.size() == 2);
        modified[0] = plb::modif::staticVariables; // lattice
        modified[1] = plb::modif::nothing; // permMask
    }

private:

    static constexpr size_t dim = Descriptor<T>::d;

    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_LS<T, Descriptor, MaskType>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;
    using ScalarFieldMask = typename Traits::template argument<3>::type;


public:
    virtual void process(Box domain, Lattice& lattice, ScalarFieldMask& maskField) override
    {
        // when node switches from liquid to solid it happens that there is
        // still a residual population because the population was relaxated before.
        // its important to redistribute that to the liquid neibors, too

        const auto offsetMask = plb::computeRelativeDisplacement(lattice, maskField);

        using Iterator = DataAccess::DomainIterator<dim>;
        const Iterator end = DataAccess::end(domain);
        for(Iterator it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& latticePos = *it;

            plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, latticePos);

            const MaskType& maskVal = DataAccess::get(maskField, latticePos + offsetMask);
            if(maskVal & PF_isPermeable)
            {

                T& conc = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
                const T& porosity = *cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);

                // find neighbors
                // 2D -> n = 4
                // 3D -> n = 6
                auto allNeighList = DotFindNeighbors::find(latticePos);
                assert(allNeighList.getN() == dim * 2);

                for(const auto& latticeNeighPos : allNeighList.dots)
                {
                    const MaskType& maskNeighVal = DataAccess::get(maskField, latticeNeighPos + offsetMask);

                    if((maskNeighVal & PF_isPermeable) == false)
                    {
                        // shifting signed integer is undefined behavior
                        const unsigned int nPermNeighbors = (unsigned int)maskNeighVal >> 8;
                        assert(nPermNeighbors <= dim * 2);

                        if(nPermNeighbors > 0)
                        {
                            plb::Cell<T,Descriptor>& neighCell = DataAccess::get(lattice, latticeNeighPos);

                            T& neighConc = *neighCell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
                            const T& neighPorosity = *neighCell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
                            const T moleToDistribute = neighConc * neighPorosity;

                            const T molePerCell = moleToDistribute / nPermNeighbors;

                            const T concToAdd = molePerCell / porosity;
                            conc += concToAdd;
                        }
                    }
                }
            }

        }


        // set solid to zero
        for(Iterator it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& latticePos = *it;

            const MaskType& maskVal = DataAccess::get(maskField, latticePos + offsetMask);
            if((maskVal & PF_isPermeable) == false)
            {
                auto& cell = DataAccess::get(lattice, latticePos);

                T& conc = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);
                conc = 0;
            }
        }
    }
};


}

#endif // DISTRIBUTEBOUNCEBACKPOPULATIONS_H
