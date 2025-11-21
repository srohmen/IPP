#ifndef COLLECTSOLIDCOORDINATESWITHPERMNEIGHBOR_H
#define COLLECTSOLIDCOORDINATESWITHPERMNEIGHBOR_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"
#include "permeabiltyflag.h"
#include "domainiterator.h"


namespace IPP
{

template<typename T, size_t dim>
class CollectSolidCoordinatesWithPermNeighbor : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_S<T,dim>::value
{
private:
    using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;

public:
    CollectSolidCoordinatesWithPermNeighbor(DotList& coordList)
        : m_coordList(coordList)
    {

    }

    virtual CollectSolidCoordinatesWithPermNeighbor<T,dim>* clone() const override
    {
        return new CollectSolidCoordinatesWithPermNeighbor<T,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulkAndEnvelope;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, InputField& permMaskField) override
    {
        const auto origin = permMaskField.getLocation();

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& pos = *it;
            T& maskVal = DataAccess::get(permMaskField, pos);

            if(     ((maskVal & PF_isPermeable) == false)
               &&   ((maskVal & PF_isInertSolid) == false) )
            {
                const unsigned int nPermNeighbors = (unsigned int)maskVal >> 8;
                assert(nPermNeighbors <= dim * 2);

                if(nPermNeighbors > 0)
                {
                    const auto result = pos + origin;
                    m_coordList.addDot(result);

                    maskVal |= PF_isInterface;
                }
                else
                {
                    // solid cells with no permeable neighbors are never interface
                    maskVal &= ~PF_isInterface;
                }
            }
            else
            {
                // permeable and inert cells are never interface
                maskVal &= ~PF_isInterface;
            }
        }
    }


private:
    DotList& m_coordList;
};

}

#endif // COLLECTSOLIDCOORDINATESWITHPERMNEIGHBOR_H
