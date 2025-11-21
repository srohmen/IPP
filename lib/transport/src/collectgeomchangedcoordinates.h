#ifndef COLLECTGEOMCHANGEDCOORDINATES_H
#define COLLECTGEOMCHANGEDCOORDINATES_H

#include "plbtypededuction.h"

#include "domainiterator.h"
#include "palabosdataaccess.h"
#include "geometryupdate.h"

namespace IPP
{


template<typename T, typename MaskType, size_t dim>
class CollectGeomChangedCoordinates : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_SS<T, MaskType, dim>::value
{
private:
    using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;

public:
    CollectGeomChangedCoordinates(const T& threshLower, const T& threshUpper,
                                  DotList& perm, DotList& nonPerm)
        : m_threshLower(threshLower)
        , m_threshUpper(threshUpper)
        , m_perm(perm)
        , m_nonPerm(nonPerm)
    {
        assert(m_threshLower <= m_threshUpper);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulkAndEnvelope;
    }

    virtual CollectGeomChangedCoordinates<T, MaskType, dim>* clone() const override
    {
        return new CollectGeomChangedCoordinates<T, MaskType, dim>(*this);
    }


private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_SS<T,MaskType,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using PorosField = typename Traits::template argument<2>::type;
    using MaskField = typename Traits::template argument<3>::type;

public:
    virtual void process(Box domain, PorosField& porosField, MaskField& permMask) override
    {
        const auto offset = porosField.getLocation();
        const auto offsetMask = plb::computeRelativeDisplacement(porosField, permMask);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& pos = *it;
            const T& value = DataAccess::get(porosField, pos);

            const MaskType& permFlag = DataAccess::get(permMask, pos + offsetMask);

            const auto posGlobal = pos + offset;

            const GeometryUpdate::GoemChange change =
                    GeometryUpdate::getChange(m_threshLower, m_threshUpper, value, permFlag);

            if(change == GeometryUpdate::GS_None)
            {
                if(permFlag & PF_isPermeable)
                {
                    m_perm.addDot(posGlobal);
                }
                else
                {
                    m_nonPerm.addDot(posGlobal);
                }
            }
            else if(change == GeometryUpdate::GS_SolidToLiquid)
            {
                m_perm.addDot(posGlobal);
            }
            else if(change == GeometryUpdate::GS_LiquidToSolid)
            {
                m_nonPerm.addDot(posGlobal);
            }
        }
    }


private:
    const T m_threshLower;
    const T m_threshUpper;
    DotList& m_perm;
    DotList& m_nonPerm;

};


}


#endif // COLLECTGEOMCHANGEDCOORDINATES_H
