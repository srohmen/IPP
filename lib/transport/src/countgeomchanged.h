#ifndef COUNTGEOMCHANGED_H
#define COUNTGEOMCHANGED_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"
#include "domainiterator.h"
#include "permeabiltyflag.h"
#include "mpimanager.h"
#include "geometryupdate.h"

namespace IPP
{

template<typename T, typename MaskType, size_t dim>
class CountGeomChangedT :
        public PlbTypeDeduction::GetReductiveBoxProcessingFunctionalXD_SS<T, MaskType, dim>::value
{
public:
    CountGeomChangedT(const T& threshLower, const T& threshUpper)
        : countId(this->getStatistics().subscribeIntSum())
        , m_threshLower(threshLower)
        , m_threshUpper(threshUpper)
    {
        assert(m_threshLower <= m_threshUpper);
    }

    virtual CountGeomChangedT<T, MaskType, dim>* clone() const override
    {
        return new CountGeomChangedT<T, MaskType, dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::nothing;
    }

private:
    using BaseClass = typename PlbTypeDeduction::
    GetReductiveBoxProcessingFunctionalXD_SS<T, MaskType, dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;

    using Box = typename Traits::template argument<1>::type;
    using ScalarField = typename Traits::template argument<2>::type;
    using MaskField = typename Traits::template argument<3>::type;

public:
    virtual void process(Box domain, ScalarField& porosity, MaskField& permMask) override
    {
        plb::BlockStatistics& statistics = this->getStatistics();

        const auto offsetMask = plb::computeRelativeDisplacement(porosity, permMask);

        using Iterator = DataAccess::DomainIterator<dim>;
        const Iterator end = DataAccess::end(domain);
        for(Iterator it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto pos = *it;
            const T& newPorosity = DataAccess::get(porosity, pos);
            const MaskType& permFlag = DataAccess::get(permMask, pos + offsetMask);

            const GeometryUpdate::GoemChange change =
                    GeometryUpdate::getChange(m_threshLower, m_threshUpper,
                                              newPorosity, permFlag);

            if(change != GeometryUpdate::GS_None)
            {
                statistics.gatherIntSum(countId, 1);
            }

        }
    }

    plb::plint getCount() const
    {
        const plb::BlockStatistics& statistics = this->getStatistics();
        const plb::plint count = statistics.getIntSum(countId);
        return count;
    }

private:
    plb::plint countId;
    const T m_threshLower;
    const T m_threshUpper;
};


}

#endif // COUNTGEOMCHANGED_H
