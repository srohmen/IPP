#ifndef COLLECTPERMMASKFLAGCOORDINATES_H
#define COLLECTPERMMASKFLAGCOORDINATES_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"
#include "domainiterator.h"

namespace IPP
{

template<typename T, size_t dim>
class CollectPermMaskFlagCoordinates : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_S<T, dim>::value
{

private:
    using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;

public:
    CollectPermMaskFlagCoordinates(const T& flag, DotList& result)
        : m_flag(flag)
        , m_result(result)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulkAndEnvelope;
    }

    virtual CollectPermMaskFlagCoordinates<T,dim>* clone() const override
    {
        return new CollectPermMaskFlagCoordinates<T,dim>(*this);
    }


private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, InputField& field) override
    {
        const auto offset = field.getLocation();

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& pos = *it;
            const T& value = DataAccess::get(field, pos);

            if((value & m_flag) == m_flag)
            {
                const auto posGlobal = pos + offset;
                m_result.addDot(posGlobal);
            }
        }
    }


private:
    const T m_flag;
    DotList& m_result;

};



}


#endif // COLLECTPERMMASKFLAGCOORDINATES_H
