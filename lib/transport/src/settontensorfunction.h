#ifndef SETTONTENSORFUNCTION_H
#define SETTONTENSORFUNCTION_H

#include "plbtypededuction.h"
#include "domainiterator.h"

namespace IPP
{

template<typename T, class Func, size_t dim>
class SetToNTensorFunction : public
        PlbTypeDeduction::GetBoxProcessingFunctionalXD_N<T, dim>::value
{
public:
    SetToNTensorFunction(const Func& f)
        : m_func(f)
    {

    }

    virtual SetToNTensorFunction<T, Func, dim>* clone() const
    {
        return new SetToNTensorFunction<T, Func, dim>(*this);
    }

    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return plb::BlockDomain::bulk;
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::staticVariables;
    }


private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_N<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using NTensorField = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, NTensorField& field)
    {
        using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;
        const Dot offset = field.getLocation();

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const Dot& pos = *it;
            T* valArr = DataAccess::get(field, pos);
            const Dot outPos = pos + offset;
            m_func(outPos, valArr);
        }
    }

private:
    const Func& m_func;
};


}

#endif // SETTONTENSORFUNCTION_H
