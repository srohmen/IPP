#ifndef SETDOTSTOSCALAR_H
#define SETDOTSTOSCALAR_H

#include "palabosdataaccess.h"
#include "plbtypededuction.h"

namespace IPP
{


template<typename T, size_t dim>
class SetDotsToScalar_S :
        public PlbTypeDeduction::GetDotProcessingFunctionalXD_S<T,dim>::value
{
public:
    SetDotsToScalar_S(const T& value)
        : m_value(value)
    {

    }

    virtual SetDotsToScalar_S<T,dim>* clone() const override
    {
        return new SetDotsToScalar_S<T,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulkAndEnvelope;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using Field = typename Traits::template argument<2>::type;

public:
    virtual void process(DotList const& dotList, Field& field) override
    {
        for(const auto& dot : dotList.dots)
        {
            DataAccess::get(field, dot) = m_value;
        }
    }

private:
    const T m_value;
};


template<typename T, size_t dim>
class SetDotsToScalar_SS :
        public PlbTypeDeduction::GetDotProcessingFunctionalXD_SS<T,T,dim>::value
{
public:
    SetDotsToScalar_SS()
    {

    }

    virtual SetDotsToScalar_SS<T,dim>* clone() const
    {
        return new SetDotsToScalar_SS<T,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::staticVariables;
    }    

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulkAndEnvelope;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_SS<T,T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using InputField1 = typename Traits::template argument<2>::type;
    using InputField2 = typename Traits::template argument<3>::type;
    static_assert(std::is_same<InputField1, InputField2>::value, "not the same types");

public:
    virtual void process(DotList const& dotList,
                         InputField1& field1,
                         InputField2& field2)
    {
        for(const auto& dot : dotList.dots)
        {
            const T value = DataAccess::get(field1, dot);
            DataAccess::get(field2, dot) = value;
        }
    }

private:
};

}

#endif // SETDOTSTOSCALAR_H
