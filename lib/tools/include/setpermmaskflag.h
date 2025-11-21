#ifndef SETPERMMASKFLAG_H
#define SETPERMMASKFLAG_H

#include "plbtypededuction.h"

#include "palabosdataaccess.h"
#include "domainiterator.h"

namespace IPP
{

struct EnableFlag
{
    template<typename T>
    static void exec(const T& flag, T& permVal)
    {
        permVal |= flag;
    }
};

struct DisableFlag
{
    template<typename T>
    static void exec(const T& flag, T& permVal)
    {
        permVal &= ~flag;
    }
};



template<typename T, size_t dim, typename Operation>
class BoxSetPermMaskFlag : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_S<T,dim>::value
{
public:
    BoxSetPermMaskFlag(const T& flag)
        : m_flag(flag)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
    }

    virtual BoxSetPermMaskFlag<T, dim, Operation>* clone() const override
    {
        return new BoxSetPermMaskFlag<T, dim, Operation>(*this);
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using ScalarField = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, ScalarField& field) override
    {
        using Iterator = DataAccess::DomainIterator<dim>;
        const Iterator end = DataAccess::end(domain);
        for(Iterator it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& pos = *it;
            T& permVal = DataAccess::get(field, pos);
            Operation::exec(m_flag, permVal);
        }
    }


private:
    const T m_flag;
};



template<typename T, size_t dim, typename Operation>
class DotSetPermMaskFlag : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD_S<T,dim>::value
{
public:
    DotSetPermMaskFlag(const T& flag)
        : m_flag(flag)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
    }

    virtual DotSetPermMaskFlag<T, dim, Operation>* clone() const override
    {
        return new DotSetPermMaskFlag<T, dim, Operation>(*this);
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using ScalarField = typename Traits::template argument<2>::type;

public:
    virtual void process(const DotList& dotList, ScalarField& field) override
    {
        for(const auto& pos : dotList.dots)
        {
            T& permVal = DataAccess::get(field, pos);
            Operation::exec(m_flag, permVal);
        }
    }


private:
    const T m_flag;
};



}

#endif // SETPERMMASKFLAG_H
