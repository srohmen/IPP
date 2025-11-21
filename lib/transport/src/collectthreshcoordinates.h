#ifndef COLLECTTHRESHCOORDINATES_H
#define COLLECTTHRESHCOORDINATES_H

#include "plbtypededuction.h"


namespace IPP
{

template<typename T, size_t dim>
class CollectThreshCoordinates : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_S<T,dim>::value
{
private:
    using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;

public:
    CollectThreshCoordinates(const T& thresh, DotList& above, DotList& below);

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override;
    virtual plb::BlockDomain::DomainT appliesTo() const override;

    virtual CollectThreshCoordinates<T,dim>* clone() const override;


private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, InputField& field) override;


private:
    const T m_thresh;
    DotList& m_above;
    DotList& m_below;

};


}

#ifdef EXPLICIT_INSTANTS
#include "externtemplatehelper.h"
EXTERN_TEMPLATE_SCALAR_DIM(class IPP::CollectThreshCoordinates)
#else
#include "collectthreshcoordinates.hh"
#endif

#endif // COLLECTTHRESHCOORDINATES_H
