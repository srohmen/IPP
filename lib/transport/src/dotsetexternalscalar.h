#ifndef DOTSETEXTERNALSCALAR_H
#define DOTSETEXTERNALSCALAR_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class DotSetExternalScalar :
        public PlbTypeDeduction::GetDotProcessingFunctionalXD_L<T,Descriptor>::value
{
public:
    DotSetExternalScalar(const int whichScalar, const T& externalScalar)
        : whichScalar(whichScalar)
        , externalScalar(externalScalar)
    {

    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_L<T,Descriptor>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;

public:

    virtual void process(DotList& dotList, Lattice& lattice) override
    {
        for(const auto& dot : dotList.dots)
        {
            plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, dot);
            *cell.getExternal(whichScalar) = externalScalar;
        }
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulk;
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
    }

    virtual DotSetExternalScalar<T,Descriptor>* clone() const override
    {
        return new DotSetExternalScalar<T,Descriptor>(*this);
    }

private:
    int whichScalar;
    T externalScalar;
};

}

#endif // DOTSETEXTERNALSCALAR_H
