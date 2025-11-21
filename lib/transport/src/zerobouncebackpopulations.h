#ifndef ZEROBOUNCEBACKPOPULATIONS_H
#define ZEROBOUNCEBACKPOPULATIONS_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"
#include "domainiterator.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class ZeroBounceBackPopulations : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_L<T,Descriptor>::value
{
public:
    ZeroBounceBackPopulations()
    {

    }

    virtual ZeroBounceBackPopulations<T, Descriptor>* clone() const
    {
        return new ZeroBounceBackPopulations<T, Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::staticVariables;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulkAndEnvelope;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_L<T,Descriptor>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, Lattice& lattice) override
    {
        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto pos = *it;

            plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, pos);

            plb::BounceBack<T,Descriptor>* toCheck =
                    dynamic_cast<plb::BounceBack<T,Descriptor>*>(&cell.getDynamics());
            if (toCheck) {
                cell.getRawPopulations().resetToZero();
                *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt) = 0.0;
            }
        }
    }
};


}

#endif // ZEROBOUNCEBACKPOPULATIONS_H
