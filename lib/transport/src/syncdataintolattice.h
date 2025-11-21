#ifndef SYNCDATAINTOLATTICE_H
#define SYNCDATAINTOLATTICE_H

#include "plbtypededuction.h"

#include <palabos/core/cell.h>

#include "array_view.h"
#include "domainiterator.h"

namespace IPP
{

template<typename T, typename U, template<typename V> class Descriptor>
class SyncDataIntoLattice : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_LS<T, Descriptor>::value
{
    static constexpr size_t dim = Descriptor<T>::d;
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

public:
    SyncDataIntoLattice(const Dot& location, const av::array_view<const T, dim>& transportScalar,
                        const av::array_view<const U, dim>& conc)
        : m_location(location)
        , m_transportScalar(transportScalar)
        , m_conc(conc)
    {

    }

    virtual SyncDataIntoLattice<T, U, Descriptor>* clone() const override
    {
        return new SyncDataIntoLattice<T, U, Descriptor>(*this);
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulk;
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_LS<T, Descriptor>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;
    using ScalarField = typename Traits::template argument<3>::type;

public:
    virtual void process(Box domain, Lattice& lattice, ScalarField& porosity) override
    {
        const auto porosOffset = plb::computeRelativeDisplacement(lattice, porosity);
        const auto latticeLocation = lattice.getLocation();

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto pos = *it;
            plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, pos);

            const T& poros = DataAccess::get(porosity, pos + porosOffset);

            // correct again for absolute positions
            const auto arrPos = pos + latticeLocation - m_location;
            const T& transScalar = DataAccess::get(m_transportScalar, arrPos);
            const T& conc = DataAccess::get(m_conc, arrPos);

            using EF = typename Descriptor<T>::ExternalField;
            *cell.getExternal(EF::porosityBeginsAt) = poros;
            *cell.getExternal(EF::transportScalarBeginsAt) = transScalar;
            *cell.getExternal(EF::scalarBeginsAt) = conc;
        }
    }



private:
    const Dot m_location;
    const av::array_view<const T, dim>& m_transportScalar;
    const av::array_view<const U, dim>& m_conc;
};

}

#endif // SYNCDATAINTOLATTICE_H
