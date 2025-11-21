#ifndef INITEQUILIBRIUMCONCENTRATION_H
#define INITEQUILIBRIUMCONCENTRATION_H

#include "plbtypededuction.h"

#include "palabosdataaccess.h"
#include "domainiterator.h"

namespace IPP
{

template<typename T1, template<class U> class Descriptor, typename T2 = T1>
class BoxInitEquilibriumConcentration_LS
        : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_LS<T1,Descriptor,T2>::value
{
public:
    BoxInitEquilibriumConcentration_LS()
    {

    }

    virtual BoxInitEquilibriumConcentration_LS<T1,Descriptor,T2>* clone() const
    {
        return new BoxInitEquilibriumConcentration_LS<T1,Descriptor,T2>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
    }

private:    
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_LS<T1,Descriptor,T2>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;
    using ScalarField = typename Traits::template argument<3>::type;

public:
    virtual void process(Box domain, Lattice& lattice, ScalarField& rhoField ) override
    {
        auto offset = plb::computeRelativeDisplacement(lattice, rhoField);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto dstPos = *it;
            plb::Cell<T1,Descriptor>& cell = DataAccess::get(lattice, dstPos);

            const auto srcPos = dstPos + offset;
            const T2 rho = DataAccess::get(rhoField, srcPos);

            const T1 rhoBar = Descriptor<T1>::rhoBar(rho);
            plb::Array< T1, Descriptor<T1>::d > j;
            j.resetToZero();
            const T1 jSqr = T1();
            for (plb::plint iPop = 0; iPop < Descriptor<T1>::q; ++iPop)
            {
                cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
            }

        }
    }

};



template<typename T1, template<class U> class Descriptor, typename T2 = T1>
class DotInitEquilibriumConcentration_LS
        : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD_LS<T1,Descriptor,T2>::value
{
public:
    DotInitEquilibriumConcentration_LS()
    {

    }

    virtual DotInitEquilibriumConcentration_LS<T1,Descriptor,T2>* clone() const override
    {
        return new DotInitEquilibriumConcentration_LS<T1,Descriptor,T2>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_LS<T1,Descriptor,T2>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;
    using ScalarField = typename Traits::template argument<3>::type;

public:
    virtual void process(DotList& dotList,
                         Lattice& lattice,
                         ScalarField& rhoField ) override
    {
        for(plb::plint iDot = 0; iDot < dotList.getN(); ++iDot)
        {
            const auto& dot = dotList.getDot(iDot);

            // TODO: needs dot any offset correction?

            plb::Cell<T1,Descriptor>& cell = DataAccess::get(lattice, dot);
            const T2 rho = DataAccess::get(rhoField, dot);
            const T1 rhoBar = Descriptor<T1>::rhoBar(rho);

            plb::Array< T1, Descriptor<T1>::d > j;
            j.resetToZero();
            const T1 jSqr = T1();

            for (plb::plint iPop = 0; iPop < Descriptor<T1>::q; ++iPop)
            {
                cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
            }

        }
    }

};


template<typename T, template<class U> class Descriptor>
class DotInitEquilibriumConcentration_L : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD_L<T,Descriptor>::value
{
public:
    DotInitEquilibriumConcentration_L(const T& rho)
        : m_rho(rho)
    {

    }

    virtual DotInitEquilibriumConcentration_L<T,Descriptor>* clone() const override
    {
        return new DotInitEquilibriumConcentration_L<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
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
        for(plb::plint iDot = 0; iDot < dotList.getN(); ++iDot)
        {
            const auto& dot = dotList.getDot(iDot);

            // TODO: needs dot any offset correction?
            const T rhoBar = Descriptor<T>::rhoBar(m_rho);

            plb::Array< T, Descriptor<T>::d > j;
            j.resetToZero();
            const T jSqr = T();

            plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, dot);
            for (plb::plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
            {
                cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
            }

        }
    }

private:
    const T m_rho;
};



}

#endif // INITEQUILIBRIUMCONCENTRATION_H
