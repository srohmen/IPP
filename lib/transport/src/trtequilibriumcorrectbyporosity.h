#ifndef TRTEQUILIBRIUMCORRECTBYPOROSITY_H
#define TRTEQUILIBRIUMCORRECTBYPOROSITY_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"
#include "domainiterator.h"

namespace IPP
{

template<typename T, template<class U> class Descriptor>
class BoxTRTEquilibriumCorrectByPorosity
        : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_L<T,Descriptor>::value
{
public:
    BoxTRTEquilibriumCorrectByPorosity() = default;

    virtual BoxTRTEquilibriumCorrectByPorosity<T,Descriptor>* clone() const override
    {
        return new BoxTRTEquilibriumCorrectByPorosity<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_L<T,Descriptor>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, Lattice& lattice)  override
    {

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto dstPos = *it;
            plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, dstPos);

            const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
            const T& poros = *pPoros;


            // density calc inside of TRT dynamics is dividing by porosity already
            // const T1 rho = cell.computeDensity();
            T rho = cell[0];

            T sum = 0.0;
            for(size_t iPop = 1; iPop < Descriptor<T>::q; ++iPop)
            {
                const T f = cell[iPop];
                rho += f;
                sum += f;
            }

            cell[0] = rho * poros - sum;
        }
    }
};


template<typename T, template<class U> class Descriptor>
class DotTRTEquilibriumCorrectByPorosity
        : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD_L<T,Descriptor>::value
{
public:
    DotTRTEquilibriumCorrectByPorosity() = default;

    virtual DotTRTEquilibriumCorrectByPorosity<T,Descriptor>* clone() const override
    {
        return new DotTRTEquilibriumCorrectByPorosity<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_L<T,Descriptor>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;

public:
    virtual void process(const DotList& dotList, Lattice& lattice) override
    {
        for(const auto& dot : dotList.dots)
        {
            plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, dot);

            const T* pPoros = cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
            const T& poros = *pPoros;


            // density calc inside of TRT dynamics is dividing by porosity already
            // const T1 rho = cell.computeDensity();
            T rho = cell[0];

            T sum = 0.0;
            for(size_t iPop = 1; iPop < Descriptor<T>::q; ++iPop)
            {
                const T f = cell[iPop];
                rho += f;
                sum += f;
            }

            cell[0] = rho * poros - sum;
        }
    }
};


}

#endif // TRTEQUILIBRIUMCORRECTBYPOROSITY_H
