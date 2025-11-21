#ifndef SURROUNDEDCOMPAREFUNCTIONAL_H
#define SURROUNDEDCOMPAREFUNCTIONAL_H

#include "plbtypededuction.h"
#include "domainiterator2d.h"
#include "domainiterator3d.h"
#include "palabosconversiontools.h"

namespace  IPP
{

template<typename BaseClass, typename Comparator>
class SurroundedCompareFunctional : public BaseClass
{
public:
    SurroundedCompareFunctional(const Comparator& comp)
        : m_comp(comp)
    {

    }


   // using BaseClass = typename PalabosTypeDeductionHelper::GetBoxProcessingFunctionalXD_SS<T,T,dim>::value;
    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;
    using OutputField = typename Traits::template argument<2>::type;


    virtual void process(Box domain, InputField& input, OutputField& output)
    {
        auto offset = plb::computeRelativeDisplacement(input, output);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto pos = *it;

            const auto outputCoord = pos + offset;
            auto& out = DataAccess::get(output, outputCoord);
            out = DataAccess::get(input, pos);

            const Box localBox = PalabosConversionTools::toBox(pos);
            const Box extendedDomain = localBox.enlarge(1);

            const auto endLocal = DataAccess::end(extendedDomain);
            for(auto itLocal = DataAccess::begin(extendedDomain);
                itLocal < endLocal; ++itLocal)
            {
                const auto posLocal = *itLocal;
                const auto& in = DataAccess::get(input, posLocal);
                // std::cout << in << " vs. " << out << std::endl;
                out = m_comp(in, out);
                // std::cout << "\t -> " << out << std::endl;
            }

        }
    }

    virtual SurroundedCompareFunctional<BaseClass, Comparator>* clone() const
    {
        return new SurroundedCompareFunctional<BaseClass, Comparator>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::staticVariables;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return plb::BlockDomain::bulk;
    }

private:
    const Comparator& m_comp;
};





template<typename T, size_t dim, typename Comparator>
class ScalarFieldSurroundedCompareFunctional : public SurroundedCompareFunctional<
        typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_SS<T,T,dim>::value,
        Comparator
        >
{
    using BaseClass =
        SurroundedCompareFunctional<
            typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_SS<T,T,dim>::value,
            Comparator
        >;

public:
    ScalarFieldSurroundedCompareFunctional(const Comparator& comp)
        : BaseClass(comp)
    {

    }


};

template<typename T, template<typename U> class Descriptor, typename Comparator>
class LatticeSurroundedCompareFunctional : public SurroundedCompareFunctional<
        typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_LS<T,Descriptor>::value,
        Comparator
        >
{
    using BaseClass =
        SurroundedCompareFunctional<
            typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_LS<T,Descriptor>::value,
            Comparator
        >;

public:
    LatticeSurroundedCompareFunctional(const Comparator& comp)
        : BaseClass(comp)
    {

    }

};


template<typename T>
struct MinCompare
{
    inline const T& operator()(const T& a, const T& b) const
    {
        return std::min(a, b);
    }
};


template<typename T, size_t dim>
class ScalarFieldSurroundedMinFunctional : public ScalarFieldSurroundedCompareFunctional<T, dim, MinCompare<T>>
{
public:
    ScalarFieldSurroundedMinFunctional()
        : ScalarFieldSurroundedCompareFunctional<T, dim, MinCompare<T>>(MinCompare<T>())
    {

    }
};

template<typename T, template<typename U> class Descriptor>
class LatticeSurroundedMinFunctional : public LatticeSurroundedCompareFunctional<T, Descriptor, MinCompare<T>>
{
public:
    LatticeSurroundedMinFunctional()
        : LatticeSurroundedCompareFunctional<T, Descriptor, MinCompare<T>>(MinCompare<T>())
    {

    }
};



template<typename T>
struct MaxCompare
{
    inline const T& operator()(const T& a, const T& b) const
    {
        return std::max(a, b);
    }
};

template<typename T, size_t dim>
class ScalarFieldSurroundedMaxFunctional : public ScalarFieldSurroundedCompareFunctional<T, dim, MaxCompare<T>>
{
public:
    ScalarFieldSurroundedMaxFunctional()
        : ScalarFieldSurroundedCompareFunctional<T, dim, MaxCompare<T>>(MaxCompare<T>())
    {

    }
};

template<typename T, template<typename U> class Descriptor>
class LatticeSurroundedMaxFunctional : public LatticeSurroundedCompareFunctional<T, Descriptor, MaxCompare<T>>
{
public:
    LatticeSurroundedMaxFunctional()
        : LatticeSurroundedCompareFunctional<T, Descriptor, MaxCompare<T>>(MaxCompare<T>())
    {

    }
};


}

#endif // SURROUNDEDCOMPAREFUNCTIONAL_H
