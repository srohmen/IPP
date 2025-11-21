#ifndef COPYPERMEABLECELLSCONC_H
#define COPYPERMEABLECELLSCONC_H

#include "palabosdataaccess.h"

namespace IPP
{

template<typename T, size_t dim>
class CopyPermeableCellsConc : public PlbTypeDeduction::
        GetScalarFieldBoxProcessingFunctionalXD<T,dim>::value
{
public:
    CopyPermeableCellsConc()
    {

    }

    virtual CopyPermeableCellsConc<T, dim>* clone() const override
    {
        return new CopyPermeableCellsConc<T, dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::nothing;
        modified[2] = plb::modif::staticVariables;
    }


    using BaseClass =
    typename PlbTypeDeduction::GetScalarFieldDotProcessingFunctionalXD<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using FieldVec = typename Traits::template argument<2>::type;

    virtual void process(Box domain, FieldVec fields) override
    {
        const auto& inputConcField = *fields[0];
        const auto& permMask = *fields[1];
        auto& outputConcField = *fields[2];

        using Iterator = DataAccess::DomainIterator<dim>;
        const Iterator end = DataAccess::end(domain);
        for(Iterator it = DataAccess::begin(domain);
            it < end; ++it)
        {
            const auto pos = *it;

            const T& permVal = DataAccess::get(permMask, pos);
            if(permVal > 0.0)
            {
                const T& conc = DataAccess::get(inputConcField, pos);
                DataAccess::get(outputConcField, pos) = conc;
            }
        }
    }
};


}

#endif // COPYPERMEABLECELLSCONC_H
