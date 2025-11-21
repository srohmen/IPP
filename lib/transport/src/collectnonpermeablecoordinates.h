#ifndef COLLECTNONPERMEABLECOORDINATES_H
#define COLLECTNONPERMEABLECOORDINATES_H

#include "plbtypededuction.h"

#include "palabosdataaccess2d.h"
#include "palabosdataaccess3d.h"

#include "externtemplatehelper.h"

namespace IPP
{

template<typename T, size_t dim>
class CollectNonPermeableCoordinates :
        public PlbTypeDeduction::GetDotProcessingFunctionalXD_S<T,dim>::value
{
private:
    using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;

public:
    CollectNonPermeableCoordinates(const T& thresh, DotList& partiallyPermCells, DotList& nonPermCells)
        : m_thresh(thresh)
        , partiallyPermCells(partiallyPermCells)
        , nonPermCells(nonPermCells)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }

    virtual CollectNonPermeableCoordinates<T,dim>* clone() const override
    {
        return new CollectNonPermeableCoordinates<T,dim>(*this);
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using ScalarField = typename Traits::template argument<2>::type;

public:
    virtual void process(const DotList& dotList, ScalarField& field) override
    {
        const auto offset = field.getLocation();
        for(const auto& dot : dotList.dots)
        {
            const auto result = dot + offset;

            const T& val = DataAccess::get(field, dot);

            if(val < m_thresh)
            {
                nonPermCells.addDot(result);
            }
            else
            {
                partiallyPermCells.addDot(result);
            }
        }

    }

private:
    const T m_thresh;
    DotList& partiallyPermCells;
    DotList& nonPermCells;
};

EXTERN_TEMPLATE_SCALAR_DIM(class CollectNonPermeableCoordinates)

}

#endif // COLLECTNONPERMEABLECOORDINATES_H
