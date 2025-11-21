#ifndef SETDOTTONSOURCE_H
#define SETDOTTONSOURCE_H

#include "plbtypededuction.h"

#include <unordered_map>

#include "celltotalsdiff.h"
#include "domainiterator.h"

namespace IPP
{

template<typename T, size_t dim>
class SetDotToNSource : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD_N<T, dim>::value
{

public:
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;
    using CellToDiffMap = std::unordered_map<Dot, const CellTotalsDiff*>;

    SetDotToNSource(const CellToDiffMap& cellToDiff)
        : m_cellToDiff(cellToDiff)

    {

    }

    virtual SetDotToNSource<T, dim>* clone() const override
    {
        return new SetDotToNSource<T, dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        assert(modified.size() == 1);
        modified[0] = plb::modif::staticVariables;
    }


private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_N<T, dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using NTensorField = typename Traits::template argument<2>::type;


public:
    virtual void process(DotList dotList, NTensorField& field) override
    {
        const Dot fieldLocation = field.getLocation();

        for(const Dot& pos : dotList.dots)
        {
            const Dot absPos = pos + fieldLocation;

            const CellTotalsDiff& cellAndDiff = *m_cellToDiff.at(absPos);

            const std::vector<double>& sourceVec = cellAndDiff.totalsDiff;
            T* nTensor = DataAccess::get(field, pos);

            for(size_t iComp = 0; iComp < sourceVec.size(); ++iComp)
            {
                nTensor[iComp] = sourceVec[iComp];
            }
        }


    }

private:
    const CellToDiffMap& m_cellToDiff;
};


}

#endif // SETDOTTONSOURCE_H
