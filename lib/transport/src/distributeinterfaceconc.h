#ifndef DISTRIBUTEINTERFACECONC_H
#define DISTRIBUTEINTERFACECONC_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"
#include "dotfindneighbors.h"

namespace IPP
{

template<typename T, size_t dim>
class DistributeInterfaceConc : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD_SS<T,T,dim>::value
{
public:
    DistributeInterfaceConc(const T& weight)
        : m_weight(weight)
    {

    }

    virtual DistributeInterfaceConc<T, dim>* clone() const override
    {
        return new DistributeInterfaceConc<T, dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::staticVariables;
    }

    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_SS<T,T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using MaskField = typename Traits::template argument<2>::type;
    using ConcField = typename Traits::template argument<3>::type;

    virtual void process(DotList dotList, MaskField& maskField, ConcField& concField) override
    {
        for(const auto& interfacePos : dotList.dots)
        {
            // assert(PalabosDataAccess::get(maskField, interfacePos) < 1.0);

            T& interfaceConc = DataAccess::get(concField, interfacePos);

            // find neighbors
            // 2D -> n = 4
            // 3D -> n = 6
            const DotList neighList = DotFindNeighbors::find(interfacePos);
            assert(neighList.getN() == dim * 2);

            for(const auto& neigh : neighList.dots)
            {
                const T& permVal = DataAccess::get(maskField, neigh);

                if(permVal > 0.0)
                {
                    T& cCurr = DataAccess::get(concField, neigh);

                    const T cDiff = interfaceConc - cCurr;
                    const T cDiffSource = m_weight * cDiff;

                    cCurr += cDiffSource;
                }
            }

            interfaceConc = 0.0;
        }
    }

private:
    const T m_weight;

};

}

#endif // DISTRIBUTEINTERFACECONC_H

