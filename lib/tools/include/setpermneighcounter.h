#ifndef SETPERMNEIGHCOUNTER_H
#define SETPERMNEIGHCOUNTER_H

#include "plbtypededuction.h"

#include "palabosdataaccess.h"
#include "dotfindneighbors.h"
#include "permeabiltyflag.h"
#include "domainiterator.h"


namespace IPP
{

template<typename T, size_t dim>
class SetPermNeighCounter : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_S<T, dim>::value
{
public:
    SetPermNeighCounter()
    {

    }

    virtual SetPermNeighCounter<T,dim>* clone() const
    {
        return new SetPermNeighCounter<T,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::staticVariables;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, InputField& permMask)
    {
        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& pos = *it;

            T& permVal = DataAccess::get(permMask, pos);

            const auto neighbors = DotFindNeighbors::find(pos);

            T counter = 0;

            for(const auto& neigh : neighbors.dots)
            {
                if(neigh == pos)
                {
                    continue;
                }

                const T& permValNeigh = DataAccess::get(permMask, neigh);

                if(permValNeigh & PF_isPermeable)
                {
                    ++counter;
                }
            }

            const T oldMask = permVal & 0xFF;

            permVal = counter << 8;
            permVal |= oldMask;
        }
    }

};
}

#endif // SETPERMNEIGHCOUNTER_H
