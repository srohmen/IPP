#ifndef DOTFINDMIN_H
#define DOTFINDMIN_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"

namespace IPP
{

// TODO: palabos data processors has some bug when applying ReductiveDotProcessor
// some statistics evaluate is missing

template<typename T, size_t dim>
class DotFindMin : public PlbTypeDeduction::GetReductiveDotProcessingFunctionalXD_S<T,dim>::value
{
public:
    DotFindMin()
     : maxScalarId(this->getStatistics().subscribeMax())
    {

    }

    virtual DotFindMin<T,dim>* clone() const override
    {
        return new DotFindMin<T,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_SS<T,T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using Field = typename Traits::template argument<2>::type;

public:
    virtual void process(const DotList& dotList, Field& field) override
    {
        plb::BlockStatistics& statistics = this->getStatistics();

        for(const auto& dot : dotList.dots)
        {
            const T value = DataAccess::get(field, dot);
            // BlockStatistics computes only maximum, no minimum. Therefore,
            //   the relation min(x) = -max(-x) is used.
            statistics.gatherMax(maxScalarId, -value);
        }

    }


    T getMinScalar() const
    {
        // The minus sign accounts for the relation min(x) = -max(-x).
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        const plb::BlockStatistics& statistics = this->getStatistics();
        double doubleMin = - statistics.getMax(maxScalarId);
        if (std::numeric_limits<T>::is_integer)
        {
            return (T) plb::util::roundToInt(doubleMin);
        }
        return (T) doubleMin;
    }


private:
    plb::plint maxScalarId;
};

}

#endif // DOTFINDMIN_H
