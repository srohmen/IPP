#ifndef PRINTSCALARS_H
#define PRINTSCALARS_H

#include "plbtypededuction.h"

#include "mpimanager.h"
#include "domainiterator.h"
#include <bitset>

namespace IPP
{

template<typename T, size_t dim>
struct PrintScalars : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_S<T, dim>::value
{

public:
    PrintScalars()
    {

    }

    virtual PrintScalars<T,dim>* clone() const override
    {
        return new PrintScalars<T,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }


private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, InputField& field) override
    {
        std::stringstream ss;
        ss << "-------------------------\n";
        ss << "rank: " << MPIManager::getInstance().getRank() << std::endl;

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto pos = *it;
            const T& value = DataAccess::get(field, pos);

            std::bitset< sizeof(T)*8 > binary(value);

            ss << pos << ": " << value << "\t" << binary << std::endl;;
        }

        ss << "-------------------------\n";
        std::cout << ss.str() << std::flush;
    }
};


}

#endif // PRINTSCALARS_H
