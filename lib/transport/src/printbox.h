#ifndef PRINTBOX_H
#define PRINTBOX_H

#include "plbtypededuction.h"
#include "palabosio.h"
#include <sstream>

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class PrintBox
        : public PlbTypeDeduction::GetBoxProcessingFunctionalXD_L<T,Descriptor>::value
{
public:
    PrintBox()
    {

    }

    virtual PrintBox<T,Descriptor>* clone() const
    {
        return new PrintBox<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_L<T,Descriptor>::value;

    using Traits =
    PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;

    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;

    virtual void process(Box domain, Lattice& lattice)
    {
        Box box = lattice.getBoundingBox();
        std::stringstream ss;
        ss << domain << "\t->\t" << box << std::endl;
        std::cerr << ss.str();
        assert(plb::contained(domain, box));
    }


};

}

#endif // PRINTBOX_H
