#ifndef GETLOCALBOUNDINGBOX_H
#define GETLOCALBOUNDINGBOX_H

#include "plbtypededuction.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class GetLocalBoundingBox
        : public PlbTypeDeduction::GetBoxProcessingFunctionalXD_L<T,Descriptor>::value
{
private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_L<T,Descriptor>::value;

    using Traits =
    PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;

    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;


public:
    GetLocalBoundingBox(Box& localBox)
        : m_box(localBox)
    {

    }

    virtual GetLocalBoundingBox<T,Descriptor>* clone() const
    {
        return new GetLocalBoundingBox<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
    }

private:
    virtual void process(Box domain, Lattice& lattice)
    {
        m_box = lattice.getBoundingBox();
    }

private:
    Box& m_box;
};

}

#endif // GETLOCALBOUNDINGBOX_H
