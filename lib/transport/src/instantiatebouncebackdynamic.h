#ifndef INSTANTIATEBOUNCEBACKDYNAMIC_H
#define INSTANTIATEBOUNCEBACKDYNAMIC_H

#include "plbtypededuction.h"

#include <palabos/core/dynamics.h>

#include "palabosdataaccess.h"




namespace IPP
{

template<size_t dim>
struct AttributeDynamics;

template<>
struct AttributeDynamics<2>
{
    template<typename T, template<typename U> class Descriptor, typename Lattice>
    static void apply(Lattice& lattice, const plb::Dot2D& dot, plb::Dynamics<T,Descriptor>* dynamic)
    {
        lattice.attributeDynamics(dot.x, dot.y, dynamic);
    }

};

template<>
struct AttributeDynamics<3>
{
    template<typename T, template<typename U> class Descriptor, typename Lattice>
    static void apply(Lattice& lattice, const plb::Dot3D& dot,
                      plb::Dynamics<T,Descriptor>* dynamic)
    {
        lattice.attributeDynamics(dot.x, dot.y, dot.z, dynamic);
    }
};

template<typename T, template<typename U> class Descriptor>
class InstantiateBounceBackDynamic_LS : public PlbTypeDeduction::
        GetDotProcessingFunctionalXD_LS<T, Descriptor>::value
{
public:
    InstantiateBounceBackDynamic_LS()
    {

    }

    virtual ~InstantiateBounceBackDynamic_LS()
    {

    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetDotProcessingFunctionalXD_LS<T, Descriptor>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using DotList = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;
    using Field = typename Traits::template argument<3>::type;


public:
    virtual void process(DotList& dotList, Lattice& lattice, Field& field)
    {
        for(plb::plint iDot = 0; iDot < dotList.getN(); ++iDot)
        {
            const auto& dot = dotList.getDot(iDot);
            const T& rho = DataAccess::get(field, dot);

            AttributeDynamics<Descriptor<T>::d>::apply(lattice, dot, new plb::BounceBack<T, Descriptor>(rho));
        }
    }

    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        // Dynamics needs to be instantiated everywhere, including envelope??
        // see getTypeOfModification
        return plb::BlockDomain::bulkAndEnvelope;
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        // TODO: why not plb::modif::dataStructure?
        modified[0] = plb::modif::allVariables;
        modified[1] = plb::modif::nothing;
    }

    virtual InstantiateBounceBackDynamic_LS<T,Descriptor>* clone() const
    {
        return new InstantiateBounceBackDynamic_LS<T,Descriptor>(*this);
    }


};



}

#endif // INSTANTIATEBOUNCEBACKDYNAMIC_H
