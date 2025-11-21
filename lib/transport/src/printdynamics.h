#ifndef PRINTDYNAMICS_H
#define PRINTDYNAMICS_H

#include "plbtypededuction.h"
#include "palabosdataaccess.h"
#include "mpimanager.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class PrintDynamics : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_L<T,Descriptor>::value
{
public:
    PrintDynamics() = default;

    virtual PrintDynamics<T,Descriptor>* clone() const override
    {
        return new PrintDynamics<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulkAndEnvelope;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_L<T,Descriptor>::value;

    using Traits =
    PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;

    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, Lattice& lattice) override
    {
        std::stringstream ss;
        ss << "-------------------------\n";
        ss << "rank: " << MPIManager::getInstance().getRank() << std::endl;

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto pos = *it;
            const plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, pos);
            const plb::Dynamics<T,Descriptor>& dyn = cell.getDynamics();
            ss << pos << ": " << dyn.getId() << std::endl;;
        }

        ss << "-------------------------\n";
        std::cout << ss.str() << std::flush;

    }


};


}

#endif // PRINTDYNAMICS_H
