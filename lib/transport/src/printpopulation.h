#ifndef PRINTPOPULATION_H
#define PRINTPOPULATION_H

#include "plbtypededuction.h"

namespace IPP
{

template<typename T, template<typename U> class Descriptor>
class PrintPopulation
        : public PlbTypeDeduction::GetOneCellIndexedFunctionalXD<T,Descriptor>::value
{
public:
    PrintPopulation()
    {

    }

    virtual PrintPopulation<T,Descriptor>* clone() const
    {
        return new PrintPopulation<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return plb::BlockDomain::bulkAndEnvelope;
    }


private:
    virtual void execute(plb::plint x, plb::plint y,
                         plb::Cell<T,Descriptor>& cell) const
    {
        std::stringstream ss;
        ss << plb::global::mpi().getRank()
           << " (" << x << "," << y << "): ";

        for(size_t iPop = 0; iPop < Descriptor<T>::q; ++iPop)
        {
            ss << cell[iPop] << " ";
        }
        ss << std::endl;

        std::cout << ss.str();

    }

    virtual void execute(plb::plint x, plb::plint y, plb::plint z,
                         plb::Cell<T,Descriptor>& cell) const
    {
        std::stringstream ss;
        ss << plb::global::mpi().getRank()
           << " (" << x << "," << y << "," << z << "): ";

        for(size_t iPop = 0; iPop < Descriptor<T>::q; ++iPop)
        {
            ss << cell[iPop] << " ";
        }
        ss << std::endl;

        std::cout << ss.str();

    }


};

}

#endif // PRINTPOPULATION_H
