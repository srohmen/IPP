#ifndef PRINTPOROSITY_H
#define PRINTPOROSITY_H

#include <palabos/atomicBlock/dataProcessingFunctional2D.h>
#include <palabos/core/cell.h>


namespace IPP
{

template<typename T, template<typename U> class Descriptor>
struct PrintPorosity : public plb::BoxProcessingFunctional2D_L<T, Descriptor>
{

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
    }

    virtual PrintPorosity<T,Descriptor>* clone() const
    {
        return new PrintPorosity<T,Descriptor>(*this);
    }

    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T, Descriptor>& lattice)
    {
        for(int y = 0; y < lattice.getNy(); ++y)
        {
            for(int x = 0; x < lattice.getNx(); ++x)
            {
                const plb::Cell<T, Descriptor>& cell = lattice.get(x, y);
                std::cout << *cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt) << "\t";
            }

            std::cout << std::endl;
        }
    }
};


}

#endif // PRINTPOROSITY_H
