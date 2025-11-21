#ifndef CONCSETINIT_H
#define CONCSETINIT_H

namespace IPP
{

template<typename PorosityFunc>
class ConcSetInit
{
public:    
    ConcSetInit(const PorosityFunc& func)
        : m_func(func)
    {

    }

    template<typename T,
             template<typename U> class Descriptor,
             template<typename U, template<typename V> class Desc> class Lattice,
             typename ScalarField
             >
    void operator()(ScalarField& concentrations,
                    const ScalarField& /*preChemistryConc*/,
                    Lattice<T,Descriptor>& diffLattice) const
    {      
        m_func.apply(diffLattice, concentrations);
    }

private:
    const PorosityFunc& m_func;

};

}

#endif // CONCSETINIT_H
