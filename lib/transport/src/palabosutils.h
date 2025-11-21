#ifndef PALABOSUTILS_H
#define PALABOSUTILS_H

#include "arraydimensionconvert.h"
#include "palabosboundaryconditionutils.h"
#include "palabosinitfunctors.h"
#include "palabosconversiontools.h"
#include "abstractboundaryconditions.h"
#include "palabosmpiutils.h"
#include "mpitools.h"
#include "palabosobjectfactory.h"
#include "mathtools.h"
#include "palaboslatticevalueaccess.h"
#include "array_view.h"
#include "serialtomulticonversion.h"
#include "palabosdataaccess.h"

// TODO: check includes
#include <palabos/multiBlock/multiDataProcessorWrapper2D.hh>
#include <palabos/multiBlock/multiDataProcessorWrapper3D.hh>

#include <palabos/dataProcessors/dataInitializerFunctional3D.h>

#include <palabos/dataProcessors/dataAnalysisWrapper2D.h>


namespace IPP
{

template<typename T, template<class U> class Descriptor>
class PrintDensity : public plb::OneCellIndexedFunctional2D<T, Descriptor>
{
public:

    virtual plb::OneCellIndexedFunctional2D<T,Descriptor>* clone() const
    {
        return new PrintDensity<T, Descriptor>(*this);
    }

    virtual void execute(plb::plint iX, plb::plint iY, plb::Cell<T,Descriptor>& cell) const
    {
        plb::pcout << iX << "\t" << iY << "\t" << cell.computeDensity() << std::endl;
    }
};

template<typename TransportTraits>
struct ConcSetUpdate
{
    using Descriptor = typename TransportTraits::DiffusionDescriptor;
    using Lattice = typename TransportTraits::DiffusionLattice;
    using ScalarField = typename TransportTraits::ScalarField;
    using ScalarFieldPtr = typename TransportTraits::ScalarFieldPtr;

    void operator()(ScalarField& concentrations,
                    ScalarField& preChemistryConc,
                    Lattice& diffLattice) const
    {
        ScalarFieldPtr difference = plb::subtract(concentrations, preChemistryConc);
        PalabosLatticeValueAccess::setScalarIfPossible<Descriptor>(diffLattice, *difference);
    }

};


}

#endif // PALABOSUTILS_H
