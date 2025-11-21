#ifndef RESULTCALC_H
#define RESULTCALC_H

#include <functional>
#include <numeric>

#include "ippexception.h"
#include "densityvelocityresults.h"
#include "palabostransportdata.h"


namespace IPP
{

namespace IPPResults
{

template<typename TransportTraits>
struct ResultCalc
{
    typedef typename TransportTraits::Scalar T;

    typedef IPPResults::DensityVelocityResults<TransportTraits> DensityVelocityResults;
    typedef IPPResults::DensityVelocityResultsVec<TransportTraits> DensityVelocityResultsVec;
    typedef IPPResults::FlowDiffResults<TransportTraits> FlowDiffResults;

    typedef PalabosTransportData<TransportTraits> PlbTransData;

    template<typename LatticeType, typename TransFunc>
    static void calcResultFields(LatticeType& lattice,
                                 const TransFunc& transformation,
                                 DensityVelocityResults& results)
    {
        const auto box = lattice.getBoundingBox();
        plb::computeDensity(lattice, *results.density, box);
        transformation.apply(*results.density);
        plb::computeVelocity(lattice, *results.velocity, box);
    }

    template<typename LatticeType>
    static void calcResultFields(LatticeType& lattice,
                                 DensityVelocityResults& results)
    {
        const auto box = lattice.getBoundingBox();
        plb::computeDensity(lattice, *results.density, box);
        plb::computeVelocity(lattice, *results.velocity, box);
    }

    static void calc(const PlbTransData& palabosData, FlowDiffResults& result)
    {
        DensityVelocityResults& flowResults = result.flowResults;
        calcResultFields(*(palabosData.flowLattice), flowResults);

        DensityVelocityResultsVec& diffResultsVec = result.diffResultsVec;
        const size_t nDiffLattices = palabosData.diffLattices.size();
        IPPCheck::assertCheck(nDiffLattices == diffResultsVec.size());

        for(size_t iComp = 0; iComp < nDiffLattices; ++iComp)
        {
            const typename TransportTraits::DiffusionLatticePtr& diffLattice =
                    palabosData.diffLattices[iComp];

            DensityVelocityResults& diffResults = diffResultsVec[iComp];
            calcResultFields(*diffLattice, diffResults);
        }
    }
};

} // end of namespace LBGeoChemResults

}

#endif // RESULTCALC_H
