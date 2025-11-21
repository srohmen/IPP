#ifndef PALABOSRESULTWRITER_H
#define PALABOSRESULTWRITER_H

#include "palabosobjectfactory.h"
#include "palabostransportdata.h"
#include "densityvelocityresults.h"
#include "PalabosConvergenceCheck.h"

namespace IPP
{

namespace ResultWriter
{

template<typename TransportTraits>
class PalabosResultWriter
{
private:

    using Scalar = typename TransportTraits::Scalar;

    using ScalarField = typename TransportTraits::ScalarField;
    using ScalarFieldPtr = typename TransportTraits::ScalarFieldPtr;
    using TensorField = typename TransportTraits::TensorField;
    using TensorFieldPtr = typename TransportTraits::TensorFieldPtr;
    using DiffusionDescriptor = typename TransportTraits::DiffusionDescriptor;
    using DiffusionLattice = typename TransportTraits::DiffusionLattice;
    using DiffusionLatticePtr = typename TransportTraits::DiffusionLatticePtr;
    using Vector = typename TransportTraits::Vector;
    using Box = typename TransportTraits::Box;
    using Dot = typename TransportTraits::Dot;
    static constexpr size_t dim = TransportTraits::dim;
    using PlbFactory = PalabosObjectFactory<dim>;


    using PlbTransData = PalabosTransportData<TransportTraits>;
    using FlowDiffResults = IPPResults::FlowDiffResults<TransportTraits>;

    using VtkOut = typename TransportTraits::VtkOut;


    using ConvCheck = PalabosConvergenceCheck<Scalar, dim>;

public:

    static void writeDebugData(const ResultsToWrite& resultsToWrite,
                               const Scalar& dx,
                               const plb::Array<Scalar,3>& offset3d,
                               const std::array<size_t,3>& dims,
                               const SimulationExchangeData& simData,
                               const PlbTransData& palabosData,
                               const PlbFactory& factory);




    static void writeResults(const ResultsToWrite& resultsToWrite,
                             const Scalar& dx,
                             const plb::Array<Scalar,3>& offset3d,
                             const std::array<size_t,3>& dims,
                             const FlowDiffResults& results,
                             const SimulationExchangeData &simData,
                             const PlbTransData &palabosData,
                             const PlbFactory& factory,
                             std::shared_ptr<AuxResultProcessing>& resultProcess,
                             ConvCheck &convergenceCheck);

    static void writeOutletFlux(const Scalar& dx,
                                const plb::Array<Scalar,3>& offset3d,
                                const std::array<size_t,3>& dims,
                                const SimulationExchangeData &simData,
                                const PlbTransData &palabosData,
                                const PlbFactory& factory,
                                const bool writeFluxImage,
                                ConvCheck &convergenceCheck);
};


}
}

#include "externtemplatehelper.h"
#ifdef EXPLICIT_INSTANTS
EXTERN_TEMPLATE_TRANSPORT_TRAITS(class IPP::ResultWriter::PalabosResultWriter)
#else
#include "palabosresultwriter.hh"
#endif

#endif // PALABOSRESULTWRITER_H
