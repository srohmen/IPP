#ifndef DENSITYVELOCITYRESULTS_H
#define DENSITYVELOCITYRESULTS_H

#include <memory>
#include <vector>

namespace IPP
{

namespace IPPResults
{

    template<typename TransportTraits>
    struct DensityVelocityResults
    {
        typename TransportTraits::ScalarFieldSharedPtr density;
        typename TransportTraits::TensorFieldSharedPtr velocity;
    };

    template<typename TransportTraits>
    using DensityVelocityResultsVec = std::vector< DensityVelocityResults<TransportTraits> >;

    template<typename TransportTraits>
    struct FlowDiffResults
    {
        DensityVelocityResults<TransportTraits> flowResults;
        DensityVelocityResultsVec<TransportTraits> diffResultsVec;
    };

} // end of namespace LBGeoChemResults

}

#endif // DENSITYVELOCITYRESULTS_H
