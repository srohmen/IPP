#ifndef TRANSPORTMODULECONFIG_H
#define TRANSPORTMODULECONFIG_H

#include <boost/filesystem/path.hpp>
#include "ippconfig.h"

namespace IPP
{

class SimulationExchangeData;

class TransportModuleConfig
{
public:
    TransportModuleConfig(const boost::mpi::communicator& globalComm,
                          const IPPConfig& conf,
                          const boost::filesystem::path& resultDir,
                          const boost::filesystem::path& checkpointDir,
                          SimulationExchangeData& simData,
                          const double& convergenceTolerance,
                          bool& isConverged);


    const boost::mpi::communicator& globalComm;

    const size_t nx, ny, nz;
    const bool is3DSimulation;
    const double spatialResolution;
    const double diffusionCoefRef;
    const double tauRef;
    const double porosRef;
    const double porosLow;
    const LatticeBoltzmannCollisionType latticeBoltzmannCollType;

    const IPPVector3D& offset;
    const AbstractFlowFunctorPtr& flowFunc;
    const boost::filesystem::path resultDir;
    const boost::filesystem::path checkpointDir;
    const double distTransPorosThresh;

    const bool isChemistryEnabled;

    const IPPConfig& ippConf;

    std::vector<BoundaryConditionData*> bcData;
    std::vector<size_t> effDiffDims;
    SimulationExchangeData& simData;

    const double& convergenceTolerance;
    bool& isConverged;

    std::shared_ptr<AuxResultProcessing> resultProcessing;

};



} // end of namespace LBGeoChem

#endif // TRANSPORTMODULECONFIG_H
