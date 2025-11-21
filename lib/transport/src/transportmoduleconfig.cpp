#include "transportmoduleconfig.h"

namespace IPP
{

TransportModuleConfig::TransportModuleConfig(const boost::mpi::communicator& globalComm,
                                             const IPPConfig& conf,
                                             const boost::filesystem::path &resultDir,
                                             const boost::filesystem::path& checkpointDir,
                                             SimulationExchangeData& simData,
                                             const double &convergenceTolerance, bool &isConverged)
    : globalComm(globalComm),
      nx(conf.nx),
      ny(conf.ny),
      nz(conf.nz),
      is3DSimulation(conf.is3DSimulation),
      spatialResolution(conf.spatialResolution),
      diffusionCoefRef(conf.diffusionCoefReference),
      tauRef(conf.tauRef),
      porosRef(conf.porosRef),
      porosLow(conf.porosLow),
      latticeBoltzmannCollType(conf.latticeBoltzmannCollType),
      offset(conf.offset),
      flowFunc(conf.flowFunc),
      resultDir(resultDir),
      checkpointDir(checkpointDir)
    , distTransPorosThresh(conf.distTransPorosThresh)
    , isChemistryEnabled(conf.isChemistryEnabled)
    , ippConf(conf)
    , simData(simData)
    , convergenceTolerance(convergenceTolerance)
    , isConverged(isConverged)
    , resultProcessing(conf.resultProcessing)
{
    for(BoundaryConditionData& data : conf.boundaryConditions->advectiveBC)
    {
        bcData.push_back(&data);
    }

    for(BoundaryConditionData& data : conf.boundaryConditions->diffusiveBC)
    {
        bcData.push_back(&data);
    }

    for(const IPPConfig::Results::DiffEffInfo& info : conf.results.diffEffInfos)
    {
        effDiffDims.push_back(info.dim);
    }
}

} // end of namespace LBGeoChem
