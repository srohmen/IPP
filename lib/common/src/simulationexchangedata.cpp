#include "simulationexchangedata.h"

#include "mpitools.h"
#include "ippexception.h"

#include "mpimanager.h"
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

namespace IPP
{


SimulationExchangeData::SimulationExchangeData(const boost::mpi::communicator& comm)
    : m_comm(comm),
      m_iteration(0),
      m_time(0),
      m_dtCurr(-1.0),
      m_decomp(nullptr)
{

}

void SimulationExchangeData::setBoundaryConditions(const AbstractBoundaryConditionsPtr&bc)
{
    m_boundaryConditions = bc;
}

const std::vector<std::string>& SimulationExchangeData::getCompNames() const
{
    return m_compNames;
}

const AbstractBoundaryConditionsPtr& SimulationExchangeData::getBoundaryConditions() const
{
    return m_boundaryConditions;
}

std::vector<double>& SimulationExchangeData::getInertComposition()
{
    return m_inertComposition;
}

std::vector<double> &SimulationExchangeData::getInertVolFrac()
{
    return m_localInertVolFrac;
}

const std::vector<double> &SimulationExchangeData::getInertVolFrac() const
{
    return m_localInertVolFrac;
}

std::vector<char> &SimulationExchangeData::getInertSolidCells()
{
    return m_inertSolidCells;
}

void SimulationExchangeData::setCompNames(const std::vector<std::string>& names)
{
    m_compNames = names;
}

std::vector<double>& SimulationExchangeData::getPostReacConc()
{
    // MPI_CHECK_SYNC;
    return m_postReacConc;
}

const std::vector<double>& SimulationExchangeData::getPostReacConc() const
{
    // MPI_CHECK_SYNC;
    return m_postReacConc;
}

std::vector<double>& SimulationExchangeData::getPostTransportConc()
{
    // MPI_CHECK_SYNC;
    return m_postTransConc;
}

const std::vector<double>& SimulationExchangeData::getPostTransportConc() const
{
    // MPI_CHECK_SYNC;
    return m_postTransConc;
}

std::vector<double> &SimulationExchangeData::getPrecipField()
{
    // MPI_CHECK_SYNC
    return m_precipField;
}

const std::vector<double> &SimulationExchangeData::getPrecipField() const
{
    // MPI_CHECK_SYNC
    return m_precipField;
}

std::vector<double>& SimulationExchangeData::getPorosity()
{
    return m_porosityPhysical;
}

const std::vector<double> &SimulationExchangeData::getPorosity() const
{
    return m_porosityPhysical;
}

std::vector<double> &SimulationExchangeData::getCapillaryPorosity()
{
    return m_porosityCapillary;
}

const std::vector<double> &SimulationExchangeData::getCapillaryPorosity() const
{
    return m_porosityCapillary;
}

std::vector<double> &SimulationExchangeData::getDiffusionCoefs()
{
    return m_localDiffusionCoef;
}

const std::vector<double> &SimulationExchangeData::getDiffusionCoefs() const
{
    return m_localDiffusionCoef;
}

std::vector<double> &SimulationExchangeData::getDistField()
{
    return m_distField;
}

const std::vector<double> &SimulationExchangeData::getDistField() const
{
    return m_distField;
}

std::vector<std::string>& SimulationExchangeData::getPhaseNames()
{
    return m_phaseNames;
}

const std::vector<std::string>& SimulationExchangeData::getPhaseNames() const
{
    return m_phaseNames;
}

std::vector<std::string> &SimulationExchangeData::getEqPhaseNames()
{
    return m_eqPhaseNames;
}

const std::vector<std::string>&SimulationExchangeData::getEqPhaseNames() const
{
    return m_eqPhaseNames;
}

std::vector<double>& SimulationExchangeData::getSaturationIndices()
{
    return m_saturationIndices;
}

const std::vector<double>& SimulationExchangeData::getSaturationIndices() const
{
    return m_saturationIndices;
}

AuxDataVec &SimulationExchangeData::getAuxData()
{
    return m_auxData;
}

const AuxDataVec&SimulationExchangeData::getAuxData() const
{
    return m_auxData;
}

void SimulationExchangeData::setIteration(const size_t iteration)
{
    m_iteration = iteration;
}

size_t SimulationExchangeData::getIteration() const
{
    return m_iteration;
}

void SimulationExchangeData::initTime(const double &time)
{
    m_time = time;
}

void SimulationExchangeData::incrementTime()
{
    m_time += m_dtCurr;
}

double SimulationExchangeData::getTime() const
{
    return m_time;
}

void SimulationExchangeData::setInertTracer(const size_t diffDim, const std::string& compName)
{
    auto it = std::find(m_compNames.begin(), m_compNames.end(), compName);
    IPPCheck::assertCheck(it != m_compNames.end());
    const size_t iComp = it - m_compNames.begin();
    m_dimToDiffSpecies.push_back(std::make_pair(diffDim, iComp));
}

void SimulationExchangeData::setDtCurr(const double &dt)
{
    m_dtCurr = dt;
}

const double& SimulationExchangeData::getDtCurr() const
{
    return m_dtCurr;
}

const FieldDecomposition *SimulationExchangeData::getDecomp() const
{
    return m_decomp;
}

void SimulationExchangeData::setDecomp(const FieldDecomposition *decomp)
{
    m_decomp = decomp;
}

PhaseNameToInfos &SimulationExchangeData::getPhaseNameToInfos()
{
    return m_phaseNameToInfos;
}

const PhaseNameToInfos& SimulationExchangeData::getPhaseNameToInfos() const
{
    return m_phaseNameToInfos;
}

const SimulationExchangeData::DimToDiffSpecies &SimulationExchangeData::getDimToDiffSpecies() const
{
    return m_dimToDiffSpecies;
}

void SimulationExchangeData::initSyncEarly()
{
    MpiTools::syncFromRoot(m_comm, m_compNames);
    MpiTools::syncFromRoot(m_comm, m_phaseNames);
    MpiTools::syncFromRoot(m_comm, m_eqPhaseNames);
    boost::mpi::broadcast(m_comm, m_phaseNameToInfos, 0);

    MPI_CHECK_SYNC;
}

void SimulationExchangeData::initSyncLate()
{
    MpiTools::syncFromRoot(m_comm, m_dimToDiffSpecies);
    MPI_CHECK_SYNC;
}


}
