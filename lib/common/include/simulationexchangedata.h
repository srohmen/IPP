#ifndef SIMULATIONEXCHANGEDATA_H
#define SIMULATIONEXCHANGEDATA_H

#include <vector>
#include <string>
#include <boost/mpi/communicator.hpp>

#include "abstractboundaryconditionsfwd.h"
#include "phasenametoinfos.h"
#include "ippvector.h"
#include "auxdata.h"
#include "celltotalsdiff.h"
#include "cellneighborinfo.h"
#include "debugdata.h"

namespace IPP
{

class FieldDecomposition;

class SimulationExchangeData
{
public:
    SimulationExchangeData(const boost::mpi::communicator& comm);

    void setBoundaryConditions(const AbstractBoundaryConditionsPtr& bc);
    const AbstractBoundaryConditionsPtr& getBoundaryConditions() const;

    // TODO: add const getter ?
    std::vector<double>& getInertComposition();
    std::vector<char> &getInertSolidCells();

    void setCompNames(const std::vector<std::string>& names);
    const std::vector<std::string>& getCompNames() const;

    std::vector<double>& getPostReacConc();
    const std::vector<double>& getPostReacConc() const;

    std::vector<double>& getPostTransportConc();
    const std::vector<double>& getPostTransportConc() const;

    std::vector<double>& getPrecipField();
    const std::vector<double>& getPrecipField() const;

    std::vector<double>& getInertVolFrac();
    const std::vector<double>& getInertVolFrac() const;

    std::vector<double>& getPorosity();
    const std::vector<double>& getPorosity() const;

    std::vector<double>& getCapillaryPorosity();
    const std::vector<double>& getCapillaryPorosity() const;

    std::vector<double>& getDiffusionCoefs();
    const std::vector<double>& getDiffusionCoefs() const;

    std::vector<double>& getDistField();
    const std::vector<double>& getDistField() const;


    std::vector<std::string>& getPhaseNames();
    const std::vector<std::string>& getPhaseNames() const;

    std::vector<std::string>& getEqPhaseNames();
    const std::vector<std::string>& getEqPhaseNames() const;

    std::vector<double>& getSaturationIndices();
    const std::vector<double>& getSaturationIndices() const;


    AuxDataVec& getAuxData();
    const AuxDataVec& getAuxData() const;

    void setIteration(const size_t iteration);
    size_t getIteration() const;

    void initTime(const double& time);
    void incrementTime();
    double getTime() const;

    // need to be executed once at the right time
    // TODO: make more clear what is synced and why
    // with changing function names
    void initSyncEarly();
    void initSyncLate();
    /////////////////////

    typedef std::vector<std::pair<size_t, size_t>> DimToDiffSpecies;
    const DimToDiffSpecies& getDimToDiffSpecies() const;

    void setInertTracer(const size_t diffDim, const std::string &compName);

    void setDtCurr(const double& m_dtCurr);
    const double& getDtCurr() const;

    const FieldDecomposition *getDecomp() const;
    void setDecomp(const FieldDecomposition *decomp);

    PhaseNameToInfos& getPhaseNameToInfos();
    const PhaseNameToInfos& getPhaseNameToInfos() const;



    std::vector<CellTotalsDiff> additionalSources;

    std::vector<double> solidFlags;
    std::vector<std::vector<double>> nucleationPhasesVolumeFractions;
    std::vector<CellNeighborInfo> neighInfos;

#ifdef IPP_DEBUG
    DebugData debugData;
#endif



private:
    // global data: update by each proc
    const boost::mpi::communicator& m_comm;
    size_t m_iteration;
    double m_time;
    double m_dtCurr;
    AbstractBoundaryConditionsPtr m_boundaryConditions;



    // local data
    const FieldDecomposition* m_decomp;
    std::vector<double> m_postReacConc;
    std::vector<double> m_postTransConc;

    std::vector<double> m_precipField;

    std::vector<double> m_localInertVolFrac;
    std::vector<double> m_porosityPhysical;
    std::vector<double> m_porosityCapillary;

    std::vector<double> m_localDiffusionCoef;
    std::vector<double> m_distField;

    std::vector<char> m_inertSolidCells;
    AuxDataVec m_auxData;


    // TODO: currently main proc data, move to local data
    std::vector<double> m_saturationIndices;

    // global data: must be initialized once on all children procs
    std::vector<std::string> m_compNames;
    std::vector<std::string> m_phaseNames;
    std::vector<std::string> m_eqPhaseNames;
    PhaseNameToInfos m_phaseNameToInfos;
    DimToDiffSpecies m_dimToDiffSpecies;

    // root proc only data
    std::vector<double> m_inertComposition;


};

}

#endif // SIMULATIONEXCHANGEDATA_H
