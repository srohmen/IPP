#ifndef PALABOSTRANSPORTMODULE_H
#define PALABOSTRANSPORTMODULE_H

#include "transportmodule.h"

#include "transportmoduleconfigfwd.h"
#include "palabosobjectfactory.h"
#include "arraydimensionconvert.h"
#include "abstractboundaryconditions.h"
#include "nonlocaloperations.h"
#include "PalabosConvergenceCheck.h"

// TODO: obsolete! remove flow results calc from transport module
// calc results on demand in result writer!
#include "resultcalc.h"


namespace IPP
{

template<typename TransportTraits>
class PalabosTransportData;

class SimulationExchangeData;

template<typename TransportTraits>
class PalabosTransportModule : public TransportModule
{
public:
    PalabosTransportModule(TransportModuleConfigPtr& config);
    virtual ~PalabosTransportModule();

    virtual void init() override final;
    virtual NonLocalOperations& getNonLocalOperations() override final;
    virtual const DomainInfos &getDomainInfos() const override final;
    virtual const FieldDecomposition* getDecomposition() const override final;


    virtual void run() override final;


    virtual void writeDebugData(const ResultsToWrite &resultsToWrite) const override final;
    virtual void writeResults(const ResultsToWrite& resultsToWrite) override final;
    virtual void writeOutletFlux() override final;


    virtual void updateRenderer(IPPRendererPtr& renderer) override;

    virtual void updateTransportScalarDref(const double& newPorosRef) override;

    
protected:
    using Scalar = typename TransportTraits::Scalar;
    static constexpr size_t dim = TransportTraits::dim;
    using Dot = typename TransportTraits::Dot;
    using DotList = typename TransportTraits::DotList;
    using Box = typename TransportTraits::Box;

    using ScalarFieldPtr = typename TransportTraits::ScalarFieldPtr;
    using ScalarField = typename TransportTraits::ScalarField;
    using ScalarFieldSharedPtr = typename TransportTraits::ScalarFieldSharedPtr;

    using TensorField = typename TransportTraits::TensorField;
    using TensorFieldPtr = typename TransportTraits::TensorFieldPtr;


    //    template <typename T>
    //    using HydrodynamicDescT = typename TransportTraits:: template HydrodynamicDescriptorT<T>;
    //    using HydrodynamicDesc = typename TransportTraits::HydrodynamicDescriptor;
    using HydrodynamicLattice = typename TransportTraits::HydrodynamicLattice;
    using HydrodynamicLatticePtr = typename TransportTraits::HydrodynamicLatticePtr;
    using HydrodynamicDynamic = typename TransportTraits::HydrodynamicDynamic;


    //    template <typename T>
    //    using DiffusionDescT = typename TransportTraits:: template DiffusionDescriptorT<T>;
    //    using DiffusionDesc = typename TransportTraits::DiffusionDescriptor;
    using DiffusionLattice = typename TransportTraits::DiffusionLattice;
    using DiffusionLatticePtr = typename TransportTraits::DiffusionLatticePtr;
    using DiffusionDynamic = typename TransportTraits::DiffusionDynamic;


    using PlbFactory = PalabosObjectFactory<dim>;
    using TransportData = PalabosTransportData<TransportTraits>;
    using ArrayView = av::array_view<Scalar, dim>;


    void initHydrodynamicLattice(const std::vector<int>& periodic,
                                 const AbstractBoundaryConditions::BCVec& advectiveBC);
    void initDiffusionLattices(const AbstractBoundaryConditions& bc);
    void initInertCellsMarker();
    void finishLatticeInit();


    void saveLatticeCheckpoint() const;
    void loadLatticeCheckpoint(size_t iteration);

    void updateFieldTransportScalar();

    using ResultValCalc = typename IPPResults::ResultCalc<TransportTraits>;
    using FlowDiffResults = typename ResultValCalc::FlowDiffResults;
    void prepareResultArrays(typename ResultValCalc::DensityVelocityResults& results);


    SimulationExchangeData& m_simData;
    TransportData m_palabosData;

    const TransportModuleConfigPtr m_conf;
    size_t nx, ny, nz;
    DomainInfos m_domainInfos;
    std::vector<int> m_periodicity;
    Box m_transportSyncEnabledBox;

    std::unique_ptr<FieldDecomposition> m_decomp;
    std::unique_ptr<PlbFactory> m_plbFactory;

    NonLocalOperations m_nonLocalOperations;

    std::vector<Scalar> m_transportScalar;

    FlowDiffResults m_results;

private:
    virtual void initModule() = 0;
    virtual void addPostStreamOperation(DiffusionLattice &lattice) = 0;
    virtual void prepareResults() = 0;


    void initSize();

    template<typename Transformation>
    void setupDiffBoundaryConditions(const AbstractBoundaryConditions& bc,
                                     const std::string& speciesName,
                                     const Transformation& trans,
                                     DiffusionLattice& lattice);


    void updateLatticeTransportProperties();

    const plb::Array<Scalar,3> m_offset;

    using OmegaCalc = typename TransportTraits::template DiffusionOmegaCalc<Scalar>;
    OmegaCalc m_omegaCalc;


    PalabosConvergenceCheck<Scalar, dim> m_convergenceCheck;

};


} // end of namespace IPP


#include "externtemplatehelper.h"
EXTERN_TEMPLATE_TRANSPORT_TRAITS(class IPP::PalabosTransportModule)

#ifndef EXPLICIT_INSTANTS
#include "palabostransportmodule.hh"
#endif

#endif // PALABOSTRANSPORTMODULE_H
