#ifndef REACTIVEPALABOSTRANSPORTMODULE_HH
#define REACTIVEPALABOSTRANSPORTMODULE_HH

#include "reactivepalabostransportmodule.h"

#include <unordered_map>

#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "ippconstants.h"
#include "fielddecomposition.h"
#include "palabosconversiontools.h"
#include "concarrayview.h"
#include "palaboslatticevalueaccess.h"
#include "plbcheckdata.h"
#include "countgeomchanged.h"
#include "calcdistance.h"
#include "simulationexchangedata.h"
#include "setdotstoscalar.h"
#include "configurebouncebackdynamics.h"
#include "dotsetexternalscalar.h"
#include "initequilibriumconcentration.h"

#include "palabosdistributeexcesstotals.h"
#include "palaboscalcinterfaceproperties.h"
#include "collectgeomchangedcoordinates.h"
#include "collectsolidcoordinateswithpermneighbor.h"

// #include "distributeinterfaceconc.h"
#include "porositybounceback.h"
#include "pushbouncebackpopulations.h"
#include "distributebouncebackpopulations.h"
#include "zerobouncebackpopulations.h"
#include "syncdataintolattice.h"
#include "getplbsize.h"
#include "multitoserialdatasync.h"
#include "serialtomultidatasync.h"
#include "plb_hash.h"
#include "createfilename.h"
#include "abstractgeometrysync.h"
#include "plbserialization.h"
#include "updatedistancefield.h"

#include "permflagsutils.h"
#include "setdottonsource.h"
#include "setntensortosource.h"
#include "setpermneighcounter.h"

#include "settontensorfunction.h"
#include "setntensorfromarray.h"

#include "cellneighborinfo.h"
#include "getneighinfo.h"

// debug helper
#include "printscalars.h"
#include "printdynamics.h"
/////////


#define HYDRODYNAMIC_DESCRIPTOR TransportTraits::template HydrodynamicDescriptorT
#define DIFFUSION_DESCRIPTOR TransportTraits::template DiffusionDescriptorT


namespace IPP
{



namespace
{



template<typename TensorField, size_t index, size_t dim>
class TensorWrapper;

template<typename TensorField, size_t index>
class TensorWrapper<TensorField, index, 2>
{
public:
    TensorWrapper(const TensorField& field)
        : m_field(field)
    {

    }

    template<typename I>
    auto& get(const I& x, const I& y)
    {
        return m_field.get(x, y)[index];
    }

    template<typename I>
    const auto& get(const I& x, const I& y) const
    {
        return m_field.get(x, y)[index];
    }

private:
    const TensorField& m_field;
};

template<typename TensorField, size_t index>
class TensorWrapper<TensorField, index, 3>
{
public:
    TensorWrapper(const TensorField& field)
        : m_field(field)
    {

    }

    template<typename I>
    auto& get(const I& x, const I& y, const I& z)
    {
        return m_field.get(x, y, z)[index];
    }

    template<typename I>
    const auto& get(const I& x, const I& y, const I& z) const
    {
        return m_field.get(x, y, z)[index];
    }

private:
    const TensorField& m_field;
};

template<typename T, size_t dim, size_t porosIndex, typename Dot, typename VolumesTensorN>
static void handleNeighVolumes(const Dot& center, const VolumesTensorN& volumeTensorField, CellNeighborInfo& result)
{
    const plb::plint nDim = volumeTensorField.getNdim();
    if(nDim > 1)
    {
        // need to calc vol fractions sum of neighs

        CellNeighborInfo::NuclPhaseVolumeFractions& fracsVec = result.phaseVolFrac;
        const size_t nAux = porosIndex + 1;
        const size_t nNuclPhases = nDim - nAux;
        fracsVec.resize(nNuclPhases, {0, 0} );

        using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;
        const DotList neighs = DotFindNeighbors::find(center);
        for(const Dot& neighPos : neighs.dots)
        {
            const T* arrVolFractions =  DataAccess::get(volumeTensorField, neighPos);

            const T& porosity = arrVolFractions[porosIndex];
            const T solidFrac = 1.0 - porosity;

            for(plb::plint d = porosIndex + 1; d < nDim; ++d)
            {
                const T& volFrac = arrVolFractions[d];

                const size_t iNuclPhase = d - nAux;
                assert(fracsVec.size() > iNuclPhase);
                CellNeighborInfo::VolumeFractions& fracs = fracsVec[iNuclPhase];

                fracs.self += static_cast<double>(volFrac);

                const double rest = solidFrac - volFrac;
                if(rest >= 0.0)
                {
                    fracs.rest += rest;
                }
                else if(rest < -1.0E-9)
                {
                    throw std::runtime_error("negative rest vol frac: " + std::to_string(rest));
                }
            }

        }

    }
}


template<typename T, typename Array, size_t dim, size_t solidFlags, size_t porosIndex>
class CalcNeighInfos : public
        PlbTypeDeduction::GetBoxProcessingFunctionalXD_N<T, dim>::value
{
public:
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

    CalcNeighInfos(Array& arr, const Dot& location, const T& thresh)
        : m_arr(arr)
        , m_location(location)
        , m_thresh(thresh)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulk;
    }

    virtual CalcNeighInfos<T, Array, dim, solidFlags, porosIndex>* clone() const override
    {
        return new CalcNeighInfos<T, Array, dim, solidFlags, porosIndex>(*this);
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_N<T, dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using VolumesTensorN = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, VolumesTensorN& volumeTensorField) override
    {
        const Dot location = volumeTensorField.getLocation();

        using FlagsWrapper = TensorWrapper<VolumesTensorN, solidFlags, dim>;
        const FlagsWrapper flagsField(volumeTensorField);

        using PorosWrapper = TensorWrapper<VolumesTensorN, porosIndex, dim>;
        const PorosWrapper porosityField(volumeTensorField);

        const GetNeighInfo::GetNeighInfo<T, PorosWrapper, FlagsWrapper, dim> getter(porosityField, m_thresh, flagsField);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const Dot& pos = *it;

            CellNeighborInfo result = getter.get(pos);
            handleNeighVolumes<T, dim, porosIndex>(pos, volumeTensorField, result);

            const Dot posArr = pos + location - m_location;
            DataAccess::get(m_arr, posArr) = result;
        }
    }

private:
    Array& m_arr;
    const Dot& m_location;
    const T m_thresh;
};

template<typename GeomFlag, typename T, typename Array, size_t dim, size_t solidFlags, size_t porosIndex>
class CalcNeighInfosAndGeomFlags : public
        PlbTypeDeduction::GetBoxProcessingFunctionalXD_SN<GeomFlag, T, dim>::value
{
public:
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

    CalcNeighInfosAndGeomFlags(Array& arr, const Dot& location, const T& thresh)
        : m_arr(arr)
        , m_location(location)
        , m_thresh(thresh)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulk;
    }

    virtual CalcNeighInfosAndGeomFlags<GeomFlag, T, Array, dim, solidFlags, porosIndex>* clone() const override
    {
        return new CalcNeighInfosAndGeomFlags<GeomFlag, T, Array, dim, solidFlags, porosIndex>(*this);
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_SN<GeomFlag, T, dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using ShapeField = typename Traits::template argument<2>::type;
    using VolumesTensorN = typename Traits::template argument<3>::type;


public:
    virtual void process(Box domain, ShapeField& shapeMask, VolumesTensorN& volumeTensorField) override
    {
        const Dot location = shapeMask.getLocation();
        const Dot offset = plb::computeRelativeDisplacement(shapeMask, volumeTensorField);

        using FlagsWrapper = TensorWrapper<VolumesTensorN, solidFlags, dim>;
        const FlagsWrapper flagsField(volumeTensorField);

        using PorosWrapper = TensorWrapper<VolumesTensorN, porosIndex, dim>;
        const PorosWrapper porosityField(volumeTensorField);

        const GetNeighInfo::GetNeighInfo<T, PorosWrapper, FlagsWrapper, dim> getter(porosityField, m_thresh, flagsField);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const Dot& pos = *it;
            const Dot posIn = pos + offset;

            CellNeighborInfo result = getter.get(posIn);

            // bitmask of GeomFlag = SSSS LLLL
            const GeomFlag geomFlag = result.shape + result.location;
            DataAccess::get(shapeMask, pos) = geomFlag;

            handleNeighVolumes<T, dim, porosIndex>(posIn, volumeTensorField, result);

            const Dot posArr = pos + location - m_location;
            DataAccess::get(m_arr, posArr) = result;
        }
    }

private:
    Array& m_arr;
    const Dot& m_location;
    const T m_thresh;
};



template<typename T, typename Array, size_t dim>
class FixConvexNeighbors : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_S<T, dim>::value
{
public:
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

    FixConvexNeighbors(Array& arr, const Dot& location)
        : m_arr(arr)
        , m_location(location)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulk;
    }

    virtual FixConvexNeighbors<T,Array,dim>* clone() const override
    {
        return new FixConvexNeighbors<T,Array,dim>(*this);
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_S<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputField = typename Traits::template argument<2>::type;

public:
    virtual void process(Box domain, InputField& geomFlags) override
    {
        const Dot fieldLocation = geomFlags.getLocation();
        const Dot offset = fieldLocation - m_location;

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const Dot& pos = *it;
            const T ownGeomFlags = DataAccess::get(geomFlags, pos);

            if((ownGeomFlags & CellNeighborInfo::S_Flat)
                    != CellNeighborInfo::S_Flat)
            {
                const Dot posArr = pos + offset;
                CellNeighborInfo& ownInfo = DataAccess::get(m_arr, posArr);

                // disable precipitation for surface if at least one neighbor is concave

                //  using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;
                //  const DotList neighs = DotFindNeighbors::find(pos);
                //  for(const Dot& neighPos : neighs.dots)

                const Box localBox = PalabosConversionTools::toBox(pos);
                const Box localDomain = localBox.enlarge(1);
                const auto localEnd = DataAccess::end(localDomain);
                for(auto localIt = DataAccess::begin(localDomain); localIt < localEnd; ++localIt)
                {
                    const Dot& neighPos = *localIt;
                    if(neighPos == pos)
                    {
                        continue;
                    }

                    const T neighGeomFlags = DataAccess::get(geomFlags, neighPos);

                    if((neighGeomFlags & CellNeighborInfo::S_Concave)
                            == CellNeighborInfo::S_Concave)
                    {
                        ownInfo.noPrecip = true;
                        break;
                    }
                }
            }
            else if((ownGeomFlags & CellNeighborInfo::S_Convex)
                    != CellNeighborInfo::S_Convex)
            {
                const Dot posArr = pos + offset;
                CellNeighborInfo& ownInfo = DataAccess::get(m_arr, posArr);

                // disable precipitation for edge/corner if at least one neighbor is flat or concave

                const Box localBox = PalabosConversionTools::toBox(pos);
                const Box localDomain = localBox.enlarge(1);
                const auto localEnd = DataAccess::end(localDomain);
                for(auto localIt = DataAccess::begin(localDomain); localIt < localEnd; ++localIt)
                {
                    const Dot& neighPos = *localIt;
                    if(neighPos == pos)
                    {
                        continue;
                    }

                    const T neighGeomFlags = DataAccess::get(geomFlags, neighPos);

                    if((neighGeomFlags & (CellNeighborInfo::S_Flat | CellNeighborInfo::S_Concave))
                            == (CellNeighborInfo::S_Flat | CellNeighborInfo::S_Concave))
                    {
                        ownInfo.noPrecip = true;
                        break;
                    }
                }
            }
        }
    }

private:
    Array& m_arr;
    const Dot& m_location;
};

}

static bool isInBoundaryConditions(const std::vector<std::string>& compNames,
                                   const AbstractBoundaryConditions&bc,
                                   const IPPVector3DLong& pos)
{
    const IPPVector3DInt posInt = convert<int>(pos);

    {
        const AbstractBoundaryConditions::BCVec& advBc = bc.getAdvectiveBoundaryConditions();

        for(const AbstractBoundaryConditions::BoundaryConditionData& data : advBc)
        {
            if(data.domain.range.contains(posInt))
            {
                return true;
            }
        }
    }

    for(const std::string& name : compNames)
    {
        const AbstractBoundaryConditions::BCVec& diffBc = bc.getDiffusiveBoundaryConditions(name);
        for(const AbstractBoundaryConditions::BoundaryConditionData& data : diffBc)
        {
            if(data.domain.range.contains(posInt))
            {
                return true;
            }
        }
    }

    return false;
}



template<typename TransportTraits>
ReactivePalabosTransportModule<TransportTraits>
::ReactivePalabosTransportModule(TransportModuleConfigPtr& config)
    : BC(config)
    , m_geometrySync(nullptr)
{
    assert(this->m_decomp);
    using CalcInterfProp = PalabosCalcInterfaceProperties<Scalar, MaskType, dim>;
    this->m_nonLocalOperations.interfaceProperties.reset(new CalcInterfProp(0.5/dim,
                                                                            *this->m_decomp,
                                                                            m_oldNodes.interfaceNodes));
}

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::setGeomSync(AbstractGeometrySync *geomSync)
{
    assert(m_geometrySync == nullptr);
    m_geometrySync = geomSync;
}


template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::initModule()
{
    MPI_CHECK_SYNC;

    this->initFields();
    this->updateFieldPorosity(); // here and later once more? why? bug?

    this->initInertCellsMarker();

    const AbstractBoundaryConditions& bc = *this->m_simData.getBoundaryConditions();
    const std::vector<int>& periodic = bc.getPeriodicBoundaryDimensions();

    if(this->m_conf->flowFunc->isEnabled())
    {
        const AbstractBoundaryConditions::BCVec& advectiveBC = bc.getAdvectiveBoundaryConditions();
        this->initHydrodynamicLattice(periodic, advectiveBC);
        this->prepareResultArrays(this->m_results.flowResults);
    }

    MPI_CHECK_SYNC;

    this->initDiffusionLattices(bc);
    MPI_CHECK_SYNC;

    this->updateFieldPorosity();
    this->initLatticePorosity();

    this->updateFieldTransportScalar();
    this->updateLatticeTransportScalar();
    MPI_CHECK_SYNC;

    UpdateDistanceField::execute<dim>(this->m_simData,
                                      this->m_conf->distTransPorosThresh);
    MPI_CHECK_SYNC;

    this->initLatticeConcentrations();
    MPI_CHECK_SYNC;

    this->initBounceBackNodes();
    MPI_CHECK_SYNC;

    this->finishLatticeInit();
    MPI_CHECK_SYNC;

    this->updatePostTransportConcField();
    MPI_CHECK_SYNC;

    if(this->m_geometrySync->needNeighMinPoros())
    {
        this->updateNeighInfos();
    }

}

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::initLatticeConcentrations()
{
    const ScopedDisableFloatingPointException fpe;

    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;

    const auto concView = ConcArrayViewTool::makeConcView<dim>(this->m_simData,
                                                               this->m_simData.getPostReacConc());

    const size_t nComps = this->m_palabosData.diffLattices.size();
    for(size_t iComp = 0; iComp < nComps; ++iComp)
    {
        const auto concSlice = concView[iComp];
        ScalarField& conc = this->m_palabosData.postReactionConcentrations[iComp];

        using Sync = SetScalarsFromArray<const double, dim>;
        plb::setToFunction(conc, conc.getBoundingBox(), Sync(concSlice, location));

        DiffusionLatticePtr diffLattice = this->m_palabosData.diffLattices[iComp];


        using InitFunc = typename TransportTraits::PorosityFunc;
        InitFunc::apply(*diffLattice, conc);
    }
}


template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::addPostStreamOperation(DiffusionLattice& lattice)
{
    using PushOperation = PushBounceBackPopulations<Scalar, DIFFUSION_DESCRIPTOR, MaskType>;
    plb::integrateProcessingFunctional(new PushOperation, lattice.getBoundingBox(),
                                       lattice, *this->m_palabosData.permMask);
}



template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::initFields()
{
    m_oldPorosity.reset(this->m_plbFactory->template createMultiScalarFieldPtr<Scalar>(-1.0));
    m_unityField.reset(this->m_plbFactory->template createMultiScalarFieldPtr<Scalar>(1.0));

    //this->m_nonLocalOperations.distExcess->init();

    using CalcInterfProp = PalabosCalcInterfaceProperties<Scalar, MaskType, dim>;
    CalcInterfProp& interfProp = static_cast<CalcInterfProp&>(*this->m_nonLocalOperations.interfaceProperties);
    interfProp.initFields(&(*this->m_palabosData.permMask), *this->m_plbFactory);

    const size_t nCells = this->m_decomp->getLocalNumberOfCells();
    const size_t nComps = this->m_simData.getCompNames().size();
    m_negConcBuffer.resize(nCells * nComps, 0.0);


    const size_t nNuclePhases = this->m_simData.nucleationPhasesVolumeFractions.size();
    const size_t nFields = nNuclePhases + 2;
    // init field with default values
    std::vector<Scalar> arr(nFields, 0.0);
    arr[0] = IPPConstants::s_isNotSteadyFlag;
    arr[1] = 1.0;
    m_capilPorosVolumes.reset(this->m_plbFactory->template createMultiNTensorFieldPtr<Scalar>(nFields, arr.data()));

}

template<typename T, typename Field>
static bool hasNoNegSolute(Field& conc, size_t icomp)
{
    if(icomp != 3)
    {
        return hasNoNeg<T>(conc);
    }
    else
    {
        return true;
    }
}

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::updatePostTransportConcField()
{
    const ScopedDisableFloatingPointException fpe;

    const size_t nComps = this->m_palabosData.postTransportConcentrations.size();
    assert(nComps == this->m_palabosData.diffLattices.size());

    for(size_t iComp = 0; iComp < nComps; ++iComp)
    {
        ScalarField& postTransConc = this->m_palabosData.postTransportConcentrations[iComp];

        const DiffusionLatticePtr& diffLattice = this->m_palabosData.diffLattices[iComp];
        plb::computeDensity(*diffLattice, postTransConc, diffLattice->getBoundingBox());
        assert(hasNoNan<Scalar>(postTransConc));
    }
}

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::collectData()
{
    this->updatePostTransportConcField();

    {
        const FieldDecomposition* decomp = this->m_simData.getDecomp();
        const IPPBox3DLong& ownDomain = decomp->getOwnDomain();
        Dot location;
        PalabosConversionTools::toPlb(ownDomain.lower, location);
        const size_t nCellsLocal = decomp->getLocalNumberOfCells();
        const IPPVector3DLong& localSize = ownDomain.getSize();

        const size_t nComps = this->m_simData.getCompNames().size();

        std::vector<double>& transConc = this->m_simData.getPostTransportConc();
        assert(transConc.size() == nComps * nCellsLocal);
        (void)nCellsLocal; // disable unused variable warning in release build

        auto transConcView = ConcArrayViewTool::makeConcView<dim>(nComps, localSize, transConc);

        for(size_t iComp = 0; iComp < nComps; ++iComp)
        {
            // postTransConc is updated in updatePostTransportConcField() already
            ScalarField& postTransConc = this->m_palabosData.postTransportConcentrations[iComp];

            using ArrayViewDouble = av::array_view<double, dim>;
            using Sync = GetScalarsToArray<Scalar, ArrayViewDouble, dim>;

            ArrayViewDouble transConcSlice = transConcView[iComp];
            plb::applyProcessingFunctional(new Sync(transConcSlice, location),
                                           postTransConc.getBoundingBox(), postTransConc);

        }


        for(size_t i = 0; i < transConc.size(); ++i)
        {
            // exclude charge
            if(i / nCellsLocal == 3)
            {
                continue;
            }

            double& conc = transConc[i];
            if(conc < 0.0)
            {
                m_negConcBuffer[i] += conc;
                conc = 0.0;
            }
        }
    }
}

template<typename T>
static T calcUpperThresh(const T& lower)
{
    const T upper = (1.0 + GeometryUpdate::s_relTreshDiff) * lower;
    return fmin(1.0, upper);
}

struct CollectDiff
{
    template<typename DotList>
    static void collect(const DotList& oldNonPermNodes,
                        const DotList& /*permNodes*/, const DotList& nonPermNodes,
                        DotList& newPermNodes, DotList& newNonPermNodes)
    {
        // there are only two possibilities: permeable or non permeable for
        // each node. switch does occur from one set into the other only.
        // new perm nodes are the negation of new non perm nodes

        // solid -> liquid
        assert(newPermNodes.dots.empty());
        boost::range::set_difference(oldNonPermNodes.dots, nonPermNodes.dots,
                                     std::back_inserter(newPermNodes.dots));

        // liquid -> solid
        assert(newNonPermNodes.dots.empty());
        boost::range::set_difference(nonPermNodes.dots, oldNonPermNodes.dots,
                                     std::back_inserter(newNonPermNodes.dots));
    }

    template<typename T, template<typename U> class Descriptor, typename MaskType,
             typename Lattice, typename MaskField>
    static void apply(Lattice& lattice, MaskField& mask)
    {
        using DistOperation = DistributeBounceBackPopulations<T, Descriptor, MaskType>;
        plb::applyProcessingFunctional(new DistOperation, lattice.getBoundingBox(),
                                       lattice, mask);
    }

    template<typename MaskType, typename T, typename Field, typename MaskField, typename DotList>
    static void collectNodes(const T& porosThreshLower, const Field& porosity, const MaskField& mask,
                             DotList& permNodes, DotList& nonPermNodes)
    {
        static constexpr size_t dim = PlbTypeDeduction::GetFieldDim<T, Field>::value;

        const T porosThreshUpper = calcUpperThresh(porosThreshLower);

        using CollectSolidCoords = CollectGeomChangedCoordinates<T, MaskType, dim>;
        plb::applyProcessingFunctional(new CollectSolidCoords(porosThreshLower, porosThreshUpper,
                                                              permNodes, nonPermNodes),
                                       porosity.getBoundingBox(),
                                       const_cast<Field&>(porosity),
                                       const_cast<MaskField&>(mask));
    }
};

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::syncDataIntoLattices()
{
    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong size = ownDomain.getSize();

    ScalarField& physicalPorosity = this->m_palabosData.getPorosity();

    const size_t nComps = this->m_simData.getCompNames().size();
    std::vector<double>& concRaw = this->m_simData.getPostReacConc();

    // unloading negative buffer
    const size_t nCells = ownDomain.getDiagonalSize();

    for(size_t i = 0; i < concRaw.size(); ++i)
    {
        // exclude charge
        if(i / nCells == 3)
        {
            continue;
        }

        concRaw[i] += m_negConcBuffer[i];
        m_negConcBuffer[i] = 0.0;
    }

    const auto concView = ConcArrayViewTool::makeConcView<dim>(nComps, size, concRaw);

    const auto bounds = av::MakeBounds<dim>::get(size[0], size[1], size[2]);
    const av::array_view<const Scalar, dim> transScalarView(this->m_transportScalar, bounds);

    Dot plbLocation;
    PalabosConversionTools::toPlb(ownDomain.lower, plbLocation);


    for(size_t iComp = 0; iComp < nComps; ++iComp)
    {
        DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];
        const auto concSlice = concView[iComp];
        ScalarField& porosity = iComp == 0 ? *m_unityField : physicalPorosity;

        using SyncData = SyncDataIntoLattice<Scalar, double, DIFFUSION_DESCRIPTOR>;
        plb::applyProcessingFunctional(new SyncData(plbLocation, transScalarView, concSlice),
                                       this->m_transportSyncEnabledBox, lattice, porosity);

    }

    // TODO: update transport scalar in hydrodynamic lattice
    if(this->m_palabosData.flowLattice)
    {
        throw std::runtime_error("implemetation for hydrodynamic lattice missing");
    }
}


template<typename T, int scalarBeginsAt, typename Box, typename Dot,
         typename ScalarFieldMask, typename Lattice>
static void handleInterfaces(const T& weight, const Box& domain, const Dot& location,
                             const ScalarFieldMask& permMask, Lattice& lattice)
{
    const Dot latticeLocation = lattice.getLocation();
    const Dot offsetPermMask = plb::computeRelativeDisplacement(lattice, permMask);

    static constexpr size_t dim = PlbTypeDeduction::GetDotDim<Dot>::value;

    const auto end = DataAccess::end(domain);

    // fill some scratch space with old concentrations

    IPPBox3DInt ippBox;
    PalabosConversionTools::convertBox(domain, ippBox);
    ArrayDimensionConvert indexConv(ippBox.getSize());

    std::vector<T> concField(indexConv.getNxyz());
    for(auto it = DataAccess::begin(domain); it < end; ++it)
    {

        const Dot& pos = *it;
        const auto& cell = DataAccess::get(lattice, pos);
        const T& cCurr = *cell.getExternal(scalarBeginsAt);

        const Dot absPos = pos + latticeLocation - location;
        IPPVector3DUInt ippPos;
        PalabosConversionTools::convertToIPP(absPos, ippPos);
        const size_t index = indexConv.calcIndex(ippPos);

        assert(concField.size() > index);
        concField[index] = cCurr;
    }

    // pulling concentrations from interface nodes
    for(auto it = DataAccess::begin(domain); it < end; ++it)
    {
        const Dot& pos = *it;

        const int mask = DataAccess::get(permMask, pos);
        if(mask & PF_isPermeable)
        {
            auto& cell = DataAccess::get(lattice, pos);
            T& cCurr = *cell.getExternal(scalarBeginsAt);

            const Dot absPos = pos + latticeLocation - location;
            IPPVector3DUInt ippPos;
            PalabosConversionTools::convertToIPP(absPos, ippPos);
            const size_t index = indexConv.calcIndex(ippPos);

            assert(concField.size() > index);
            const T& cOrig = concField[index];


            using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;
            const DotList neighList = DotFindNeighbors::find(pos);
            assert(neighList.getN() == dim * 2);

            for(const Dot& neigh : neighList.dots)
            {
                const int maskNeigh = DataAccess::get(permMask, neigh + offsetPermMask);

                if( ((maskNeigh & PF_isPermeable) == false) && (maskNeigh & PF_isInterface) )
                {
                    const auto& cellNeigh = DataAccess::get(lattice, neigh);
                    const T& cInterface = *cellNeigh.getExternal(scalarBeginsAt);

                    const T cDiff = cInterface - cOrig;
                    const T cDiffSource = weight * cDiff;

                    cCurr += cDiffSource;
                }
            }

        }

    }


    // setting interface node conc to zero
    for(auto it = DataAccess::begin(domain); it < end; ++it)
    {
        const Dot pos = *it;
        const int mask = DataAccess::get(permMask, pos);

        if( ((mask & PF_isPermeable) == false) && (mask & PF_isInterface) )
        {
            auto& cell = DataAccess::get(lattice, pos);
            T& cCurr = *cell.getExternal(scalarBeginsAt);
            cCurr = 0;
        }

    }

}

template<typename T, int scalarBeginsAt, int porosityBeginsAt>
struct UpdateSourceTermRespectingPoros
{
    template<typename Box, typename ScalarField, typename Dot, typename ArrayView, typename Lattice>
    static void calc(const Box& domain, const ScalarField& oldPorosity,
                     const ScalarField& postTransConc, const Dot& location,
                     ArrayView& postReacConc, Lattice& lattice)
    {
        const auto latticeLocation = lattice.getLocation();
        const Dot offsetPoros = plb::computeRelativeDisplacement(lattice, oldPorosity);
        const Dot offsetPostTransConc = plb::computeRelativeDisplacement(lattice, postTransConc);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const Dot& pos = *it;

            auto& cell = DataAccess::get(lattice, pos);

            const T& oldConc = DataAccess::get(postTransConc, pos + offsetPostTransConc);
            const T& oldPoros = DataAccess::get(oldPorosity, pos + offsetPoros);
            const T& newPoros = *cell.getExternal(porosityBeginsAt);

            T& source = *cell.getExternal(scalarBeginsAt);

            const Dot& absPos = pos + latticeLocation - location;
            DataAccess::get(postReacConc, absPos) = static_cast<typename ArrayView::value_type>(source);

            source = source * newPoros - oldConc * oldPoros;
        }
    }

};

template<typename T, int scalarBeginsAt>
struct UpdateSourceTermIgnorePoros
{
    template<typename Box, typename ScalarField, typename Dot,
             typename ArrayView, typename Lattice>
    static void calc(const Box& domain, const ScalarField& /*oldPorosity*/,
                     const ScalarField& postTransConc, const Dot& location,
                     ArrayView& postReacConc, Lattice& lattice)
    {
        const Dot latticeLocation = lattice.getLocation();
        const Dot offsetPostTransConc = plb::computeRelativeDisplacement(lattice, postTransConc);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const Dot& pos = *it;

            auto& cell = DataAccess::get(lattice, pos);

            const T& oldConc = DataAccess::get(postTransConc, pos + offsetPostTransConc);
            T& source = *cell.getExternal(scalarBeginsAt);

            const Dot& absPos = pos + latticeLocation - location;
            DataAccess::get(postReacConc, absPos) = static_cast<typename ArrayView::value_type>(source);

            source = source - oldConc;
        }
    }


};

template<typename T, int scalarBeginsAt, typename Box, typename Dot, typename Lattice>
static void handleTransportScalars(const Box& domain, const Dot& location,
                                   Lattice& lattice)
{
    IPPBox3DInt ippBox;
    PalabosConversionTools::convertBox(domain, ippBox);
    ArrayDimensionConvert indexConv(ippBox.getSize());

    const Dot latticeLocation = lattice.getLocation();

    std::vector<T> oldTransScalar(indexConv.getNxyz());

    const auto end = DataAccess::end(domain);
    for(auto it = DataAccess::begin(domain); it < end; ++it)
    {
        const Dot& pos = *it;
        const auto& cell = DataAccess::get(lattice, pos);
        const T& transScalar = *cell.getExternal(scalarBeginsAt);

        const Dot absPos = pos + latticeLocation - location;
        IPPVector3DUInt ippPos;
        PalabosConversionTools::convertToIPP(absPos, ippPos);
        const size_t index = indexConv.calcIndex(ippPos);

        assert(oldTransScalar.size() > index);
        oldTransScalar[index] = transScalar;
    }


    for(auto it = DataAccess::begin(domain); it < end; ++it)
    {
        const Dot& pos = *it;

        const Dot absPos = pos + latticeLocation - location;
        IPPVector3DUInt ippPos;
        PalabosConversionTools::convertToIPP(absPos, ippPos);
        const size_t posIndex = indexConv.calcIndex(ippPos);

        assert(oldTransScalar.size() > posIndex);
        T transScalar = oldTransScalar[posIndex];

        const auto neighList = DotFindNeighbors::find(pos);

        // T harmMeanNeigh = 0;
        // size_t nValidNeigh = 0;

        for(const Dot& neigh : neighList.dots)
        {
            if(plb::contained(neigh, domain))
            {
                // ++nValidNeigh;

                const Dot absPosNeigh = neigh + latticeLocation - location;
                IPPVector3DUInt ippPosNeigh;
                PalabosConversionTools::convertToIPP(absPosNeigh, ippPosNeigh);
                const size_t neighIndex = indexConv.calcIndex(ippPosNeigh);

                assert(oldTransScalar.size() > neighIndex);
                const T& transScalarNeigh = oldTransScalar[neighIndex];
                assert(transScalarNeigh > 0.0);

                transScalar = std::min(transScalar, transScalarNeigh);
                // harmMeanNeigh += 1.0 / transScalarNeigh;
            }
        }

        // assert(harmMeanNeigh > 0);
        // harmMeanNeigh = nValidNeigh / harmMeanNeigh;
        // transScalar = std::min(transScalar, harmMeanNeigh);

        auto& cell = DataAccess::get(lattice, pos);
        *cell.getExternal(scalarBeginsAt) = transScalar;
    }
}

template<typename T,
         typename Array,
         template<typename U> class Descriptor,
         size_t dim,
         typename SourceFunc>
class UpdateLattice : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD<dim>::value
{
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

public:
    UpdateLattice(const T& weight, const Dot& location,
                  Array& postReacConc)
        : m_weight(weight)
        , m_location(location)
        , m_postReacConc(postReacConc)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        assert(modified.size() == 4);

        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
        modified[2] = plb::modif::nothing;
        modified[3] = plb::modif::nothing;
    }

    virtual UpdateLattice<T, Array, Descriptor, dim, SourceFunc>* clone() const override
    {
        return new UpdateLattice<T, Array, Descriptor, dim, SourceFunc>(*this);
    }

private:
    using Lattice = typename PlbTypeDeduction::GetBlockLattice_XD<T, Descriptor>::value;
    using ScalarField = typename PlbTypeDeduction::GetScalarField_XD<T, dim>::value;
    using ScalarFieldInt = typename PlbTypeDeduction::GetScalarField_XD<int, dim>::value;

    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD<dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::processGenericBlocks)>;
    using Box = typename Traits::template argument<1>::type;
    using AtomicBlockVec = typename Traits::template argument<2>::type;


public:
    virtual void processGenericBlocks(Box domain, AtomicBlockVec fields) override
    {
        assert(fields.size() == 4);

        Lattice& lattice = *dynamic_cast<Lattice*>(fields[0]);
        ScalarFieldInt& permMask = *dynamic_cast<ScalarFieldInt*>(fields[3]);
        static constexpr int scalarIndex = Descriptor<T>::ExternalField::scalarBeginsAt;
        handleInterfaces<T, scalarIndex>(m_weight, domain, m_location, permMask, lattice);

        ScalarField& oldPorosity = *dynamic_cast<ScalarField*>(fields[1]);
        ScalarField& postTransConc = *dynamic_cast<ScalarField*>(fields[2]);
        SourceFunc::template calc(domain, oldPorosity, postTransConc,
                                  m_location, m_postReacConc, lattice);


        static constexpr int transportScalarIndex =
                Descriptor<T>::ExternalField::transportScalarBeginsAt;
        handleTransportScalars<T, transportScalarIndex>(domain, m_location, lattice);

    }



private:
    const T m_weight;
    const Dot m_location;
    Array& m_postReacConc;

};


template<typename T, template<typename U> class Descriptor, typename MaskType>
class UpdateConc : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_LS<T, Descriptor, MaskType>::value
{
    static constexpr size_t dim = Descriptor<T>::d;
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

public:
    UpdateConc(const T& weight, const Dot& location)
        : m_weight(weight)
        , m_location(location)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        assert(modified.size() == 2);

        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
    }

    virtual UpdateConc<T, Descriptor, MaskType>* clone() const override
    {
        return new UpdateConc<T, Descriptor, MaskType>(*this);
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_LS<T, Descriptor, MaskType>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;
    using ScalarFieldMask = typename Traits::template argument<3>::type;


public:
    virtual void process(Box domain, Lattice& lattice, ScalarFieldMask& permMask) override
    {
        static constexpr int scalarIndex = Descriptor<T>::ExternalField::scalarBeginsAt;

        handleInterfaces<T, scalarIndex>(m_weight, domain, m_location, permMask, lattice);
    }

private:
    const T m_weight;
    const Dot m_location;

};


template<typename T,
         typename U,
         template<typename V> class Descriptor,
         typename MaskType, typename SourceFunc>
class UpdateSourceAndTransportScalar : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD<Descriptor<T>::d>::value
{
private:
    static constexpr size_t dim = Descriptor<T>::d;
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

public:
    UpdateSourceAndTransportScalar(const Dot& location, av::array_view<U, dim>& postReacConc)
        : m_location(location)
        , m_postReacConc(postReacConc)
    {

    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        assert(modified.size() == 4);

        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
        modified[2] = plb::modif::nothing;
        modified[3] = plb::modif::nothing;
    }

    virtual UpdateSourceAndTransportScalar<T, U, Descriptor, MaskType, SourceFunc>* clone() const override
    {
        return new UpdateSourceAndTransportScalar<T, U, Descriptor, MaskType, SourceFunc>(*this);
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD<dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::processGenericBlocks)>;
    using Box = typename Traits::template argument<1>::type;
    using AtomicBlockVec = typename Traits::template argument<2>::type;

public:
    virtual void processGenericBlocks(Box domain, AtomicBlockVec fields) override
    {
        assert(fields.size() == 4);

        using Lattice = typename PlbTypeDeduction::GetBlockLattice_XD<T, Descriptor>::value;
        using ScalarField = typename PlbTypeDeduction::GetScalarField_XD<T, dim>::value;

        Lattice &lattice = *dynamic_cast<Lattice *>(fields[0]);
        ScalarField& oldPorosity = *dynamic_cast<ScalarField*>(fields[1]);
        ScalarField& postTransConc = *dynamic_cast<ScalarField*>(fields[2]);

        SourceFunc::calc(domain, oldPorosity, postTransConc,
                         m_location, m_postReacConc, lattice);

        static constexpr int transportScalarIndex =
                Descriptor<T>::ExternalField::transportScalarBeginsAt;
        handleTransportScalars<T, transportScalarIndex>(domain, m_location, lattice);
    }

private:
    const Dot m_location;
    av::array_view<U, dim>& m_postReacConc;
};


class GeometrySync
{
public:
    GeometrySync(AbstractGeometrySync& geomSync,
                 std::vector<double>& conc)
        : m_geomSync(geomSync)
        , m_conc(conc)
    {

    }

    template<typename ...T>
    void run(T... args)
    {
        m_geomSync.run(args..., m_conc);
    }

private:
    AbstractGeometrySync& m_geomSync;
    std::vector<double>& m_conc;
};

template<typename DotList>
static void syncGeomData(const DotList& newPermNodes,
                         const DotList& newInterfaceNodes,
                         const DotList& newNonPermNonInterfaceNodes,
                         GeometrySync& geometrySync)
{
    std::vector<IPPVector3DLong> newPerm;
    PalabosConversionTools::convertTo3D(newPermNodes, newPerm);

    std::vector<IPPVector3DLong> newInterface;
    PalabosConversionTools::convertTo3D(newInterfaceNodes, newInterface);

    std::vector<IPPVector3DLong> newNonPermNonInterface;
    PalabosConversionTools::convertTo3D(newNonPermNonInterfaceNodes, newNonPermNonInterface);

    geometrySync.run(newPerm, newInterface, newNonPermNonInterface);
}

template<typename Box, typename DotList>
static void filter(const Box& box,
                   const DotList& input,
                   DotList& output)
{
    for(const auto& pos : input.dots)
    {
        if(plb::contained(pos, box))
        {
            output.dots.push_back(pos);
        }
    }
}

template<typename DotList, typename Box>
static void syncGeomData(const DotList& newPermNodes,
                         const DotList& newInterfaceNodes,
                         const DotList& newNonPermNonInterfaceNodes,
                         Box& ownBox,
                         GeometrySync& geometrySync)
{
    DotList newPermNodesFilt;
    filter(ownBox, newPermNodes, newPermNodesFilt);

    DotList newInterfaceNodesFilt;
    filter(ownBox, newInterfaceNodes, newInterfaceNodesFilt);

    DotList newNonPermNonInterfaceNodesFilt;
    filter(ownBox, newNonPermNonInterfaceNodes, newNonPermNonInterfaceNodesFilt);

    syncGeomData(newPermNodesFilt, newInterfaceNodesFilt,
                 newNonPermNonInterfaceNodesFilt, geometrySync);
}

template<typename DotList>
struct CopyAllNodes
{
    CopyAllNodes(const DotList& boundaryNodes)
        : m_boundaryNodes(boundaryNodes)
    {

    }

    void collect(const DotList & /*oldNonPermNodes*/,
                 const DotList & /*permNodes*/,
                 const DotList &nonPermNodes,
                 DotList & /*newPermNodes*/,
                 DotList &newNonPermNodes) const
    {
        // constraints:
        // 1) boundary nodes must never be newPermNodes! otherwise the
        //      boundary dynamics getting overwritten
        // 2) do not claim cells as newly permeable on init, because this
        //      leads the reaction module to react as it has to init this cells with water or something
        //  boost::range::set_difference(permNodes.dots, m_boundaryNodes.dots,
        //                               std::back_inserter(newPermNodes.dots));

        newNonPermNodes = nonPermNodes;
    }

    template<typename T, template<typename U> class Descriptor, typename MaskType,
             typename Lattice, typename MaskField>
    static void apply(Lattice& /*lattice*/, MaskField& /*mask*/)
    {

    }

    template<typename MaskType, typename T, typename Field, typename MaskField>
    static void collectNodes(const T& porosThresh, const Field& porosity, const MaskField& /*mask*/,
                             DotList& permNodes, DotList& nonPermNodes)
    {
        static constexpr size_t dim = PlbTypeDeduction::GetFieldDim<T, Field>::value;

        using CollectSolidCoords = CollectThreshCoordinates<T, dim>;
        plb::applyProcessingFunctional(new CollectSolidCoords(porosThresh, permNodes, nonPermNodes),
                                       porosity.getBoundingBox(), const_cast<Field&>(porosity));
    }

private:
    const DotList& m_boundaryNodes;
};

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::addAdditionalSource()
{
    MPI_CHECK_SYNC;

    const std::vector<CellTotalsDiff>& diffPerCell = this->m_simData.additionalSources;

    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& origin = ownDomain.lower;
    const ArrayDimensionConvert indexConv(ownDomain.getSize());

    // collect all sources
    using SetScalars = SetDotToNSource<Scalar, dim>;
    typename SetScalars::CellToDiffMap cellToDiff;

    DotList allDots;

    for(const CellTotalsDiff& source : diffPerCell)
    {
        const size_t iCell = source.iCell;
        const IPPVector3DInt localCoord = indexConv.calcCoordinate(iCell);
        const IPPVector3DLong globalCoord = origin + localCoord;

        Dot plbCoordGlobal;
        PalabosConversionTools::toPlb(globalCoord, plbCoordGlobal);

        allDots.addDot(plbCoordGlobal);
        cellToDiff[plbCoordGlobal] = &source;
    }

    const size_t nComps = this->m_palabosData.diffLattices.size();

    // TODO: maybe construct only once
    auto tensorSource = this->m_plbFactory->template createMultiNTensorField<Scalar>(nComps);

    plb::applyProcessingFunctional(new SetScalars(cellToDiff),
                                   allDots, tensorSource);


    using MultiBlock = typename PlbTypeDeduction::GetMultiBlock_XD<dim>::value;
    std::vector<MultiBlock*> fields;
    fields.push_back(&(*this->m_palabosData.permMask));
    fields.push_back(&tensorSource);
    for(size_t iComp = 0; iComp < nComps; ++iComp)
    {
        DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];
        fields.push_back(&lattice);
    }

    av::array_view<double, dim+1> concView =
            ConcArrayViewTool::makeConcView<dim>(this->m_simData,
                                                 this->m_simData.getPostReacConc());

    Dot plbLocation;
    PalabosConversionTools::toPlb(ownDomain.lower, plbLocation);

    using SetSource = SetNTensorToSource<Scalar, double, DIFFUSION_DESCRIPTOR, MaskType>;
    plb::applyProcessingFunctional(new SetSource(plbLocation, concView),
                                   this->m_transportSyncEnabledBox, fields);

    MPI_CHECK_SYNC;
}


template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::updateNeighInfos()
{    
    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    const IPPVector3DLong ownSize = ownDomain.getSize();
    const av::bounds<dim> bounds = av::MakeBounds<dim>::get(ownSize[0], ownSize[1], ownSize[2]);


    {
        using ArrayView = av::array_view<const double, dim>;
        using ArrayViewVec = std::vector<ArrayView>;
        using Setter = SetNTensorFromArray<Scalar, ArrayViewVec, dim>;

        ArrayViewVec arrViewVec;

        const std::vector<double>& solidFlags = this->m_simData.solidFlags;
        arrViewVec.emplace_back(ArrayView(solidFlags, bounds));

        const std::vector<double>& capilPoros = this->m_simData.getCapillaryPorosity();
        arrViewVec.emplace_back(ArrayView(capilPoros, bounds));


        const std::vector<std::vector<double>>& nuclPhases = this->m_simData.nucleationPhasesVolumeFractions;
        for(size_t i = 0; i < nuclPhases.size(); ++i)
        {
            const std::vector<double>& nuclPhase = nuclPhases[i];
            arrViewVec.emplace_back(ArrayView(nuclPhase, bounds));
        }

        using Sync = SetToNTensorFunction<Scalar, Setter, dim>;
        plb::applyProcessingFunctional(new Sync(Setter(arrViewVec, location, arrViewVec.size())),
                                       this->m_transportSyncEnabledBox, *m_capilPorosVolumes);
    }


    const Dot dotOrigin = PalabosConversionTools::ToPlb<dim>::convertVector(location);

    using CellNeighArrView = av::array_view<CellNeighborInfo, dim>;
    AbstractGeometrySync::NeighborInfos& neighInfos = this->m_geometrySync->getNeighInfos();
    std::vector<CellNeighborInfo>& neighInfosResult = neighInfos.results;
    CellNeighArrView neighInfosView(neighInfosResult, bounds);

    using CalcInfos = CalcNeighInfos<Scalar, CellNeighArrView, dim, 0, 1>;
    plb::applyProcessingFunctional(new CalcInfos(neighInfosView, dotOrigin,
                                                 neighInfos.porosityThresh),
                                   this->m_transportSyncEnabledBox,
                                   *m_capilPorosVolumes);
}


template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::updatePostReactionState()
{
    MPI_CHECK_SYNC;

    MPI_Request request;
    bool needAdditionalSource = false;
    {
        bool needAdditionalSourceLocal = this->m_simData.additionalSources.empty() == false;
        MPI_Iallreduce(&needAdditionalSourceLocal, &needAdditionalSource, 1,
                       MPI_CXX_BOOL, MPI_LOR, MPIManager::getInstance().getCommunicator(),
                       &request);

    }

    // update intermediate fields
    BENCH(this->updateFieldPorosity());
    BENCH(this->updateFieldTransportScalar());

    // push all data into the lattices
    BENCH(this->syncDataIntoLattices());

    MPI_CHECK_SYNC;

    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    Dot plbLocation;
    PalabosConversionTools::toPlb(location, plbLocation);

    const bool hasGeomChange = hasGeomChanged();
    if(hasGeomChange == false)
    {
        MPI_CHECK_SYNC;

        using MultiBlock = typename PlbTypeDeduction::GetMultiBlock_XD<dim>::value;

        std::vector<MultiBlock*> fields(4, nullptr);
        fields[1] = &(*m_oldPorosity);
        fields[3] = &(*this->m_palabosData.permMask);

        std::vector<double>& concRaw = this->m_simData.getPostReacConc();
        const auto concView = ConcArrayViewTool::makeConcView<dim>(this->m_simData, concRaw);

        for(size_t iComp = 0; iComp < this->m_palabosData.diffLattices.size(); ++iComp)
        {
            DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];
            fields[0] = &lattice;
            fields[2] = &(this->m_palabosData.postTransportConcentrations[iComp]);

            auto concSlice = concView[iComp];
            using Array = decltype(concSlice);

            using DiffusionDesc = typename DIFFUSION_DESCRIPTOR<Scalar>;
            static constexpr int scalarIndex = DiffusionDesc::ExternalField::scalarBeginsAt;


            if(iComp == 0)
            {
                using SourceFunc = UpdateSourceTermIgnorePoros<Scalar, scalarIndex>;
                using Update = UpdateLattice<Scalar, Array, DIFFUSION_DESCRIPTOR, dim, SourceFunc>;
                plb::applyProcessingFunctional(new Update(0.5/dim, plbLocation, concSlice),
                                               this->m_transportSyncEnabledBox, fields);
            }
            else
            {
                static constexpr int porosIndex = DiffusionDesc::ExternalField::porosityBeginsAt;
                using SourceFunc = UpdateSourceTermRespectingPoros<Scalar, scalarIndex, porosIndex>;
                using Update = UpdateLattice<Scalar, Array, DIFFUSION_DESCRIPTOR, dim, SourceFunc>;
                plb::applyProcessingFunctional(new Update(0.5/dim, plbLocation, concSlice),
                                               this->m_transportSyncEnabledBox, fields);

            }

            MPI_CHECK_SYNC;
        }
    }
    else
    {
        // geometry has changed
        MPI_CHECK_SYNC;

        if(this->m_geometrySync->needDistField())
        {
            UpdateDistanceField::execute<dim>(this->m_simData,
                                              this->m_conf->distTransPorosThresh);
        }

        // update geometry if needed. this must be performed after post reac conc
        // is updated since distribution at new non perm nodes which
        // were perm nodes before must be performig the pushing operation

        const av::array_view<double, dim+1> concView =
                ConcArrayViewTool::makeConcView<dim>(this->m_simData, this->m_simData.getPostReacConc());

        for(size_t iComp = 0; iComp < this->m_palabosData.diffLattices.size(); ++iComp)
        {
            DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];

            using UpdateCon = UpdateConc<Scalar, DIFFUSION_DESCRIPTOR, MaskType>;
            plb::applyProcessingFunctional(new UpdateCon(0.5/dim, plbLocation),
                                           this->m_transportSyncEnabledBox,
                                           lattice, *this->m_palabosData.permMask);
        }


        this->updateGeometry(CollectDiff());

        using MultiBlock = typename PlbTypeDeduction::GetMultiBlock_XD<dim>::value;
        std::vector<MultiBlock*> fields(4, nullptr);
        fields[1] = &(*m_oldPorosity);
        fields[3] = &(*this->m_palabosData.permMask);

        for(size_t iComp = 0; iComp < this->m_palabosData.diffLattices.size(); ++iComp)
        {
            DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];

            fields[0] = &lattice;
            fields[2] = &(this->m_palabosData.postTransportConcentrations[iComp]);


            using DiffusionDesc = typename DIFFUSION_DESCRIPTOR<Scalar>;
            static constexpr int scalarIndex = DiffusionDesc::ExternalField::scalarBeginsAt;

            if(iComp > 0)
            {
                static constexpr int porosIndex = DiffusionDesc::ExternalField::porosityBeginsAt;

                using SourceFunc = UpdateSourceTermRespectingPoros<Scalar, scalarIndex, porosIndex>;
                using UpdateS = UpdateSourceAndTransportScalar<Scalar, double, DIFFUSION_DESCRIPTOR, MaskType, SourceFunc>;

                auto concSlice = concView[iComp];
                plb::applyProcessingFunctional(new UpdateS(plbLocation, concSlice),
                                               this->m_transportSyncEnabledBox, fields);
            }
            else
            {
                using SourceFunc = UpdateSourceTermIgnorePoros<Scalar, scalarIndex>;
                using UpdateS = UpdateSourceAndTransportScalar<Scalar, double, DIFFUSION_DESCRIPTOR, MaskType, SourceFunc>;

                auto concSlice = concView[iComp];
                plb::applyProcessingFunctional(new UpdateS(plbLocation, concSlice),
                                               this->m_transportSyncEnabledBox, fields);
            }

            MPI_CHECK_SYNC;

        }
    }

    MPI_CHECK_SYNC;


    if(this->m_geometrySync->needNeighMinPoros())
    {
        // const auto start = std::chrono::steady_clock::now();

        this->updateNeighInfos();

        // const auto duration = CalcDuration::calc(start);
        // plb::pcout << "\t-> updating neighbor infos in TM took:\t"
        //            << duration.count() << " ms" << std::endl;
    }


    MPI_CHECK_SYNC;

    // this check communication and check should be almost without costs
    // communication has finished very likely in the meantime
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    if(needAdditionalSource)
    {
        MPI_CHECK_SYNC;
        this->addAdditionalSource();
    }

    MPI_CHECK_SYNC;
}


template<typename TransportTraits>
bool ReactivePalabosTransportModule<TransportTraits>::hasGeomChanged()
{
    const double porosThresholdLower = this->m_conf->porosLow;
    const double porosThresholdUpper = calcUpperThresh(porosThresholdLower);

    ScalarField& poros = this->m_palabosData.getPorosity();

    MPI_CHECK_SYNC;


    using CountGeomChanged = CountGeomChangedT<Scalar, MaskType, dim>;
    CountGeomChanged geomChange(porosThresholdLower, porosThresholdUpper);
    plb::applyProcessingFunctional(geomChange, this->m_transportSyncEnabledBox,
                                   poros, *this->m_palabosData.permMask);

    const plb::plint changedCount = geomChange.getCount();
    if(changedCount > 0)
    {
        plb::pcout << "geometry change: " << geomChange.getCount() << std::endl;
        MPI_CHECK_SYNC;
        return true;
    }
    else
    {
        MPI_CHECK_SYNC;
        return false;
    }
}


template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::initBounceBackNodes()
{
    MPI_CHECK_SYNC;

    assert(this->m_simData.getDimToDiffSpecies().empty());


    const AbstractBoundaryConditions& bc = *this->m_simData.getBoundaryConditions();
    const std::vector<std::string>& compNames = this->m_simData.getCompNames();

    const Box& domain = this->m_palabosData.permMask->getBoundingBox();
    const auto end = DataAccess::end(domain);
    for(auto it = DataAccess::begin(domain); it < end; ++it)
    {
        const auto pos = *it;

        IPPVector3DLong vec;
        PalabosConversionTools::convertToIPP(pos, vec);

        if(isInBoundaryConditions(compNames, bc, vec))
        {
            m_boundaryNodes.dots.push_back(pos);
        }
    }

    using SetValue = SetDotsToScalar_S<MaskType, dim>;
    plb::applyProcessingFunctional(new SetValue(PF_isPermeable), m_boundaryNodes,
                                   *this->m_palabosData.permMask);




    const CopyAllNodes<DotList> func(m_boundaryNodes);
    PermNodesInfos currNodes, newNodes;
    this->collectBounceBackInfos(func, this->m_palabosData.getPorosity(),
                                 m_oldNodes, currNodes, newNodes);

    this->updateLatticesBounceBackNodes<CopyAllNodes<DotList>>(newNodes.permNodes,
                                                               newNodes.nonPermNodes);


    DotList nonPermNonInterfaceNodesNew;
    boost::range::set_difference(currNodes.nonPermNodes.dots,
                                 m_oldNodes.nonPermNodes.dots,
                                 std::back_inserter(nonPermNonInterfaceNodesNew.dots));


    // nodes have to be filtered for own domain
    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const auto ownBox = PalabosConversionTools::ToPlb<dim>::convertBox(ownDomain);

    std::vector<double> dummy(this->m_decomp->getLocalNumberOfCells());
    GeometrySync geomSync(*m_geometrySync, dummy);
    syncGeomData(newNodes.permNodes, newNodes.interfaceNodes, nonPermNonInterfaceNodesNew, ownBox, geomSync);

    m_oldNodes.swap(currNodes);
}

template<typename CollectFunc, typename MaskType,
         typename T, typename Field, typename MaskField,
         typename DotList, typename PermInfos>
static void collectPermeableNodes(const Field& porosity,
                                  const MaskField& mask,
                                  const T& porosThreshLower,
                                  const DotList& inertCells,
                                  PermInfos& nodeInfos)
{

    CollectFunc::template collectNodes<MaskType>(porosThreshLower, porosity, mask,
                                                 nodeInfos.permNodes, nodeInfos.nonPermNodes);

    MPI_CHECK_SYNC;

    std::sort(nodeInfos.permNodes.dots.begin(), nodeInfos.permNodes.dots.end());
    nodeInfos.permNodes.dots.erase(std::unique(nodeInfos.permNodes.dots.begin(),
                                               nodeInfos.permNodes.dots.end()),
                                   nodeInfos.permNodes.dots.end());


    // appending static inert cells
    nodeInfos.nonPermNodes.dots.reserve( nodeInfos.nonPermNodes.dots.size() + inertCells.dots.size() );
    nodeInfos.nonPermNodes.dots.insert( nodeInfos.nonPermNodes.dots.end(),
                                        inertCells.dots.begin(),
                                        inertCells.dots.end() );

    std::sort(nodeInfos.nonPermNodes.dots.begin(), nodeInfos.nonPermNodes.dots.end());
    nodeInfos.nonPermNodes.dots.erase(std::unique(nodeInfos.nonPermNodes.dots.begin(),
                                                  nodeInfos.nonPermNodes.dots.end()),
                                      nodeInfos.nonPermNodes.dots.end());



}

template<typename MaskType, typename DotList, typename Field>
static void collectInterfaceNodes(const DotList& inertNodes, const Field& permMask,
                                  DotList& result)
{
    // collect non-inert solid nodes which have at least one permeable neighbor as interface nodes

    static constexpr size_t dim = PlbTypeDeduction::GetDotListDim<DotList>::value;
    using CollectFunc = CollectSolidCoordinatesWithPermNeighbor<MaskType, dim>;
    DotList interfaceNodes;
    plb::applyProcessingFunctional(new CollectFunc(interfaceNodes),
                                   permMask.getBoundingBox(), const_cast<Field&>(permMask));

    // must be sorted for set_difference
    std::sort(interfaceNodes.dots.begin(), interfaceNodes.dots.end());
    interfaceNodes.dots.erase(std::unique(interfaceNodes.dots.begin(),
                                          interfaceNodes.dots.end()),
                              interfaceNodes.dots.end());

    assert(result.dots.empty());
    boost::range::set_difference(interfaceNodes.dots, inertNodes.dots,
                                 std::back_inserter(result.dots));
}

template<typename TransportTraits>
template<typename CollectNewFunc>
void ReactivePalabosTransportModule<TransportTraits>::collectBounceBackInfos(const CollectNewFunc& func,
                                                                             const ScalarField& porosity,
                                                                             const PermNodesInfos& oldNodes,
                                                                             PermNodesInfos& currNodes,
                                                                             PermNodesInfos& newNodes)
{
    auto& mask = *this->m_palabosData.permMask;

    const Scalar porosThresh = this->m_conf->porosLow;
    const DotList& inertCells = this->m_palabosData.getInertSolidCells();
    collectPermeableNodes<CollectNewFunc, MaskType>(porosity, mask, porosThresh, inertCells, currNodes);

    assert(checkMinPoros(currNodes.permNodes));

    func.collect(oldNodes.nonPermNodes,
                 currNodes.permNodes, currNodes.nonPermNodes,
                 newNodes.permNodes, newNodes.nonPermNodes);

    PermFlagUtils::updatePermMask<MaskType>(currNodes.permNodes, currNodes.nonPermNodes, mask);

    using SetCounter = SetPermNeighCounter<MaskType, dim>;
    plb::applyProcessingFunctional(new SetCounter, mask.getBoundingBox(), mask);

    collectInterfaceNodes<MaskType>(inertCells, mask, currNodes.interfaceNodes);

    boost::range::set_difference(currNodes.interfaceNodes.dots, oldNodes.interfaceNodes.dots,
                                 std::back_inserter(newNodes.interfaceNodes.dots));
}

template<typename TransportTraits>
template<typename CollectDistributeNewFunc>
void ReactivePalabosTransportModule<TransportTraits>::updateGeometry(const CollectDistributeNewFunc &func)
{
    MPI_CHECK_SYNC;

    assert(this->m_simData.getDimToDiffSpecies().empty());

    // we need bounce back althoug we have TRT for very low diffusivity:
    // when porosity is very low it can be lower than minPoros for

    PermNodesInfos currNodes, newNodes;
    this->collectBounceBackInfos(func, this->m_palabosData.getPorosity(),
                                 m_oldNodes, currNodes, newNodes);

    this->updateLatticesBounceBackNodes<CollectDistributeNewFunc>(newNodes.permNodes,
                                                                  newNodes.nonPermNodes);


    DotList nonPermNonInterfaceNodesOld;
    boost::range::set_difference(m_oldNodes.nonPermNodes.dots,
                                 m_oldNodes.interfaceNodes.dots,
                                 std::back_inserter(nonPermNonInterfaceNodesOld.dots));

    DotList nonPermNonInterfaceNodesCurr;
    boost::range::set_difference(currNodes.nonPermNodes.dots,
                                 currNodes.interfaceNodes.dots,
                                 std::back_inserter(nonPermNonInterfaceNodesCurr.dots));

    DotList nonPermNonInterfaceNodesNew;
    boost::range::set_difference(nonPermNonInterfaceNodesCurr.dots,
                                 nonPermNonInterfaceNodesOld.dots,
                                 std::back_inserter(nonPermNonInterfaceNodesNew.dots));


    static constexpr size_t dim = PlbTypeDeduction::GetDotListDim<DotList>::value;

    // nodes have to be filtered for own domain
    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    Dot plbLocation;
    PalabosConversionTools::toPlb(location, plbLocation);
    const auto ownBox = PalabosConversionTools::ToPlb<dim>::convertBox(ownDomain);

    std::vector<double>& concRaw = this->m_simData.getPostReacConc();
    const av::array_view<double, dim+1> concView = ConcArrayViewTool::makeConcView<dim>(this->m_simData, concRaw);

    for(size_t iComp = 0; iComp < this->m_palabosData.diffLattices.size(); ++iComp)
    {
        DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];

        using Descriptor = typename DIFFUSION_DESCRIPTOR<Scalar>;
        ScalarFieldPtr conc = plb::computeExternalScalar(lattice, Descriptor::ExternalField::scalarBeginsAt);
        auto concSlice = concView[iComp];
        using GetConc = GetScalarsToArray<Scalar, decltype(concSlice), dim>;
        plb::applyProcessingFunctional(new GetConc(concSlice, plbLocation),
                                       conc->getBoundingBox(), *conc);
    }

    GeometrySync geomSync(*m_geometrySync, concRaw);
    syncGeomData(newNodes.permNodes, newNodes.interfaceNodes, nonPermNonInterfaceNodesNew, ownBox, geomSync);

    for(size_t iComp = 0; iComp < this->m_palabosData.diffLattices.size(); ++iComp)
    {
        DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];


        using SetFunc = SetScalarsFromArray<const double, dim>;

        const auto concSlice = concView[iComp];
        ScalarField concField = this->m_plbFactory->template createMultiScalarField<Scalar>();
        plb::setToFunction(concField, concField.getBoundingBox(),
                           SetFunc(concSlice, location));

        using Descriptor = typename DIFFUSION_DESCRIPTOR<Scalar>;
        plb::setExternalScalar(lattice, lattice.getBoundingBox(),
                               Descriptor::ExternalField::scalarBeginsAt, concField);
    }

    m_oldNodes.swap(currNodes);
}

template<typename TransportTraits>
template<typename DistributeFunc>
void ReactivePalabosTransportModule<TransportTraits>::updateLatticesBounceBackNodes(const DotList& newPermNodes,
                                                                                    const DotList& newNonPermNodes)
{

    HydrodynamicLatticePtr flowLattice = this->m_palabosData.flowLattice;
    if(flowLattice)
    {
        using BounceBack = plb::BounceBack<Scalar, HYDRODYNAMIC_DESCRIPTOR>;
        plb::defineDynamics(*flowLattice, newNonPermNodes,
                            new BounceBack(0.0));

        using Dynamics = plb::Dynamics<Scalar, HYDRODYNAMIC_DESCRIPTOR>;
        const Dynamics& backgroundDyn = flowLattice->getBackgroundDynamics();
        plb::defineDynamics(*flowLattice, newPermNodes, backgroundDyn.clone());

        using SetDens = DotInitEquilibriumConcentration_L<Scalar, HYDRODYNAMIC_DESCRIPTOR>;
        const AbstractFlowFunctor& flowFunc = *this->m_conf->flowFunc;
        double rho;
        IPPVector3D u;
        flowFunc.getDefault(rho, u);
        plb::applyProcessingFunctional(new SetDens(rho), newPermNodes, *flowLattice);
    }

    std::vector<DiffusionLatticePtr>& diffLattices = this->m_palabosData.diffLattices;
    for(size_t iComp = 0; iComp < diffLattices.size(); ++iComp)
    {
        DiffusionLattice& lattice = *diffLattices[iComp];

        using BounceBack = plb::BounceBack<Scalar, DIFFUSION_DESCRIPTOR>;
        plb::defineDynamics(lattice, newNonPermNodes,
                            new BounceBack);

        using Dynamics = plb::Dynamics<Scalar, DIFFUSION_DESCRIPTOR>;
        const Dynamics& backgroundDyn = lattice.getBackgroundDynamics();
        plb::defineDynamics(lattice, newPermNodes, backgroundDyn.clone());
    }

    // fix the concentrations of new solid nodes and their neighbor
    for(size_t iComp = 0; iComp < diffLattices.size(); ++iComp)
    {
        DiffusionLattice &lattice = *diffLattices[iComp];

        // do not push out water concentration
        if(iComp > 0)
        {
            DistributeFunc::template apply<Scalar, DIFFUSION_DESCRIPTOR,
                    MaskType>(lattice, *this->m_palabosData.permMask);
        }

        using ZeroOperation = ZeroBounceBackPopulations<Scalar, DIFFUSION_DESCRIPTOR>;
        plb::applyProcessingFunctional(new ZeroOperation, lattice.getBoundingBox(), lattice);
    }

}


template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::saveCheckpoint() const
{
    this->saveLatticeCheckpoint();

    const std::string folderName = CreateFileName::create(this->m_simData.getIteration());
    const boost::filesystem::path relativePath = this->m_conf->checkpointDir / folderName;

    const boost::filesystem::path permPath = relativePath / ("perm_mask" +
                                                             IPPConstants::s_dumpExtension);
    plb::saveBinaryBlock(*this->m_palabosData.permMask, permPath.string());

}

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::loadCheckpoint(const size_t iteration)
{

    const std::string folderName = CreateFileName::create(this->m_simData.getIteration());
    const boost::filesystem::path relativePath = this->m_conf->checkpointDir / folderName;

    const boost::filesystem::path permPath = relativePath / ("perm_mask" +
                                                             IPPConstants::s_dumpExtension);
    if(boost::filesystem::exists(permPath))
    {
        plb::loadBinaryBlock(*this->m_palabosData.permMask, permPath.string());
    }
    else
    {
        plb::pcout << "WARNING: no permMask file existing: " << permPath.string() << std::endl
                   << "try to continue with initial permMask..." << std::endl;
    }


    this->loadLatticeCheckpoint(iteration);
    this->updatePostTransportConcField();


    this->updateFieldPorosity();
    this->updateFieldTransportScalar();


    ScalarField& porosity = this->m_palabosData.getPorosity();
    PermNodesInfos currNodes, newNodes;
    this->collectBounceBackInfos(CollectDiff(), porosity,
                                 m_oldNodes, currNodes, newNodes);

    this->updateLatticesBounceBackNodes<CopyAllNodes<DotList>>(newNodes.permNodes,
                                                               newNodes.nonPermNodes);
    m_oldNodes.swap(currNodes);
}


template<typename TransportTraits>
bool ReactivePalabosTransportModule<TransportTraits>::checkMinPoros(const DotList& toCheck)
{
    using SetValue = SetDotsToScalar_S<Scalar, dim>;
    ScalarField maskTmp = this->m_plbFactory->template createMultiScalarField<Scalar>(0);
    const int maskVal = 1;
    plb::applyProcessingFunctional(new SetValue(maskVal), toCheck, maskTmp);

    const auto mask = plb::copyConvert<Scalar, int>(maskTmp);

    ScalarField& porosity = this->m_palabosData.getPorosity();
    const Scalar minVal = plb::computeMin(porosity, *mask, maskVal);

    const Scalar fac = minVal / this->m_conf->porosLow;

    if(toCheck.dots.size() > 0)
    {

        if(fac < 0.9)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    else
    {
        return true;
    }
}

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::prepareResults()
{
    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;

    const auto concView = ConcArrayViewTool::makeConcView<dim>(this->m_simData,
                                                               this->m_simData.getPostReacConc());

    for(size_t iComp = 0; iComp < this->m_palabosData.postReactionConcentrations.size(); ++iComp)
    {
        const auto concSlice = concView[iComp];
        ScalarField& conc = this->m_palabosData.postReactionConcentrations[iComp];

        using Sync = SetScalarsFromArray<const double, dim>;
        plb::setToFunction(conc, conc.getBoundingBox(), Sync(concSlice, location));

    }
}


template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::initLatticePorosity()
{
    ScalarField& physicalPorosity = this->m_palabosData.getPorosity();

    const size_t nComps = this->m_palabosData.diffLattices.size();
    for(size_t iComp = 0; iComp < nComps; ++iComp)
    {
        DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];

        ScalarField& porosity = iComp == 0 ? *m_unityField : physicalPorosity;

        using DiffusionDesc = typename DIFFUSION_DESCRIPTOR<Scalar>;
        plb::setExternalScalar(lattice, lattice.getBoundingBox(),
                               DiffusionDesc::ExternalField::porosityBeginsAt, porosity);
    }
}

template<typename T, template<typename U> class Descriptor>
class UpdateTransportScalar : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD_LS<T,Descriptor>::value
{
public:
    UpdateTransportScalar() = default;


    virtual UpdateTransportScalar<T,Descriptor>* clone() const override
    {
        return new UpdateTransportScalar<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::nothing;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulk;
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD_LS<T,Descriptor>::value;

    using Traits =
    PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;

    using Box = typename Traits::template argument<1>::type;
    using Lattice = typename Traits::template argument<2>::type;
    using ScalarField = typename Traits::template argument<3>::type;

public:
    virtual void process(Box domain, Lattice& lattice, ScalarField& transScalarField) override
    {
        static constexpr size_t dim = Descriptor<T>::d;
        using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

        const Dot offset = plb::computeRelativeDisplacement(lattice, transScalarField);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const Dot& pos = *it;

            plb::Cell<T, Descriptor>& cell = DataAccess::get(lattice, pos);
            T& transScalar = *cell.getExternal(Descriptor<T>::ExternalField::transportScalarBeginsAt);

            using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;
            const DotList neighList = DotFindNeighbors::find(pos);
            assert(neighList.getN() == dim * 2);

            transScalar = DataAccess::get(transScalarField, pos + offset);

            for(const Dot& neigh : neighList.dots)
            {
                const T& transScalarNeigh = DataAccess::get(transScalarField, neigh + offset);
                transScalar = std::min(transScalar, transScalarNeigh);
            }
        }
    }

};


template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::updateLatticeTransportScalar()
{
    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    const IPPVector3DLong ownSize = ownDomain.getSize();

    const auto bounds = av::MakeBounds<dim>::get(ownSize[0], ownSize[1], ownSize[2]);
    const av::array_view<const Scalar, dim> transScalarView(this->m_transportScalar, bounds);

    ScalarField transScalar = this->m_plbFactory->template createMultiScalarField<Scalar>(2.0);

    using Sync = SetScalarsFromArray<const Scalar, dim>;
    plb::setToFunction(transScalar, transScalar.getBoundingBox(),
                       Sync(transScalarView, location));

    const size_t nComps = this->m_palabosData.diffLattices.size();
    for(size_t iComp = 0; iComp < nComps; ++iComp)
    {
        DiffusionLattice& lattice = *this->m_palabosData.diffLattices[iComp];

        using Update = UpdateTransportScalar<Scalar, DIFFUSION_DESCRIPTOR>;
        plb::applyProcessingFunctional(new Update, lattice.getBoundingBox(),
                                       lattice, transScalar);
    }
}

template<typename TransportTraits>
void ReactivePalabosTransportModule<TransportTraits>::updateFieldPorosity()
{
    using ArrView = av::array_view<const double, dim>;

    const IPPBox3DLong& ownDomain = this->m_decomp->getOwnDomain();
    const IPPVector3DLong& location = ownDomain.lower;
    const IPPVector3DLong size = ownDomain.getSize();
    const av::bounds<dim> bounds = av::MakeBounds<dim>::get(size[0], size[1], size[2]);

    ScalarField& physicalPorosity = this->m_palabosData.getPorosity();

    const std::vector<double>& localPoros = this->m_simData.getPorosity();
    const ArrView porosView(localPoros, bounds);

    m_oldPorosity->swap(physicalPorosity);
    plb::setToFunction(physicalPorosity, this->m_transportSyncEnabledBox,
                       SetScalarsFromArray<const double, dim>(porosView, location));
}

} // end of namespace IPP


#endif // REACTIVEPALABOSTRANSPORTMODULE_HH
