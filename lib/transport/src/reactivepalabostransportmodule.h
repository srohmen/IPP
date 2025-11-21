#ifndef REACTIVEPALABOSTRANSPORTMODULE_H
#define REACTIVEPALABOSTRANSPORTMODULE_H

#include "palabostransportmodule.h"

namespace IPP
{

class AbstractGeometrySync;

template<typename TransportTraits>
class ReactivePalabosTransportModule : public PalabosTransportModule<TransportTraits>
{
    using BC = PalabosTransportModule<TransportTraits>;

public:
    ReactivePalabosTransportModule(TransportModuleConfigPtr& config);

    virtual void setGeomSync(AbstractGeometrySync* geomSync) override;

    virtual void collectData() override;
    virtual void updatePostReactionState() override;

    virtual void saveCheckpoint() const override;
    virtual void loadCheckpoint(const size_t iteration) override;
private:

    static constexpr size_t dim = BC::dim;
    using typename BC::Scalar;
    using typename BC::Dot;
    using typename BC::DotList;
    using typename BC::Box;
    using typename BC::ArrayView;
    using typename BC::ScalarField;
    using typename BC::ScalarFieldPtr;
    using typename BC::ScalarFieldSharedPtr;
    template<typename T>
    using ScalarFieldT = PlbTypeDeduction::GetScalarFieldT_XD<dim>;
    using typename BC::DiffusionLattice;
    using typename BC::DiffusionLatticePtr;
//    using typename BC::DiffusionDesc;

//    template <typename T>
//    using DiffusionDescT = typename TransportTraits:: template DiffusionDescriptorT<T>;

//    template <typename T>
//    using HydrodynamicDescT = typename TransportTraits:: template HydrodynamicDescriptorT<T>;

    using typename BC::PlbFactory;
    using typename BC::HydrodynamicLatticePtr;

    using typename BC::TransportData;
    using MaskType = typename TransportData::MaskType;


    virtual void initModule() override;
    virtual void addPostStreamOperation(DiffusionLattice &lattice) override;
    virtual void prepareResults() override;

    void initLatticeConcentrations();
    void updateLatticeTransportScalar();
    void initLatticePorosity();
    void updateFieldPorosity();
    void syncDataIntoLattices();
    bool hasGeomChanged();

    void initBounceBackNodes();

    template<typename CollectDistributeNewFunc>
    void updateGeometry(const CollectDistributeNewFunc& func);

    template<typename DistributeFunc>
    void updateLatticesBounceBackNodes(const DotList& newPermNodes,
                                       const DotList& newNonPermNodes);


    void updateNeighInfos();

    void addAdditionalSource();



    void initFields();    
    void updatePostTransportConcField();


    void calcSourceTerm(ScalarField& nOld,
                        ScalarField& nNew,
                        const size_t iComp,
                        ScalarField& source);



    struct PermNodesInfos
    {
        DotList permNodes;
        DotList nonPermNodes;
        DotList interfaceNodes;

        void clear()
        {
            permNodes.dots.clear();
            nonPermNodes.dots.clear();
            interfaceNodes.dots.clear();
        }

        void swap(PermNodesInfos& other)
        {
            permNodes.dots.swap(other.permNodes.dots);
            nonPermNodes.dots.swap(other.nonPermNodes.dots);
            interfaceNodes.dots.swap(other.interfaceNodes.dots);
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & permNodes;
            ar & nonPermNodes;
            ar & interfaceNodes;
        }

    };

    template<typename CollectNewFunc>
    void collectBounceBackInfos(const CollectNewFunc &func,
                                const ScalarField &porosity,
                                const PermNodesInfos &oldNodes,
                                PermNodesInfos &currNodes,
                                PermNodesInfos &newNodes);

    void collectPermNodes(const Scalar& porosThresh,
                          const DotList& oldNonPerm, DotList& nonPermNodes,
                          DotList& newPermNodes, DotList& newNonPermNodes);


    void collectNewBounceBackNodes(const DotList &nonPermNodes,
                                   DotList &newNonPermNodes,
                                   DotList &newPermNodes);

    bool checkMinPoros(const DotList &toCheck);


    // needed for source term calculation
    std::unique_ptr<ScalarField> m_oldPorosity;
    std::unique_ptr<ScalarField> m_unityField;

    using NTensorField = typename PlbTypeDeduction::GetMultiNTensorField_XD<Scalar,dim>::value;
    std::unique_ptr<NTensorField> m_capilPorosVolumes;


    PermNodesInfos m_oldNodes;

    AbstractGeometrySync* m_geometrySync;

    // only used on init and checkpoint loading
    DotList m_boundaryNodes;

    std::vector<double> m_negConcBuffer;


};



}

#ifndef EXPLICIT_INSTANTS
#include "reactivepalabostransportmodule.hh"
#endif

#endif // REACTIVEPALABOSTRANSPORTMODULE_H
