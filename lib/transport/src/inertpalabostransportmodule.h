#ifndef INERTPALABOSTRANSPORTMODULE_H
#define INERTPALABOSTRANSPORTMODULE_H

#include "palabostransportmodule.h"

namespace IPP
{

template<typename TransportTraits>
class InertPalabosTransportModule : public PalabosTransportModule<TransportTraits>
{
    using BC = PalabosTransportModule<TransportTraits>;
    using BC::PalabosTransportModule;

    using BC::dim;
    using typename BC::DotList;
    using typename BC::Box;
    using typename BC::Scalar;
    using typename BC::ScalarField;
    using typename BC::ScalarFieldPtr;
    using typename BC::ScalarFieldSharedPtr;
    using typename BC::PlbFactory;
    using typename BC::DiffusionLattice;
    using typename BC::DiffusionLatticePtr;
//    using typename BC::DiffusionDesc;
//    template <typename T>
//    using DiffusionDescT = typename TransportTraits:: template DiffusionDescriptorT<T>;

    using typename BC::TransportData;

public:
    virtual void setGeomSync(AbstractGeometrySync* geomSync) override;

    virtual void updatePostReactionState() override;
    virtual void collectData() override;

    virtual void saveCheckpoint() const override;
    virtual void loadCheckpoint(const size_t iteration) override;

private:
    virtual void initModule() override;
    virtual void addPostStreamOperation(DiffusionLattice &lattice) override;
    virtual void prepareResults() override;

    void initBounceBackNodes();
    void collectPermNodes(DotList &permNodes, DotList& nonPermNodes);


    void updateFieldPorosity();
    void initLatticeTransportProperties();

    void initLatticeConcentrations();
};




}

#ifndef EXPLICIT_INSTANTS
#include "inertpalabostransportmodule.hh"
#endif

#endif // INERTPALABOSTRANSPORTMODULE_H
