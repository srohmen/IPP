#ifndef PALABOSTRANSPORTDATA_H
#define PALABOSTRANSPORTDATA_H

#include <string>
#include <memory>

#include "array_view.h"
#include "transportmoduleconfig.h"
#include "palabosobjectfactory.h"
#include "setupblockperiodicity.h"

namespace IPP
{


template<typename TransportTraits>
class PalabosTransportData
{
public:

    static constexpr size_t dim = TransportTraits::dim;
    using Scalar = typename TransportTraits::Scalar;
    using PlbFactory = PalabosObjectFactory<dim>;

    void init(const PlbFactory& factory, const size_t nComp)
    {
        porosity.reset(factory.template createMultiScalarFieldPtr<Scalar>(1.0));

        assert(postTransportConcentrations.empty());
        assert(postReactionConcentrations.empty());

        postTransportConcentrations.resize(nComp, factory.template createMultiScalarField<Scalar>(-1.0));
        postReactionConcentrations.resize(nComp, factory.template createMultiScalarField<Scalar>(-1.0));
    }

    using ScalarField = typename TransportTraits::ScalarField;

    ScalarField& getPorosity()
    {
        return *porosity;
    }

    const ScalarField& getPorosity() const
    {
        return *porosity;
    }

    using DotList = typename TransportTraits::DotList;

    DotList& getInertSolidCells()
    {
        return inertCells;
    }

    // used in each iteration, but updated only when geometry changed
    using MaskType = int;
    using ScalarFieldMask = typename PlbTypeDeduction::GetMultiScalarField_XD<MaskType, dim>::value;
    std::unique_ptr<ScalarFieldMask> permMask;


public:
    std::vector<std::string> componentNames;

    using HydrodynamicLatticePtr = typename TransportTraits::HydrodynamicLatticePtr;
    HydrodynamicLatticePtr flowLattice;

    using DiffusionLatticePtr = typename TransportTraits::DiffusionLatticePtr;
    using DiffusionLatticePtrVec = std::vector<DiffusionLatticePtr>;
    DiffusionLatticePtrVec diffLattices;

    std::vector<ScalarField> postTransportConcentrations;
    std::vector<ScalarField> postReactionConcentrations;

private:
    std::unique_ptr<ScalarField> porosity;
    DotList inertCells;

};

}

#endif // PALABOSTRANSPORTDATA_H
