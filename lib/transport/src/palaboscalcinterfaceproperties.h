#ifndef PALABOSCALCINTERFACEPROPERTIES_H
#define PALABOSCALCINTERFACEPROPERTIES_H

#include <memory>

#include "abstractcalcinterfaceproperties.h"
#include "plbtypededuction.h"
#include "palabosobjectfactory.h"
#include "ippvector.h"

namespace IPP
{

class FieldDecomposition;

template<typename T, typename MaskType, size_t dim>
class PalabosCalcInterfaceProperties : public AbstractCalcInterfaceProperties
{
    using DotList = typename PlbTypeDeduction::GetDotListXD<dim>::value;
    using MaskField = typename PlbTypeDeduction::GetMultiScalarField_XD<MaskType, dim>::value;

public:
    PalabosCalcInterfaceProperties(const T& weight,
                                   const FieldDecomposition& decomp,
                                   const DotList& interfaceNodes);

    using Factory = PalabosObjectFactory<dim>;
    void initFields(const MaskField* permableMask, const Factory& factory);

    virtual void setnComps(const size_t nComps) override;

    virtual void execute(const std::vector<double>& conc,
                         const std::vector<double>& porosity,
                         std::vector<double>& avgConc,
                         std::vector<double>& avgPoros) const override;

private:
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;
    using ScalarField = typename PlbTypeDeduction::GetMultiScalarField_XD<T,dim>::value;
    using NTensorField = typename PlbTypeDeduction::GetMultiNTensorField_XD<T,dim>::value;
    using MultiBlock = typename PlbTypeDeduction::GetMultiBlock_XD<dim>::value;

    void averagePorosity(const IPPVector3DLong& location,
                         const Dot& dotOrigin,
                         const IPPVector3DLong& ownSize,
                         const std::vector<double>& porosity,
                         ScalarField& avgPorosField,
                         std::vector<double>& avgPoros) const;

    void averageConc(const IPPVector3DLong& location,
                     const Dot& dotOrigin,
                     const IPPVector3DLong& ownSize,
                     const ScalarField& avgPorosField,
                     const std::vector<double>& conc,
                     std::vector<double>& avgConc) const;


    const T m_weight;
    const FieldDecomposition& m_decomp;
    const DotList& m_interfaceNodes;
    const MaskField* m_permableMask;

    size_t m_nComps;


    std::unique_ptr<ScalarField> m_porosityField;
    std::unique_ptr<ScalarField> m_avgPorosField;
    std::unique_ptr<NTensorField> m_concField;



};

}

#ifdef EXPLICIT_INSTANTS
#include "transporttypes.h"
extern template class IPP::PalabosCalcInterfaceProperties<IPP::TransportFloatingPointType, int, 2>;
extern template class IPP::PalabosCalcInterfaceProperties<IPP::TransportFloatingPointType, int, 3>;
#else
#include "palaboscalcinterfaceproperties.hh"
#endif


#endif // PALABOSCALCINTERFACEPROPERTIES_H
