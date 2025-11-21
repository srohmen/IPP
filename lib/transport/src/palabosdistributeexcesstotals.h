#ifndef PALABOSDISTRIBUTEEXCESSTOTALS_H
#define PALABOSDISTRIBUTEEXCESSTOTALS_H

#include "abstractdistributeexcesstotals.h"

#include "plbtypededuction.h"
#include "palabosobjectfactory.h"

namespace IPP
{

class FieldDecomposition;

template<typename T, typename Field>
class PalabosDistributeExcessTotals : public AbstractDistributeExcessTotals
{    
private:
    static constexpr size_t dim = PlbTypeDeduction::GetFieldDim<T, Field>::value;
    using Factory = PalabosObjectFactory<dim>;

public:
    PalabosDistributeExcessTotals(const FieldDecomposition& decomp,
                                  const Factory& factory);

    virtual void init() override;

    virtual void setnComps(const std::size_t nComps) override;
    virtual void updatePorosity(const std::vector<double>& porosity) override;

    virtual void execute(const std::vector<CellTotalsDiff>& diffPerCell,
                         std::vector<double>& distributedSource) override;


private:
    const FieldDecomposition& m_decomp;
    const Factory& m_factory;

    size_t m_nComps;
    std::unique_ptr<Field> m_porosity;
};

}


#ifdef EXPLICIT_INSTANTS
#include <palabos/multiBlock/multiDataField2D.h>
#include <palabos/multiBlock/multiDataField3D.h>
extern template class IPP::PalabosDistributeExcessTotals<double, plb::MultiScalarField2D<double> >;
extern template class IPP::PalabosDistributeExcessTotals<double, plb::MultiScalarField3D<double> >;
#else
#include "palabosdistributeexcesstotals.hh"
#endif


#endif // PALABOSDISTRIBUTEEXCESSTOTALS_H
