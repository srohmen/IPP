#include "phreeqccalcporosity.h"

#include <cassert>
#include "ippexception.h"

#include "array_view.h"
#include "arrayviewiterator.h"
#include "arrayviewiteratorimpl.h"
#include "array_view_strm.h"

#include "abstractporositycalc.h"
#include "simpleporosityinfos.h"

namespace IPP
{

void PhreeqcCalcPorosity::calc(const PhaseNameToInfos& phaseInfos,
                               const std::vector<double>& solidFractions,
                               const std::vector<double>& precipMol,
                               const std::vector<double>& inertVolFrac,
                               const AbstractPorosityCalc& porosCalc,
                               std::vector<SimplePorosityInfosPtr>& porosInfos)
{
    const size_t nDomains = inertVolFrac.size();
    const size_t nPhases = phaseInfos.size();
    porosInfos.resize(nDomains); // check with move
    assert(solidFractions.size() % nPhases == 0);

    const av::bounds<2> solidFracBounds = {static_cast<std::ptrdiff_t>(nDomains),
                                           static_cast<std::ptrdiff_t>(nPhases)};
    const av::array_view<const double, 2> solidFracView(solidFractions, solidFracBounds);



    const av::bounds<2> amountBounds = {static_cast<std::ptrdiff_t>(nDomains),
                                        static_cast<std::ptrdiff_t>(nPhases)};
    const av::array_view<const double, 2> amountView(precipMol, amountBounds);



    for(size_t iDomain = 0; iDomain < nDomains; ++iDomain)
    {
        AbstractPorosityCalc::PorosityCalcInput input;
        input.inertFraction = inertVolFrac[iDomain];


        auto solidFracDomainViewTmp =
                solidFracView.section({ (std::ptrdiff_t)iDomain, 0}, {1, (std::ptrdiff_t)nPhases } );
        input.volFractions.begin = ArrayViewIterator(new ArrayViewIteratorImpl(solidFracDomainViewTmp, std::begin(solidFracDomainViewTmp.bounds())));
        input.volFractions.end = ArrayViewIterator(new ArrayViewIteratorImpl(solidFracDomainViewTmp, std::end(solidFracDomainViewTmp.bounds())));

        auto amountDomainViewTmp =
                amountView.section({ (std::ptrdiff_t)iDomain, 0 }, { 1, (std::ptrdiff_t)nPhases } );
        input.phaseAmounts.begin = ArrayViewIterator(new ArrayViewIteratorImpl(amountDomainViewTmp, std::begin(amountDomainViewTmp.bounds()) ));
        input.phaseAmounts.end = ArrayViewIterator(new ArrayViewIteratorImpl(amountDomainViewTmp, std::end(amountDomainViewTmp.bounds())   ));


        SimplePorosityInfosPtr result = porosCalc.calc(input);
        porosInfos[iDomain] = std::move(result);
    }
}

}
