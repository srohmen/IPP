#include "configurephreeqccellsactivity.h"

#include "localaccessphreeqcrm.h"

#include "isinsidevolume.h"
#include "significantchangeprediction.h"
#include "anyelementunbuffered.h"
#include "handleinterfaceswitch.h"

#include "arrayviewiterator.h"

#include "mpitools.h"
#include "phreeqcconstants.h"
#include "indexhelper.h"

#include "scopedfloatingpointexception.h"

namespace IPP
{

namespace TestChains
{

class TestChain
{
public:
    virtual ~TestChain()
    {

    }

    virtual bool hasSignificantChange(const size_t iCell) const = 0;
};



class IsConcDiffHighTest : public TestChain
{
public:
    IsConcDiffHighTest(const std::vector<double>& lastPhreeqcConc,
                       const std::vector<double>& currConc,
                       const double& tolerance,
                       const size_t nCells,
                       const size_t nComps)
        : m_lastPhreeqcConc(lastPhreeqcConc)
        , m_currConc(currConc)
        , m_tolerance(tolerance)
        , m_nCells(nCells)
        , m_nComps(nComps)

    {

    }

    virtual bool hasSignificantChange(const size_t iCell) const
    {
        for(size_t iComp = 0; iComp < m_nComps; ++iComp)
        {
            if(iComp == PhreeqcConstants::s_primaryCompsBegin - 1)
            {
                continue;
            }

            const size_t index = iComp * m_nCells + iCell;

            const double& currConc = m_currConc[index];
            assert(currConc >= 0.0);

            const double& oldConc = m_lastPhreeqcConc[index];
            assert(oldConc >= 0.0);

            const double diff = currConc - oldConc;
            const double absDiff = std::abs(diff);

            if(absDiff > 1.0E-9)
            {
                const double diffRel = absDiff / oldConc;
                if(diffRel > m_tolerance)
                {
                    return true;
                }
            }

        }

        return false;
    }

private:
    const std::vector<double>& m_lastPhreeqcConc;
    const std::vector<double>& m_currConc;

    const double m_tolerance;
    const size_t m_nCells;
    const size_t m_nComps;

};


class IsResidualConcLowTest : public TestChain
{
public:
    IsResidualConcLowTest(const std::vector<double>& lastPhreeqcConc,
                          const std::vector<double>& postTransConc,
                          const std::vector<double>& cummuConcDiff,
                          const double& thresh,
                          const size_t nCells,
                          const size_t nComps)
        : m_lastPhreeqcConc(lastPhreeqcConc)
        , m_postTransportConc(postTransConc)
        , m_cummuConcDiff(cummuConcDiff)
        , m_thresh(thresh)
        , m_nCells(nCells)
        , m_nComps(nComps)
    {

    }

    bool hasSignificantChange(const size_t iCell) const
    {
        // prevent undershooting
        const size_t chargeIndex = PhreeqcConstants::s_primaryCompsBegin - 1;

        for(size_t iComp = 0; iComp < m_nComps; ++iComp)
        {
            const size_t index = IndexHelper::toPhreeqcRawIndex(m_nCells, iCell, iComp);
            const double& cummDiff = m_cummuConcDiff[index];


            const double& postTransConc = m_postTransportConc[index];
            const double realConc = postTransConc + cummDiff;

            const double& oldConc = m_lastPhreeqcConc[index];
            const double realConcDiffPrediction = oldConc * m_thresh;

            if(iComp != chargeIndex)
            {
                if(realConc <= realConcDiffPrediction)
                {
                    return true;
                }
            }
            else
            {
                const double variance = (realConc - oldConc) / oldConc;
                const double chargeThresh = 1.0;
                if(std::abs(variance) > chargeThresh)
                {
                    return true;
                }


                // prevent sign switch
                if(std::signbit(realConc) != std::signbit(oldConc))
                {
                    //                        std::cout << "sign switch in charge: " << iCell
                    //                                  << "\t" << oldConc << " -> " << realConc  << std::endl;
                    return true;
                }

            }
        }

        return false;
    }

private:
    const std::vector<double>& m_lastPhreeqcConc;
    const std::vector<double>& m_postTransportConc;
    const std::vector<double>& m_cummuConcDiff;
    const double m_thresh;
    const size_t m_nCells;
    const size_t m_nComps;

};


class SignificantChangeTest : public TestChain
{
public:
    SignificantChangeTest(const PhaseIdPerElement& phaseIdPerElement,
                          const std::vector<double>& lastPhreeqcPhaseAmount,
                          const size_t nPhases,
                          const SignificantChangePrediction& changePredict,
                          const size_t nCells,
                          const size_t nComps,
                          const std::vector<double>& lastPhreeqcTotalAmount,
                          const std::vector<double>& concDiffCummulative,
                          const std::vector<double>& lastPhreeqcConc,
                          const double& thresh)
        : m_anyUnbuffered(phaseIdPerElement)
        , m_lastPhreeqcPhaseAmount(lastPhreeqcPhaseAmount)
        , m_nPhases(nPhases)
        , m_changePredict(changePredict)
        , m_nCells(nCells)
        , m_nComps(nComps)
        , m_lastPhreeqcTotalAmount(lastPhreeqcTotalAmount)
        , m_concDiffCummulative(concDiffCummulative)
        , m_lastPhreeqcConc(lastPhreeqcConc)
        , m_thresh(thresh)
    {

    }

    virtual ~SignificantChangeTest()
    {

    }

    virtual bool hasSignificantChange(const size_t iCell) const
    {

        const AnyElementUnbuffered::Iterator phasesBegin =
                m_lastPhreeqcPhaseAmount.begin() + iCell * m_nPhases;

        const AnyElementUnbuffered::Iterator phasesEnd =
                phasesBegin + m_nPhases;

        {
            const size_t offset = PhreeqcConstants::s_primaryCompsBegin;
            const size_t nPrimaryComps = m_nComps - offset;
            const size_t primaryIndexBegin = iCell * nPrimaryComps;

            const std::vector<double>::const_iterator currTotalsBegin =
                    m_lastPhreeqcTotalAmount.begin() + primaryIndexBegin;
            const std::vector<double>::const_iterator currTotalsEnd = currTotalsBegin + nPrimaryComps;
            std::vector<double> currTotals(nPrimaryComps);
            std::copy(currTotalsBegin, currTotalsEnd, currTotals.begin());


            // filter only primary comp
            // cummulative must be transposed
            std::vector<double> predictedTotalsMin(nPrimaryComps);
            std::vector<double> predictedTotalsMax(nPrimaryComps);
            for(size_t iComp = 0; iComp < nPrimaryComps; ++iComp)
            {
                const size_t iCompAll = iComp + PhreeqcConstants::s_primaryCompsBegin;
                const size_t srcIndex = IndexHelper::toPhreeqcRawIndex(m_nCells, iCell, iCompAll);

                assert(m_concDiffCummulative.size() > srcIndex);
                const double& cummDiff = m_concDiffCummulative[srcIndex];

                const double& eqConc = m_lastPhreeqcConc[srcIndex];
                const double concDiffAbs = m_thresh * eqConc;

                predictedTotalsMin[iComp] = cummDiff - concDiffAbs;
                predictedTotalsMax[iComp] = cummDiff + concDiffAbs;
            }


            typedef SignificantChangePrediction::Interpolation Interpol;
            const Interpol::Coord& accur = m_changePredict.getAccuracy();
            assert(accur.size() == nPrimaryComps);

            for(size_t iComp = 0; iComp < nPrimaryComps; ++iComp)
            {
                typedef Interpol::Coord::value_type T;
                T& min = predictedTotalsMin[iComp];
                T& max = predictedTotalsMax[iComp];

                const T& curr = currTotals[iComp];

                min += curr;
                max += curr;

                // increase conservatism
                const Interpol::Coord::value_type& acc = accur[iComp];
                assert(acc >= 0.0);
                min -= acc;
                max += acc;
            }


            std::vector<double> currPhases(m_nPhases);
            std::copy(phasesBegin, phasesEnd, currPhases.begin());

            const bool isSignificant = m_changePredict.isSignificant(currPhases,
                                                                     currTotals,
                                                                     predictedTotalsMin,
                                                                     predictedTotalsMax);

            return isSignificant;
        }
    }

private:
    const AnyElementUnbuffered m_anyUnbuffered;
    const std::vector<double>& m_lastPhreeqcPhaseAmount;
    const size_t m_nPhases;

    const SignificantChangePrediction& m_changePredict;
    const size_t m_nCells;
    const size_t m_nComps;
    const std::vector<double>& m_lastPhreeqcTotalAmount;
    const std::vector<double>& m_concDiffCummulative;

    const std::vector<double>& m_lastPhreeqcConc;
    const double m_thresh;
};



}


ConfigurePhreeqcCellsActivity::ConfigurePhreeqcCellsActivity(const double& tolerance,
                                                             const size_t forceCalcIt,
                                                             const SignificantChangePrediction *changePredict)
    : m_isGloballyEnabled(false)
    , m_isGloballyDisabled(false)
    , m_tolerance(tolerance)
    , m_forceCalcIt(forceCalcIt)
    , m_changePredict(changePredict)
    , m_hasSignificancyChecks(false)
{

}

ConfigurePhreeqcCellsActivity::~ConfigurePhreeqcCellsActivity()
{
    for(TestChains::TestChain* test : m_chain)
    {
        delete test;
    }
    m_chain.clear();
}

void ConfigurePhreeqcCellsActivity::init(const size_t nCells, const size_t nComps,
                                         const size_t nPrimaryComps, const size_t nPhases,
                                         const std::vector<double>& currConc,
                                         const PhaseIdPerElement& phaseIdPerElement,
                                         std::vector<char>& isNonPermNonInterfaceCell,
                                         std::vector<char>& isInterfaceCell)
{
    m_interfaceSwitch.reset(new HandleInterfaceSwitch(isNonPermNonInterfaceCell, isInterfaceCell));

    m_staticDisabledCells.resize(nCells, false);
    m_lastPhreeqcEqConc.resize(nCells * nComps, -1.0);

    m_lastPhreeqcPhaseAmount.resize(nCells * nPhases, -1.0);
    m_lastPhreeqcTotalAmount.resize(nCells * nPrimaryComps, -1.0);
    m_concDiffCummulative.resize(nCells * nComps, 0.0);

    using namespace TestChains;



    assert(m_chain.empty());

    // order controls order to check
    if(m_tolerance > 0.0)
    {
        IsConcDiffHighTest* concTest =
                new IsConcDiffHighTest(m_lastPhreeqcEqConc,
                                       currConc,
                                       m_tolerance,
                                       nCells, nComps);
        m_chain.push_back(concTest);

        m_hasSignificancyChecks = true;
    }

    if(m_changePredict)
    {
        // TODO: each timestep transport can reduce 1/6 (D2Q5) or 1/4 (D3Q7) at maximum if tau = 1
        // check tau dependency!
        // conservative threshold is a residual conc of 1/6 or 1/4 of old eq-conc, respectively
        // FIXME: setting to a conservative value for 2d and 3d simulations

        const double residualThresh = 0.3;


        IsResidualConcLowTest* concTest = new IsResidualConcLowTest(m_lastPhreeqcEqConc,
                                                                    currConc,
                                                                    m_concDiffCummulative,
                                                                    residualThresh,
                                                                    nCells, nComps);
        m_chain.push_back(concTest);


        SignificantChangeTest* predictTest =
                new SignificantChangeTest(phaseIdPerElement,
                                          m_lastPhreeqcPhaseAmount,
                                          nPhases,
                                          *m_changePredict,
                                          nCells,
                                          nComps,
                                          m_lastPhreeqcTotalAmount,
                                          m_concDiffCummulative,
                                          m_lastPhreeqcEqConc,
                                          residualThresh);
        m_chain.push_back(predictTest);

        m_hasSignificancyChecks = true;
    }



}

HandleInterfaceSwitch& ConfigurePhreeqcCellsActivity::getInterfaceSwitch()
{
    return *m_interfaceSwitch;
}

void ConfigurePhreeqcCellsActivity::setStaticDisabledCells(const std::vector<size_t> &cellIDs)
{
    MPI_CHECK_SYNC;

    const size_t nCells = m_staticDisabledCells.size();
    m_staticDisabledCells.clear();
    m_staticDisabledCells.resize(nCells, false);

    for(const size_t id : cellIDs)
    {
        m_staticDisabledCells[id] = true;
    }
}

void ConfigurePhreeqcCellsActivity::setGloballyEnabledState(const bool isEnabled)
{
    m_isGloballyEnabled = isEnabled;
}

void ConfigurePhreeqcCellsActivity::setGloballyDisabledState(const bool isDisabled)
{
    m_isGloballyDisabled = isDisabled;
}

bool ConfigurePhreeqcCellsActivity::mustRecalc(const size_t it, const size_t iCell)
{
    if(m_isGloballyDisabled)
    {
        return false;
    }

    // pretests
    assert(m_staticDisabledCells.size() > iCell);
    if(m_staticDisabledCells[iCell])
    {
        return false;
    }

    if(m_interfaceSwitch->mustRun(iCell) == false)
    {
        return false;
    }


    if(m_isGloballyEnabled)
    {
        return true;
    }

    if(m_forceCalcIt > 0)
    {
        if(it % m_forceCalcIt == 0)
        {
            return true;
        }
    }

    const bool significantChange = hasSignificantChange(iCell);
    return significantChange;
}


bool ConfigurePhreeqcCellsActivity::hasSignificantChange(const size_t iCell)
{
    if(m_hasSignificancyChecks)
    {
        // if any test returns true there was significant change
        for(size_t iTest = 0; iTest < m_chain.size(); ++iTest)
        {
            const TestChains::TestChain* test = m_chain[iTest];
            if(test->hasSignificantChange(iCell) == true)
            {
                return true;
            }

        }
    }
    else
    {
        // if there is not significancy check recalc must be done everytime
        assert(m_chain.empty());
        return true;
    }

    return false;
}


void ConfigurePhreeqcCellsActivity::run(const size_t it,
                                        const std::vector<char>& oldEnabled,
                                        std::vector<size_t>& newEnabled,
                                        std::vector<size_t>& newDisabled)
{
    assert(newEnabled.empty());
    assert(newDisabled.empty());

    const ScopedDisableFloatingPointException fpe;

    const size_t nCells = m_staticDisabledCells.size();

    for(size_t iCell = 0; iCell < nCells; ++iCell)
    {
        const bool mustRecalculate = mustRecalc(it, iCell);

        if(mustRecalculate)
        {
            if(oldEnabled[iCell] == false)
            {
                newEnabled.push_back(iCell);
            }
        }
        else
        {
            if(oldEnabled[iCell])
            {
                newDisabled.push_back(iCell);
            }
        }
    }
}

std::vector<double> &ConfigurePhreeqcCellsActivity::getLastPhreeqcEqConc()
{
    return m_lastPhreeqcEqConc;
}

std::vector<double>&ConfigurePhreeqcCellsActivity::getLastPhreeqcPhaseAmount()
{
    return m_lastPhreeqcPhaseAmount;
}

std::vector<double>& ConfigurePhreeqcCellsActivity::getLastPhreeqcTotalAmount()
{
    return m_lastPhreeqcTotalAmount;
}

std::vector<double> &ConfigurePhreeqcCellsActivity::getConcDiffCummulative()
{
    return m_concDiffCummulative;
}


}
