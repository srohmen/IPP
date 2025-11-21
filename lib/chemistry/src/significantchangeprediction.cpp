#include "significantchangeprediction.h"


namespace IPP
{

SignificantChangePrediction::SignificantChangePrediction(const double &porosThresh)
    : m_changePredictData()
    , m_porosThresh(porosThresh)
{

}

bool SignificantChangePrediction::isSignificant(const std::vector<double>& currPhases,
                                                const std::vector<double>& currTotals,
                                                const std::vector<double>& predictedTotalsMin,
                                                const std::vector<double>& predictedTotalsMax) const

{
    assert(currTotals.size() == predictedTotalsMin.size());
    assert(predictedTotalsMax.size() == predictedTotalsMin.size());


    // check for valid range
    if(m_changePredictData.contains(currTotals) == false
            || m_changePredictData.contains(predictedTotalsMin) == false
            || m_changePredictData.contains(predictedTotalsMax) == false)
    {
        // prediction is not reliable in this domain. fall back to signifcantly change
        return true;
    }



    Interpolation::Element oldPhases;
    m_changePredictData.evaluate(currTotals, oldPhases);


    assert(currPhases.size() == oldPhases.size());
    for(size_t iPhase = 0; iPhase < currPhases.size(); ++iPhase)
    {
        if(
                (currPhases[iPhase] == 0.0 && oldPhases[iPhase] != 0.0)
                || (currPhases[iPhase] != 0.0 && oldPhases[iPhase] == 0.0)
            )
        {
            // prediction is not reliable in this domain. fall back to signifcantly change
            return true;
        }
    }


    if(checkDiff(oldPhases, predictedTotalsMin))
    {
        return true;
    }
    else if(checkDiff(oldPhases, predictedTotalsMax))
    {
        return true;
    }

    return false;
}

SignificantChangePrediction::Interpolation &SignificantChangePrediction::getLookup()
{
    return m_changePredictData;
}

const SignificantChangePrediction::Interpolation::Coord &SignificantChangePrediction::getAccuracy() const
{
    return m_changePredictData.getSpacing();
}

bool SignificantChangePrediction::checkDiff(const Interpolation::Element& oldPhases,
                                            const std::vector<double>& predictedTotals) const
{
    Interpolation::Element newPhases;
    m_changePredictData.evaluate(predictedTotals, newPhases);

    assert(oldPhases.size() == newPhases.size());

    for(size_t iPhase = 0; iPhase < oldPhases.size(); ++iPhase)
    {
        const double& oldVal = oldPhases[iPhase];
        const double& newVal = newPhases[iPhase];
        if( (oldVal == 0.0 && newVal > 0.0) || (oldVal > 0.0 && newVal == 0.0) )
        {
            return true;
        }
    }

    return false;

}

}
