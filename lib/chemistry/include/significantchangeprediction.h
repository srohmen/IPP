#ifndef SIGNIFICANTCHANGEPREDICTION_H
#define SIGNIFICANTCHANGEPREDICTION_H

#include "multidimnearestneighborinterpolation.h"

namespace IPP
{

class SignificantChangePrediction
{
public:
    SignificantChangePrediction(const double& porosThresh);


    bool isSignificant(const std::vector<double> &currPhases,
                       const std::vector<double>& currTotals,
                       const std::vector<double>& predictedTotalsMin,
                       const std::vector<double> &predictedTotalsMax) const;


    typedef MultiDimNearestNeighborInterpolation<double> Interpolation;
    Interpolation& getLookup();
    const Interpolation::Coord& getAccuracy() const;

private:
    bool checkDiff(const Interpolation::Element& oldPhases,
                   const std::vector<double>& predictedTotals) const;


    Interpolation m_changePredictData;

    const double m_porosThresh;

};

}

#endif // SIGNIFICANTCHANGEPREDICTION_H
