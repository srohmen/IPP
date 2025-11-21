#include "goldensectionsearch.h"

namespace IPP
{

bool GoldenSectionSearch::isConverged(const double& x0, const double& x1, const double& tolerance)
{
    double xm = 0.5 * std::abs( x1 + x0 );

    if ( xm <= 1.0 )
    {
        return ( std::abs( x1 - x0 ) < tolerance ) ? true : false;
    }

    return ( std::abs( x1 - x0 ) < tolerance * xm ) ? true : false;
}

bool GoldenSectionSearch::isConverged(const std::vector<double>& x0Vec,
                                      const std::vector<double>& x1Vec,
                                      const double& tolerance)
{
    for(size_t i = 0; i < x0Vec.size(); ++i)
    {
        const double& x0 = x0Vec[i];
        const double& x1 = x1Vec[i];
        const bool isCurrConverged = isConverged(x0, x1, tolerance);
        if (isCurrConverged == false)
        {
            return false;
        }
    }

    return true;
}

bool GoldenSectionSearch::isConverged(std::vector<size_t>& toCheck,
                                      const std::vector<double>& x0Vec,
                                      const std::vector<double>& x1Vec,
                                      const double& tolerance)
{
    std::vector<size_t>::const_iterator it = toCheck.begin();
    while(it != toCheck.end())
    {
        const size_t iCell = *it;
        const double& x0 = x0Vec[iCell];
        const double& x1 = x1Vec[iCell];
        const bool isCurrConverged = isConverged(x0, x1, tolerance);
        if(isCurrConverged)
        {
            it = toCheck.erase(it);
        }
        else
        {
            ++it;
        }
    }

    return toCheck.empty();
}

} // end of namespace IPP
