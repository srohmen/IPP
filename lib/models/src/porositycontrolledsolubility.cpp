#include "porositycontrolledsolubility.h"
#include <cassert>
#include <cmath>

namespace IPP
{

namespace PorosityControlledSolubility
{

const double R = 8.314;

Cylinder::Cylinder(const double& temp, const double& surfTens,
                   const double& molarVol,
                   const double &porosInit, const double& rInit)
    : m_geomFactor( (rInit * rInit) / porosInit)
    , m_k(surfTens * molarVol / (R*temp))
{
    assert(m_geomFactor >= 0.0);
    assert(m_k >= 0.0);
}

double Cylinder::calcSI(const double &porosity) const
{
    assert(porosity >= 0.0);

    // dense sphere packing
    const double critPoros = 0.74;

    if(porosity > critPoros)
    {
        return 0.0;
    }

    const double r = std::sqrt(porosity * m_geomFactor);

    const double ln_f = m_k / r;
    const double f = std::exp(ln_f);

    const double slope = 1.0 / critPoros;
    const double linear = slope * porosity;
    const double fAdapt = f * (1.0 - linear) + linear;

    // values lower 1 should be excluded by critPoros check above
    assert(fAdapt >= 1.0);

    const double SI_target = std::log10(fAdapt);

    return SI_target;
}


}

}
