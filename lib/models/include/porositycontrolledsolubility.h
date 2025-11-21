#ifndef POROSITYCONTROLLEDSOLUBILITY_H
#define POROSITYCONTROLLEDSOLUBILITY_H

namespace IPP
{

namespace PorosityControlledSolubility
{

class Cylinder
{
public:
    Cylinder(const double& temp, const double& surfTens,
             const double& molarVol,
             const double &porosInit, const double& rInit);
    double calcSI(const double& porosity) const;

private:
    const double m_geomFactor;
    const double m_k;
};

}

}

#endif // POROSITYCONTROLLEDSOLUBILITY_H
