#include "nucleationpropcalc.h"
#include <cmath>

namespace IPP
{

static double calcB(const double& kT, const double& V_molecule,
                    const double& interfacialTensions)
{

    const double beta = 16.0 * M_PI / 3.0;
    const double B = (beta * V_molecule * V_molecule *
                     interfacialTensions * interfacialTensions * interfacialTensions)
                     / (kT * kT);
    return B;
}

static const double kb = 1.38060445034872E-23;

NucleationPropCalc::NucleationPropCalc(const double& dt,
                                       const double& T,
                                       const double& D,
                                       const double& V_molecule,
                                       const double& N0,
                                       const double& interfacialTensions)
    : m_dt(dt)
    , m_kT(kb * T)
    , m_N0(N0)
    , m_B(calcB(m_kT, V_molecule, interfacialTensions))
    , m_c(D * std::sqrt(m_kT / interfacialTensions))
{

}

double NucleationPropCalc::calcLnSurvivalProp(const double &lnS,
                                              const double &N1,
                                              const double &volSurf) const
{
    const double nuclRate = calcNucleationRate(lnS, N1);
    const double ln_JVt = -nuclRate * volSurf * m_dt;
    return ln_JVt;
}

double NucleationPropCalc:: calcNucleationProp(const double& lnS,
                                              const double& N1,
                                              const double& volSurf) const
{
    const double nuclRate = calcNucleationRate(lnS, N1);
    const double ln_JVt = -nuclRate * volSurf * m_dt;
    const double nuclProp = 1.0 - std::exp(ln_JVt);
    return nuclProp;
}

double NucleationPropCalc::calcNucleationRate(const double& lnS, const double& N1) const
{
    const double J0 = m_N0 * N1 * m_c * lnS;

    const double dG = m_B / (lnS * lnS);
    const double J = J0 * std::exp(-dG / m_kT);

    return J;
}

}
