#ifndef NUCLEATIONPROPCALC_H
#define NUCLEATIONPROPCALC_H

namespace IPP
{

class NucleationPropCalc
{
public:
    NucleationPropCalc(const double& dt,
                       const double& T, const double& D,
                       const double& V_molecule,
                       const double& N0,
                       const double& interfacialTensions);

    double calcLnSurvivalProp(const double& lnS, const double& N1, const double& volSurf) const;
    double calcNucleationProp(const double& lnS, const double& N1, const double& volSurf) const;

private:
    double calcNucleationRate(const double& lnS, const double& N1) const;

    const double& m_dt;
    const double m_kT;
    const double m_N0;

    const double m_B;
    const double m_c;

};

}

#endif // NUCLEATIONPROPCALC_H
