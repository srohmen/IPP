#ifndef CALCEFFECTIVEDIFFCOEF_H
#define CALCEFFECTIVEDIFFCOEF_H

#include <vector>

class CalcEffectiveDiffCoef
{
public:
    CalcEffectiveDiffCoef();
    CalcEffectiveDiffCoef(const std::vector<double>& m_zVec,
                          const std::vector<double>& m_DVec);
    void init(const std::vector<double>& m_zVec,
              const std::vector<double>& m_DVec);

    std::vector<double>& getC0();

    std::vector<double>& getC1();

    double calc(const std::size_t iSpecies) const;

    void calc(std::vector<double>& D_effVec) const;


    void init();


private:
    std::vector<double> m_zVec;
    std::vector<double> m_DVec;

    std::vector<double> m_cSrcVec;
    std::vector<double> m_cDstVec;
    std::vector<double> m_cDiffVec;

    std::vector<double> m_tVec;
    double m_nominatorSum;
};
#endif // CALCEFFECTIVEDIFFCOEF_H
