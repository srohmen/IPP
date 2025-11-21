#include "calceffectivediffcoef.h"

#include <assert.h>
#include <cstddef>
#include <algorithm>
#include <functional>


static const double s_maxD = 5.0E-8; // * DiffusionCoeffConstants::s_diffCoeffConvFac;

template<typename T>
class ClampMin
{
public:
    ClampMin(const T& minVal)
        : minVal(minVal)
    {

    }

    void operator()(T& curr)
    {
        curr = std::max(minVal, curr);
    }

private:
    const T& minVal;
};


CalcEffectiveDiffCoef::CalcEffectiveDiffCoef()
    : m_nominatorSum(0.0)
{

}

CalcEffectiveDiffCoef::CalcEffectiveDiffCoef(const std::vector<double> &zVec, const std::vector<double> &DVec)
    : m_zVec(zVec),
      m_DVec(DVec),
      m_nominatorSum(0.0)
{
    this->init(zVec, DVec);
}

void CalcEffectiveDiffCoef::init(const std::vector<double>& zVec, const std::vector<double>& DVec)
{
    m_zVec = zVec;
    m_DVec = DVec;
    m_cSrcVec.resize(zVec.size(), 0.0);
    m_cDstVec.resize(zVec.size(), 0.0);
}

std::vector<double>&CalcEffectiveDiffCoef::getC0()
{
    return m_cSrcVec;
}

std::vector<double>&CalcEffectiveDiffCoef::getC1()
{
    return m_cDstVec;
}

double CalcEffectiveDiffCoef::calc(const std::size_t iSpecies) const
{
    const double& D = m_DVec[iSpecies];
    const double& z = m_zVec[iSpecies];

    if(z == 0)
    {
        return D;
    }
    else
    {

        const double& t = m_tVec[iSpecies];
        const double& cDiff = m_cDiffVec[iSpecies];

        const double denom = z * cDiff;

        // clamp at maximum diffCoeff

        if(denom != 0.0)
        {
            const double div = t * m_nominatorSum / denom;
            const double newD = D - div;

            if(newD <= 0.0)
            {
                return 0.0;
            }

            return std::min(newD, s_maxD);
        }
        else
        {
            return s_maxD;
        }
    }
}

void CalcEffectiveDiffCoef::calc(std::vector<double>& D_effVec) const
{
    const size_t nSpecies = m_zVec.size();
    D_effVec.resize(nSpecies);
    for(size_t i = 0; i < nSpecies; ++i)
    {
        D_effVec[i] = this->calc(i);
    }
}

void CalcEffectiveDiffCoef::init()
{
    // remove invalid zero/negative concentrations
    std::for_each(m_cSrcVec.begin(), m_cSrcVec.end(), ClampMin<double>(1.0E-10));


    const size_t nSpecies = m_zVec.size();

    m_cDiffVec.resize(nSpecies);
    for(size_t j = 0; j < nSpecies; ++j)
    {
        const double& cSrc = m_cSrcVec[j];
        const double& cDst = m_cDstVec[j];
        const double cDiff = cDst - cSrc;
        m_cDiffVec[j] = cDiff;
    }

    const double F = 9.65E+7;
    const double F2 = F*F;
    const double R = 8.314;
    const double T = 300; // TODO: set temp + alter diff coeffs
    const double fac = F2 / (R*T);
    std::vector<double> kVec(nSpecies);
    double kSum = 0.0;
    for(size_t j = 0; j < nSpecies; ++j)
    {
        const double& cSrc = m_cSrcVec[j];
        const double& zj = m_zVec[j];
        const double& Dj = m_DVec[j];

        const double kj = fac * cSrc * zj*zj * Dj;
        kVec[j] = kj;

        kSum += kj;
    }

    m_tVec.resize(nSpecies);
    for(size_t j = 0; j < nSpecies; ++j)
    {
        m_tVec[j] = kVec[j] / kSum;
    }


    m_nominatorSum = 0.0;
    for(size_t j = 0; j < nSpecies; ++j)
    {
        const double& zj = m_zVec[j];
        const double& Dj = m_DVec[j];
        const double& cDiff = m_cDiffVec[j];
        m_nominatorSum += zj * Dj * cDiff;
    }

}
