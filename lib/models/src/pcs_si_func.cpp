#include "pcs_si_func.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include "ippstream.h"

namespace IPP
{

template<typename PCSFunc>
PCS_SI_Func<PCSFunc>::PCS_SI_Func(const double& temp, const double& porosInit, const double& rInit,
                                  const double& surfTensFallbackConc,
                                  const std::map<std::string, double> &surfaceTensions)
    : temp(temp)
    , porosInit(porosInit)
    , rInit(rInit)
    , surfTensFallbackConc(surfTensFallbackConc)
    , surfaceTensions(surfaceTensions)
{

}

static double calcSurfTens(const double& temp, const double Vm, const double& conc)
{
    const double beta = 0.514;
    const double R = 8.314;
    const double NA = 6.022E+23;
    const double kb = R / NA;
    const double kT = temp * kb;
    const double beta_kT = beta * kT;
    const double moleculeVol = Vm / NA;
    const double surfTens = beta_kT / std::pow(moleculeVol, 2.0 / 3.0) * std::log( 1.0 / (conc * Vm) );

    return surfTens;
}

template<typename PCSFunc>
void PCS_SI_Func<PCSFunc>::init(const std::vector<std::string> &phaseNames,
                                const PhaseNameToInfos &phaseInfos,
                                const std::vector<double>* poros)
{
    for(size_t iPhase = 0; iPhase < phaseNames.size(); ++iPhase)
    {
        const std::string& name = phaseNames[iPhase];

        // convert to m3/mol
        const double& Vm = phaseInfos.at(name).molarVolume * 1.0E-6;

        double surfTens;
        if(surfaceTensions.find(name) != surfaceTensions.end())
        {
            surfTens = surfaceTensions.at(name);
        }
        else
        {
            surfTens = calcSurfTens(temp, Vm, surfTensFallbackConc);
        }

        pcout << "surface tensions: " << name << "\t" << surfTens << std::endl;

        m_pcsFuncs.emplace_back(temp, surfTens, Vm, porosInit, rInit);
    }

    m_poros = poros;
}

template<typename PCSFunc>
void PCS_SI_Func<PCSFunc>::evaluate(const size_t iCell, const iterator &begin, const iterator &end) const
{
    const double& poros = (*m_poros)[iCell];
    size_t iPhase = 0;
    for(iterator it = begin; it != end; ++it, ++iPhase)
    {
        assert(m_pcsFuncs.size() > iPhase);
        const PCSFunc& func = m_pcsFuncs[iPhase];
        const double SI = func.calcSI(poros);
        *it = SI;
    }
}

template class PCS_SI_Func<PorosityControlledSolubility::Cylinder>;

}
