#include <iostream>

#include <chrono>

#include "tennisjenningscshporositycalc.h"
#include "degradedcshdiffusioncalc.h"
#include "cshgeldiffusioncoefcalc.h"
#include "moritanakadiffusioncalc.h"

#include "array_view.h"
#include "arrayviewiterator.h"
#include "arrayviewiteratorimpl.h"

#include "cshporosityinfos.h"

#include <fenv.h>

using namespace IPP;


static void calcPortSeries(const PhaseNameToInfos& phaseInfos,
                const TennisJenningsCSHPorosityCalc& porosCalc,
                const DegradedCSHDiffusionCalc& diffCalc)
{
    const std::ptrdiff_t nSolids = phaseInfos.size();
    const av::bounds<2> amountBounds = {1, nSolids};


    static const double intrGelPorLD = 0.373;
    // static const double intrGelPorHD = 0.237;

    const double refVol = 999.999999;

    const size_t iMax = 20;
    for(size_t i = 0; i <= iMax; ++i)
    {
        const double vFracCSHGel = (double)i / iMax;
        const double vFracCSHDry = vFracCSHGel * (1.0 - intrGelPorLD);
        const double vCSH = refVol * vFracCSHDry; // cm3
        const double nCSH = vCSH / phaseInfos.at("CSH1.6").molarVolume;

        const double vFracPort = 1.0 - vFracCSHGel;
        const double vPort = refVol * vFracPort;
        const double nPort = vPort / phaseInfos.at("Portlandite").molarVolume;

        std::vector<double> phaseAmounts(nSolids, 0.0);
        phaseAmounts[0] = 0;
        phaseAmounts[1] = 0;
        phaseAmounts[2] = nCSH;
        phaseAmounts[3] = nPort;


        av::array_view<const double, 2> amountView(phaseAmounts, amountBounds);
        av::strided_array_view<const double, 2> view = amountView.section({0, 0}, {1, nSolids});

        AbstractPorosityCalc::PorosityCalcInput input;
        input.phaseAmounts.begin = ArrayViewIterator(new ArrayViewIteratorImpl(view, std::begin(view.bounds())));
        input.phaseAmounts.end = ArrayViewIterator(new ArrayViewIteratorImpl(view,  std::end(view.bounds())));


        const AbstractPorosityCalc::SimplePorosityInfosPtr porosInfos = porosCalc.calc(input);
        const double diffusionCoef = diffCalc.calc(*porosInfos);

        const CSHPorosityInfos& cshPorosInfos = static_cast<const CSHPorosityInfos&>(*porosInfos);

        std::cout << vFracPort
                  << "\t" << cshPorosInfos.porosityTotal
                  << "\t" << cshPorosInfos.porosityCapillary
                  << "\t" << diffusionCoef
                  << std::endl;

    }

}


static void calcCSHSeries(const PhaseNameToInfos& phaseInfos,
                const TennisJenningsCSHPorosityCalc& porosCalc,
                const DegradedCSHDiffusionCalc& diffCalc)
{
    const std::ptrdiff_t nSolids = phaseInfos.size();
    const av::bounds<2> amountBounds = {1, nSolids};


    static const double intrGelPorLD = 0.373;
    // static const double intrGelPorHD = 0.237;

    const size_t iMax = 100;
    for(size_t i = 0; i <= iMax; ++i)
    {
        const double vFracCSHGel = (double)i / iMax;
        const double vFracCSHDry = vFracCSHGel * (1.0 - intrGelPorLD);
        const double vCSH = 1000 * vFracCSHDry; // cm3
        const double nCSH = vCSH / phaseInfos.at("CSH1.6").molarVolume;

//        const double vFracCH = vMax - vFracCSH;
//        const double nCH = vFracCH / phaseInfos.at("Portlandite").molarVolume;

        std::vector<double> phaseAmounts(nSolids, 0.0);
        phaseAmounts[0] = nCSH;
        phaseAmounts[1] = 0;
        phaseAmounts[2] = 0;
        phaseAmounts[3] = 0.0;


        av::array_view<const double, 2> amountView(phaseAmounts, amountBounds);
        av::strided_array_view<const double, 2> view = amountView.section({0, 0}, {1, nSolids});

        AbstractPorosityCalc::PorosityCalcInput input;
        input.phaseAmounts.begin = ArrayViewIterator(new ArrayViewIteratorImpl(view, std::begin(view.bounds())));
        input.phaseAmounts.end = ArrayViewIterator(new ArrayViewIteratorImpl(view,  std::end(view.bounds())));


        const AbstractPorosityCalc::SimplePorosityInfosPtr porosInfos = porosCalc.calc(input);
        const double diffusionCoef = diffCalc.calc(*porosInfos);

        const CSHPorosityInfos& cshPorosInfos = static_cast<const CSHPorosityInfos&>(*porosInfos);

        std::cout << vFracCSHGel
                  << "\t" << cshPorosInfos.porosityTotal
                  << "\t" << cshPorosInfos.porosityCapillary
                  << "\t" << diffusionCoef
                  << std::endl;

    }

}

static void calcSingle(const PhaseNameToInfos& phaseInfos,
                const TennisJenningsCSHPorosityCalc& porosCalc,
                const DegradedCSHDiffusionCalc& diffCalc)
{
    const std::ptrdiff_t nSolids = phaseInfos.size();

    std::vector<double> phaseAmounts(nSolids, 0.0);
    phaseAmounts[0] = 0.0;
    phaseAmounts[1] = 0.0;
    phaseAmounts[2] = 627 / phaseInfos.at("CSH1.6").molarVolume;
    phaseAmounts[3] = 0.0 / phaseInfos.at("Portlandite").molarVolume;


    const av::bounds<2> amountBounds = {1, nSolids};
    av::array_view<const double, 2> amountView(phaseAmounts, amountBounds);

    av::strided_array_view<const double, 2> view = amountView.section({0, 0}, {1, nSolids});

    AbstractPorosityCalc::PorosityCalcInput input;
    input.phaseAmounts.begin = ArrayViewIterator(new ArrayViewIteratorImpl(view, std::begin(view.bounds())));
    input.phaseAmounts.end = ArrayViewIterator(new ArrayViewIteratorImpl(view,  std::end(view.bounds())));

    auto start = std::chrono::steady_clock::now();

    const AbstractPorosityCalc::SimplePorosityInfosPtr porosInfos = porosCalc.calc(input);
    const double diffusionCoef = diffCalc.calc(*porosInfos);

    const CSHPorosityInfos& cshPorosInfos = static_cast<const CSHPorosityInfos&>(*porosInfos);

    std::cout << "massFracLD:\t" << porosCalc.getMassFractionCSH_LD() << std::endl;
    std::cout << "porosityTotal:\t"  << cshPorosInfos.porosityTotal << std::endl;
    std::cout << "porosityCapillary:\t" << cshPorosInfos.porosityCapillary << std::endl;
    std::cout << "volNonPerm:\t" << cshPorosInfos.volNonPerm << std::endl;
    std::cout << "volSatCSH:\t" << cshPorosInfos.volSatCSH << std::endl;
    std::cout << "volFracCSH_HD:\t" << cshPorosInfos.volFracCSH_HD << std::endl;
    std::cout << "ratioCaSi:\t" << cshPorosInfos.ratioCaSi << std::endl;
    std::cout << "diffCoef:\t" << diffusionCoef << std::endl;


    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "\t-> calc took:\t" << duration.count() << " ms" << std::endl;
}

int main(/*int argc, char* argv[]*/)
{
    feenableexcept(
                FE_INVALID
                | FE_DIVBYZERO
                | FE_OVERFLOW
                | FE_UNDERFLOW
                );

    const double D0 = 2.2E-9;
    const double wc = 0.607661846030389;
    const double degreeHydration = 0.95;

    TennisJenningsCSHPorosityCalc porosCalc(wc, degreeHydration);
    DegradedCSHDiffusionCalc diffCalc(new CSHGelDiffusionCoefCalc,
                                  new MoriTanakaDiffusionCalc,
                                  D0);

    std::cout << "frac LD:\t" << porosCalc.getMassFractionCSH_LD() << std::endl;


    std::vector<std::string> compNames;
    compNames.push_back("H2O");
    compNames.push_back("H");
    compNames.push_back("O");
    compNames.push_back("Charge");
    compNames.push_back("Ca");
    compNames.push_back("Si");


    PhaseNameToInfos phaseInfos;
    phaseInfos["CSH0.8"].gfw = 132.688;
    phaseInfos["CSH0.8"].molarVolume = 59.290;
    phaseInfos["CSH0.8"].stoich.resize(compNames.size(), 0.0);
    phaseInfos["CSH0.8"].stoich[4] = 0.8;
    phaseInfos["CSH0.8"].stoich[5] = 1.0;

    phaseInfos["CSH1.2"].gfw = 164.486;
    phaseInfos["CSH1.2"].molarVolume = 71.95;
    phaseInfos["CSH1.2"].stoich.resize(compNames.size(), 0.0);
    phaseInfos["CSH1.2"].stoich[4] = 1.2;
    phaseInfos["CSH1.2"].stoich[5] = 1.0;

    phaseInfos["CSH1.6"].gfw = 196.285;
    phaseInfos["CSH1.6"].molarVolume = 84.68;
    phaseInfos["CSH1.6"].stoich.resize(compNames.size(), 0.0);
    phaseInfos["CSH1.6"].stoich[4] = 1.6667;
    phaseInfos["CSH1.6"].stoich[5] = 1.0;

    phaseInfos["Portlandite"].gfw = 74.0918;
    phaseInfos["Portlandite"].molarVolume = 32.9;
    phaseInfos["Portlandite"].stoich.resize(compNames.size(), 0.0);
    phaseInfos["Portlandite"].stoich[4] = 1.0;



    std::vector<std::string> phaseNames;
    phaseNames.push_back("CSH0.8");
    phaseNames.push_back("CSH1.2");
    phaseNames.push_back("CSH1.6");
    phaseNames.push_back("Portlandite");

    std::vector<std::string> permPhases;
    permPhases.push_back("CSH0.8");
    permPhases.push_back("CSH1.2");
    permPhases.push_back("CSH1.6");



    porosCalc.setCompNames(&compNames);
    porosCalc.setPhaseInfos(&phaseInfos);
    porosCalc.setPhaseNames(&phaseNames);

    for(const std::string& str : permPhases)
    {
        porosCalc.add_CSH_Phase(str);
    }

    calcSingle(phaseInfos, porosCalc, diffCalc);
    calcCSHSeries(phaseInfos, porosCalc, diffCalc);
    calcPortSeries(phaseInfos, porosCalc, diffCalc);


    return 0;
}
