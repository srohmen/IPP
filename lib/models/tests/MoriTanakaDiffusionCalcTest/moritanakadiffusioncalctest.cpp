#include <iostream>

#include <chrono>
#include <cmath>

#include "moritanakadiffusioncalc.h"
#include "cshgeldiffusioncoefcalc.h"

using namespace IPP;


int main(/*int argc, char* argv[]*/)
{
    const double D0 = 2.2e-9;

    double poros = 1.0;
    const double f = 0.9;
    for(size_t i = 0; i < 200; ++i)
    {
        MoriTanakaDiffusionCalc::InclusionsData data;
        data.push_back({ 1.0 - poros, 0.0 });

        MoriTanakaDiffusionCalc diffCalc;
        const double diffMT = diffCalc.calcDiffusionCoeffient(D0, data);
        const double diffSpheres = CSHGelDiffusionCoefCalc::calc(D0, poros);

        const double porosThresh = 0.15;
        const double t = std::min(1.0, poros / porosThresh);
        const double intermA = diffMT * t + (1.0 - t) * diffSpheres;

        const double k = 10;
        const double u = 1.0 - std::exp(-poros * k);
        const double intermB = diffMT * u + (1.0 - u) * diffSpheres;

        std::cout << poros << "\t" << u << std::endl;
        std::cout << poros << "\t" << diffMT << "\t" << diffSpheres << "\t" << intermA << "\t"
                  << intermB << std::endl;

        poros *= f;

    }

    return 0;
}
