#include "latticeboltzmannunitconversion.h"

#include "ippconfig.h"
#include "ginzburgtrtconstants.h"

namespace IPP
{

namespace LatticeBoltzmannUnitConversion
{

double calcTimeStepSRT(const double& dx, const double& D, const double& tau)
{
    return (dx*dx / (D * 3)) * (tau - 0.5);
}

double calcTimeStepPTRT(const double& dx, const double& D, const double& tau,
                        const double& porosLow, const double porosRef, const bool is3D)
{
    const double cPhi = GinzburgTRTConstants::calc_cPhi(porosLow, is3D);
    const double D_LB = (cPhi / porosRef) * (tau - 0.5);
    const double D_conv = D / D_LB;
    const double dt = (dx*dx) / D_conv;
    return dt;
}

double calcTimeStep(const IPPConfig& conf)
{
    const double& dx = conf.spatialResolution;
    const double& D = conf.diffusionCoefReference;
    const double& tau = conf.tauRef;

    switch(conf.latticeBoltzmannCollType)
    {
        case(LBCT_DVSRT):
        {
            return calcTimeStepSRT(dx, D, tau);
            break;
        }

        case(LBCT_PTRT):
        {
            return calcTimeStepPTRT(dx, D, tau, conf.porosLow, conf.porosRef, conf.is3DSimulation);
            break;
        }


        default:
        {
            throw std::runtime_error("unknow lattice boltzmann collission type: "
                                     + std::to_string(conf.latticeBoltzmannCollType));
            break;
        }


    }

}


}
}
