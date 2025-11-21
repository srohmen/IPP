#include "porosityhelper.h"

#include <chrono>

#include "porositycalcutils.h"
#include "geometrytools.h"
#include "reactionexchangedata.h"
#include "simulationexchangedata.h"
#include "mpitools.h"
#include "phreeqccalcdiffcoefs.h"
#include "ippstream.h"
#include "fielddecomposition.h"
#include "bench_tools.h"


namespace IPP
{

void PorosityHelper::calc(const ReactionExchangeData& reacData,
                          CalcPorosityOutput& output, const bool verbose)
{
    const auto calcPorosTimer = std::chrono::steady_clock::now();

    const SimulationExchangeData& simData = reacData.simData;

    const size_t nDomains = simData.getDecomp()->getLocalNumberOfCells();
    const std::vector<std::string>& phaseNames = simData.getPhaseNames();
    const size_t nPhases = phaseNames.size();

    const std::vector<double>& inertVolFrac = simData.getInertVolFrac();

    const PorosityHelper::CalcPorosityInput input =
    {
        nDomains,
        nPhases,
        reacData.invTotalCellVolume,
        reacData.phaseMolarVolumes,
        reacData.isPermeablePhase,
        inertVolFrac,
        simData.getPrecipField()
    };


    PorosityHelper::calcPorosity(input, output);


    // PhreeqcRMHelper::printNegativePorosWarnings(porosityPhysical);

    if(verbose)
    {
        const auto durationTrans = CalcDuration::calc(calcPorosTimer);
        pcout << "\t-> calc porosity took:\t" << durationTrans.count()
              << " ms" << std::endl;
    }

}

void PorosityHelper::calcPorosity(const CalcPorosityInput& input,
                                  CalcPorosityOutput& output)
{

    output.solidFractions.resize(input.nDomains * input.nPhases);

    // TODO: prevent copy of array. we need that in the current state due to transposing
    output.localPrecipMol = input.precipField;
    GeometryTools::transpose(output.localPrecipMol, input.nDomains, input.nPhases);


    const PorosityCalc porosCalc(input.nDomains,
                                 input.phaseMolarVolumes,
                                 input.isPermeablePhase,
                                 input.inertVolFrac);


    const av::bounds<2> bounds = { (ptrdiff_t)input.nDomains, (ptrdiff_t)input.nPhases};
    PorosityCalc::ConstCellPhasesView precipMolView(output.localPrecipMol, bounds);
    PorosityCalc::CellPhasesView solidFracView(output.solidFractions, bounds);
    porosCalc.calcPorosity(precipMolView, solidFracView, output.porosityTotal);
}





}
