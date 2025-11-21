#include "ippcoupling.h"

#include "transportmodule.h"
#include "reactionmodule.h"

#include "ippstream.h"
#include "mpitools.h"
#include "bench_tools.h"


namespace IPP
{

static void runCollectData(TransportModule& trans, const bool verbose)
{    
    auto start = std::chrono::steady_clock::now();

    trans.collectData();

    if(verbose)
    {
        auto duration = CalcDuration::calc(start);
        pcout << "\t-> getting post transport concentration from TM took:\t"
              << duration.count() << " ms" << std::endl;
    }

}

static void runPreCondition(ReactionModule& reac, const bool verbose)
{
    auto start = std::chrono::steady_clock::now();

    reac.updatePreCondition();

    if(verbose)
    {
        auto duration = CalcDuration::calc(start);
        pcout << "\t-> updating pre-condition in RM took:\t"
              << duration.count() << " ms" << std::endl;
    }
}


static void runPostCondition(ReactionModule& reac, const bool verbose)
{
    auto start = std::chrono::steady_clock::now();

    reac.updatePostCondition();

    if(verbose)
    {
        auto duration = CalcDuration::calc(start);
        pcout << "\t-> updating post-condition in RM took:\t"
              << duration.count() << " ms" << std::endl;
    }
}


static void runFinish(TransportModule& trans, const bool verbose)
{
    MPI_CHECK_SYNC;

    auto start = std::chrono::steady_clock::now();
    trans.updatePostReactionState();

    if(verbose)
    {
        auto duration = CalcDuration::calc(start);
        pcout << "\t-> setting post reaction state in TM took:\t"
              << duration.count() << " ms" << std::endl;
    }

    MPI_CHECK_SYNC;
}

static void runConcUpdate(ReactionModule& reac, const bool verbose)
{
    MPI_CHECK_SYNC;

    auto start = std::chrono::steady_clock::now();

    reac.updateComponentConcentrations();

    if(verbose)
    {
        auto duration = CalcDuration::calc(start);
        pcout << "\t-> setting concentrations in RM took:\t"
              << duration.count() << " ms" << std::endl;
    }

    MPI_CHECK_SYNC;
}

void IPPCoupling::runReactionModule(TransportModule& trans, ReactionModule& reac, const bool verbose)
{
    MPI_CHECK_SYNC;

    runCollectData(trans, verbose);
    runPreCondition(reac, verbose); // includes porosity update, must be done before conc update
    runConcUpdate(reac, verbose);

    reac.run();

    // must be executed after runcells, because this updates target SI and dissolve only state
    // in loading checkpoint this should not be reset to wrong state
    runPostCondition(reac, verbose);

    runFinish(trans, verbose);

    MPI_CHECK_SYNC;
}


}
