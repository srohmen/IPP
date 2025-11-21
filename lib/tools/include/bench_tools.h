#ifndef BENCH_TOOLS_H
#define BENCH_TOOLS_H

#include <chrono>
#include <sstream>

#ifdef IPP_ENABLE_BENCH_BARRIER
#include "mpimanager.h"
#include <mpi.h>
#endif

#define BENCH_FORCE(func) \
    { \
    typedef std::chrono::microseconds TimeT; \
    auto start = std::chrono::steady_clock::now(); \
    func; \
    auto duration = std::chrono::duration_cast<TimeT> \
            (std::chrono::steady_clock::now() - start); \
    std::stringstream ss;   \
    ss << #func << "\t" << duration.count() << " us" << std::endl; \
    std::cout << ss.str(); \
    }

#ifdef IPP_ENABLE_BENCHMARK
    #define BENCH BENCH_FORCE
#else
    #define BENCH(func) \
        func;
#endif

namespace IPP
{

namespace CalcDuration
{

template<typename T>
auto calc(const T& start)
{
#ifdef IPP_ENABLE_BENCH_BARRIER
    MPI_Barrier(MPIManager::getInstance().getCommunicator());
#endif

    auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>
            (std::chrono::steady_clock::now() - start );
    return duration;
}

}

}

#endif // BENCH_TOOLS_H

