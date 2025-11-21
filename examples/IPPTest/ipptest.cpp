#include <assert.h>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>
#include <ipp.h>
#include <configreader.h>
#include <jsonboundaryconditiongen.h>
#include <jsonsaturationindexcalcfactory.h>
#include <jsonporositymodelfactory.h>
#include <jsonmultiscalediffusionfactory.h>
#include <vtkrenderer.h>

#ifdef IPP_ENABLE_VALGRIND
#include <valgrind/callgrind.h>
#endif

#include <cfenv>

#include "totalbudgettest.h"

#ifdef IPP_ENABLE_CORE_DUMP
#ifdef __linux__
    #include <sys/resource.h>
#endif
#endif

int main(int argc, char *argv[])
{
#ifdef ENABLE_FPE
    feenableexcept(FE_DIVBYZERO | FE_OVERFLOW);
#endif

#ifdef IPP_ENABLE_CORE_DUMP

#ifdef __linux__
    struct rlimit core_limits;
    core_limits.rlim_cur = core_limits.rlim_max = RLIM_INFINITY;
    setrlimit(RLIMIT_CORE, &core_limits);
#endif

#endif

    namespace po = boost::program_options;

    po::options_description desc("Allowed options");

    // clang-format off
    desc.add_options()
        ("help,h", "produce help message")
        ("verbose,v","silence a lot of status messages")
        ("input,i", po::value<std::string>(), "input json file")
        ("checkpoint,c", po::value<size_t>(), "starting checkpoint iteration")
        ("output,o", po::value<std::string>(), "output root")
        ("phases,p", po::value<std::string>(), "phases vti file")
        ("molarity,m", po::value<std::string>(), "components concentration vti file");
    // clang-format on

    po::variables_map vm;

    po::store(boost::program_options::command_line_parser(argc, argv)
                  .options(desc)
                  .style(boost::program_options::command_line_style::unix_style
                         | boost::program_options::command_line_style::allow_long_disguise)
                  .run(),
              vm);

    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    std::string inputFileName;
    if (vm.count("input")) {
        inputFileName = vm["input"].as<std::string>();
    } else {
        std::cout << "missing input file. usage: " << argv[0] << " -i <input.json>" << std::endl;
        return 1;
    }

    IPP::ConfigReader::Factories factories;
    factories.bcGen.reset(new IPP::JSONBoundaryConditionGen);
    factories.si_calcFactory.reset(new IPP::JSONSaturationIndexCalcFactory);
    factories.porosModelFactory.reset(new IPP::JSONPorosityModelFactory);
    factories.multiScaleFactory.reset(new IPP::JSONMultiScaleDiffusionFactory);

    IPP::IPPConfigPtr conf = IPP::ConfigReader::read(argc, argv, vm, inputFileName, factories);
    assert(conf.get());

    if (vm.count("output")) {
        const std::string outputFolder = vm["output"].as<std::string>();
        conf->outputDir = outputFolder;
    }

    IPP::IPP ipp;
    ipp.init(conf);

    if (vm.count("checkpoint")) {
        // std::cout << "using checkpoint!" << std::endl;
        const size_t iteration = vm["checkpoint"].as<size_t>();
        // std::cout << "loading checkpoint at iteration: " << iteration << std::endl;
        ipp.loadCheckpoint(iteration);
    }

#ifdef IPP_ENABLE_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_DUMP_STATS;
#endif

    ipp.execute();

#ifdef IPP_ENABLE_CALLGRIND
    CALLGRIND_STOP_INSTRUMENTATION;
    CALLGRIND_DUMP_STATS;
#endif

    return 0;
}
