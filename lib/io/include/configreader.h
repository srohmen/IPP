#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <string>
#include <boost/program_options.hpp>
#include "ippconfig.h"

namespace IPP
{

class BCGenerator;
class SaturationIndexCalcFactory;
class PorosityModelFactory;
class MultiScaleDiffusionFactory;

namespace ConfigReader
{

struct Factories
{
    Factories()
        : bcGen()
        , si_calcFactory()
        , porosModelFactory()
        , multiScaleFactory()
    {

    }

    std::shared_ptr<BCGenerator> bcGen;
    std::shared_ptr<SaturationIndexCalcFactory> si_calcFactory;
    std::shared_ptr<PorosityModelFactory> porosModelFactory;
    std::shared_ptr<MultiScaleDiffusionFactory> multiScaleFactory;
};

IPPConfigPtr read(int argc, char* argv[], const boost::program_options::variables_map &vm,
                  const std::string& fileName, const Factories &factories);

IPPConfigPtr read(int argc, char* argv[], const boost::program_options::variables_map& vm,
                  const std::string& fileName);
}

} // end of namespace IPP

#endif // CONFIGREADER_H
