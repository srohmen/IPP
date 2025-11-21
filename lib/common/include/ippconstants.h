#ifndef IPPCONSTANTS_H
#define IPPCONSTANTS_H

#include <string>

namespace IPPConstants
{

static const std::string s_chemOutput = "chem";
static const std::string s_resultOutput = "results";
static const std::string s_checkpointOutput = "checkpoints";
static const std::string s_dumpExtension = ".dat";

static const double s_isSteadyFlag = 1.0;
static const double s_isNotSteadyFlag = 0.0;

} // end of namespace IPPConstants

#endif // IPPCONSTANTS_H
