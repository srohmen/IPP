#include "createfilename.h"

#include <iomanip>

namespace CreateFileName
{

std::string create(const size_t number)
{
    // TODO: check if this is enough or make it configurable
    const size_t width = 9;

    std::stringstream fNameStream;
    fNameStream << std::setfill('0') << std::setw(width) << number;
    return fNameStream.str();
}

std::string create(const std::string& prefix, const size_t number)
{
    std::stringstream fNameStream;
    fNameStream << prefix << create(number);
    return fNameStream.str();
}


}
