#ifndef CREATEFILENAME_H
#define CREATEFILENAME_H

#include <string>

namespace CreateFileName
{
    std::string create(const size_t number);
    std::string create(const std::string& prefix, const size_t number);
}

#endif // CREATEFILENAME_H
