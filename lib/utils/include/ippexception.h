#ifndef IPPEXCEPTION_H
#define IPPEXCEPTION_H

#include <string>
#include <stdexcept>
#include <iostream>

namespace IPP
{

class IPPException : public std::runtime_error
{
public:
    IPPException(const std::string& arg)
        : std::runtime_error(arg)
        , iDomain(-1)
    {

    }

    size_t iDomain;
};

namespace IPPCheck
{

inline void assertCheck(const bool condition, const std::string& errStr = "unknown error")
{
    if(condition == false)
    {
        throw std::runtime_error(errStr);
    }
}

inline void warnCheck(const bool condition, const std::string& errStr = "unknown warning")
{
    if(condition == false)
    {
        std::cerr << "WARNING: " << errStr << std::endl;
    }
}

}

} // end of namespace IPP

#endif // IPPEXCEPTION_H

