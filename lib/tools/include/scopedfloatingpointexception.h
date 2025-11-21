#ifndef SCOPEDFLOATINGPOINTEXCEPTION_H
#define SCOPEDFLOATINGPOINTEXCEPTION_H

#include <iostream>
#include <cfenv>

namespace IPP
{

class ScopedDisableFloatingPointException
{
public:
    ScopedDisableFloatingPointException(const ScopedDisableFloatingPointException &) = delete;
    ScopedDisableFloatingPointException(const ScopedDisableFloatingPointException &&) = delete;
    ScopedDisableFloatingPointException& operator=(const ScopedDisableFloatingPointException &) = delete;



    ScopedDisableFloatingPointException()
    #ifdef ENABLE_FPE
        :  oldFlags(fedisableexcept(0xFFFFFFFF))
    #endif
    {

    }



#ifdef ENABLE_FPE
    ~ScopedDisableFloatingPointException()
    {
        feenableexcept(oldFlags);
    }
#endif

private:
#ifdef ENABLE_FPE
    fexcept_t oldFlags;
#endif
};


}

#endif // SCOPEDFLOATINGPOINTEXCEPTION_H
