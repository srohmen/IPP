#ifndef DEBUGDATA_H
#define DEBUGDATA_H

#include <vector>

namespace IPP
{

#ifdef IPP_DEBUG

struct DebugData
{    
#ifdef IPP_DEBUG_ENABLED
    using EnabledField = std::vector<char>;
    using EnabledFieldVec = std::vector<EnabledField>;

    EnabledFieldVec enabledCellsSteps;
#endif


};

#endif

}

#endif // DEBUGDATA_H
