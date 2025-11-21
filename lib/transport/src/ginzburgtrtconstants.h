#ifndef GINZBURGTRTCONSTANTS_H
#define GINZBURGTRTCONSTANTS_H

#include <cstddef>

namespace GinzburgTRTConstants
{

inline double calc_cPhi(const double& minPoros, const bool is3D)
{
    const double cFact = 1.0 / 3.0; // is3D ? 1.0 / 3.0 : 1.0 / 3.0;
    return minPoros * cFact;
}

}

#endif // GINZBURGTRTCONSTANTS_H
