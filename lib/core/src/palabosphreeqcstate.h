#ifndef PALABOSPHREEQCSTATE_H
#define PALABOSPHREEQCSTATE_H

#include "transreacstate.h"

namespace IPP
{

class PhreeqcReactionModuleData;

struct PalabosPhreeqcState : public TransReacState
{
    PalabosPhreeqcState();
    virtual ~PalabosPhreeqcState();

    std::unique_ptr<PhreeqcReactionModuleData> phreeqcData;
};


}


#endif // PALABOSPHREEQCSTATE_H
