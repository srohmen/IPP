#include "palabosphreeqcstate.h"

#include "phreeqcreactionmoduledata.h"

namespace IPP
{

PalabosPhreeqcState::PalabosPhreeqcState()
    : phreeqcData(new PhreeqcReactionModuleData)
{

}

PalabosPhreeqcState::~PalabosPhreeqcState()
{

}

}
