#include "phreeqcreactionmoduledata.h"

#include "cellid_si_func.h"
#include "cellid_dissolveonly_func.h"

#include "reactionmoduleoptimization.h"

namespace IPP
{


PhreeqcReactionModuleData::PhreeqcReactionModuleData()
    : siFunc(nullptr)
    , dissolveOnlyFunc(nullptr)
    , optim(nullptr)
{

}

PhreeqcReactionModuleData::~PhreeqcReactionModuleData()
{
    delete siFunc;
    siFunc = nullptr;

    delete dissolveOnlyFunc;
    dissolveOnlyFunc = nullptr;

    delete optim;
    optim = nullptr;

}

}
