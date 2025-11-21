#include "transreacstate.h"

#include "transportmodule.h"
#include "reactionmodule.h"
#include "simulationexchangedata.h"

namespace IPP
{

TransReacState::TransReacState()
    : transModule(nullptr)
    , reacModule(nullptr)
    , data(nullptr)
    , isConverged(false)
{

}

TransReacState::~TransReacState()
{
    delete transModule;
    transModule = nullptr;

    delete reacModule;
    reacModule = nullptr;
}


}
