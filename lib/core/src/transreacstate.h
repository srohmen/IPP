#ifndef TRANSREACSTATE_H
#define TRANSREACSTATE_H

#include <memory>

#include "ipprendererfwd.h"

namespace IPP
{

class TransportModule;
class ReactionModule;
class SimulationExchangeData;

struct TransReacState
{
    TransReacState();
    virtual ~TransReacState();

    TransportModule* transModule;
    ReactionModule* reacModule;
    std::shared_ptr<SimulationExchangeData> data;
    IPPRendererPtr renderer;

    double fluxConvergenceTol;
    bool isConverged;
};

}

#endif // TRANSREACSTATE_H
