#ifndef IPPCOUPLING_H
#define IPPCOUPLING_H


namespace IPP
{

class TransportModule;
class ReactionModule;

namespace IPPCoupling
{
    void runReactionModule(TransportModule& trans, ReactionModule &reac, const bool verbose);    
}

}


#endif // IPPCOUPLING_H
