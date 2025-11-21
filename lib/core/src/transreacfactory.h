#ifndef TRANSREACFACTORY_H
#define TRANSREACFACTORY_H

namespace IPP
{

struct TransReacState;
class IPPConfig;

namespace TransReacFactory
{
    TransReacState* create(const IPPConfig& conf);
}

}

#endif // TRANSREACFACTORY_H
