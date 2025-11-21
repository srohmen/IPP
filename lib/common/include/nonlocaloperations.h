#ifndef NONLOCALOPERATIONS_H
#define NONLOCALOPERATIONS_H

#include <memory>

namespace IPP
{

class AbstractDistributeExcessTotals;
class AbstractCalcInterfaceProperties;

class NonLocalOperations
{
public:
    NonLocalOperations();
    ~NonLocalOperations();

    void setnComps(const size_t nComps);

    std::unique_ptr<AbstractCalcInterfaceProperties> interfaceProperties;
};

}


#endif // NONLOCALOPERATIONS_H
