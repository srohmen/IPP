#ifndef SIMPLEPOROSITYINFOS_H
#define SIMPLEPOROSITYINFOS_H

#include <memory>

namespace IPP
{

struct SimplePorosityInfos
{
    virtual ~SimplePorosityInfos()
    {

    }

    double porosityTotal;
    double porosityCapillary;
};

typedef std::unique_ptr<SimplePorosityInfos> SimplePorosityInfosPtr;

}

#endif // SIMPLEPOROSITYINFOS_H
