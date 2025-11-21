#ifndef CSHPOROSITYINFOS_H
#define CSHPOROSITYINFOS_H

#include "simpleporosityinfos.h"

namespace IPP
{

struct CSHPorosityInfos : public SimplePorosityInfos
{
    double volNonPerm;
    double volSatCSH;
    double volFracCSH_HD;
    double ratioCaSi;
};

}

#endif // CSHPOROSITYINFOS_H
