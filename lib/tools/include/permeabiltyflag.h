#ifndef PERMEABILTYFLAG_H
#define PERMEABILTYFLAG_H

namespace IPP
{

enum PermeabiltyFlag
{
    PF_isPermeable  = (1u << 0),
    PF_isInterface  = (1u << 1),
    PF_isInertSolid = (1u << 2)
};

}

#endif // PERMEABILTYFLAG_H
