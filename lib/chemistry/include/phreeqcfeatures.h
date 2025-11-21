#ifndef PHREEQCFEATURES_H
#define PHREEQCFEATURES_H

#include <unordered_map>
#include <vector>

namespace IPP
{

enum PhreeqcFeatures
{
    PF_SOLUTION             = (1u << 0),
    PF_EQUILIBRIUM_PHASES   = (1u << 1),
    PF_EXCHANGE             = (1u << 2),
    PF_SURFACE              = (1u << 3),
    PF_GAS_PHASE            = (1u << 4),
    PF_SOLID_SOLUTION       = (1u << 5),
    PF_KINETICS             = (1u << 6),

};

static const std::vector<PhreeqcFeatures> s_allFeatures =
{{
     PF_SOLUTION,
     PF_EQUILIBRIUM_PHASES,
     PF_EXCHANGE,
     PF_SURFACE,
     PF_GAS_PHASE,
     PF_SOLID_SOLUTION,
     PF_KINETICS

}};


typedef std::unordered_map<size_t, unsigned int> CellIdToPhreeqcFeatureMask;

}

#endif // PHREEQCFEATURES_H
