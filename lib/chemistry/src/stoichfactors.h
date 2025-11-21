#ifndef STOICHFACTORS_H
#define STOICHFACTORS_H

#include <vector>
#include <cstddef>

namespace IPP
{

typedef std::pair<size_t, size_t> SpeciesIDtoStoichFactor;
typedef std::vector<SpeciesIDtoStoichFactor> SpeciesIDtoStoichFactorVec;
typedef std::vector<SpeciesIDtoStoichFactorVec> StoichFactorsPerCompID;

}

#endif // STOICHFACTORS_H
