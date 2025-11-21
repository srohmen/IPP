#include "nonlocaloperations.h"

#include <cassert>
#include "abstractdistributeexcesstotals.h"
#include "abstractcalcinterfaceproperties.h"

namespace IPP
{

NonLocalOperations::NonLocalOperations() = default;
NonLocalOperations::~NonLocalOperations() = default;

void NonLocalOperations::setnComps(const size_t nComps)
{
    assert(interfaceProperties);
    interfaceProperties->setnComps(nComps);
}
}
