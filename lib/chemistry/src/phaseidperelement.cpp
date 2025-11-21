#include "phaseidperelement.h"

#include <algorithm>
#include <cassert>

namespace IPP
{

PhaseIdPerElement::PhaseIdPerElement()
{

}

void PhaseIdPerElement::init(const size_t nElements)
{
    lookup.resize(nElements);
}

void PhaseIdPerElement::add(const size_t elementId, const size_t phaseId)
{
    assert(lookup.size() > elementId);
    PhaseIds& ids = lookup[elementId];
    ids.push_back(phaseId);

    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()),
              ids.end());
}

size_t PhaseIdPerElement::getNumberElements() const
{
    return lookup.size();
}

const PhaseIdPerElement::PhaseIds&PhaseIdPerElement::getPhaseIds(const size_t elementId) const
{
    assert(lookup.size() > elementId);
    return lookup[elementId];
}

}
