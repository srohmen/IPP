#include "anyelementunbuffered.h"

#include <cassert>

#include "phreeqcconstants.h"


namespace IPP
{

static PhaseIdPerElement stripToPrimary(const PhaseIdPerElement& phaseIdPerElement)
{
    // strip (H2O), H, O and Charge (0-2)/(0-3) from lookup table
    const size_t nElementsTotal = phaseIdPerElement.getNumberElements();
    const size_t nPrimaryElements = nElementsTotal - PhreeqcConstants::s_primaryCompsBegin;

    PhaseIdPerElement stripped;
    stripped.init(nPrimaryElements);

    for(size_t iElementTotal = PhreeqcConstants::s_primaryCompsBegin;
        iElementTotal < nElementsTotal; ++iElementTotal)
    {
        const PhaseIdPerElement::PhaseIds& ids = phaseIdPerElement.getPhaseIds(iElementTotal);

        const size_t iPrimaryElement = iElementTotal - PhreeqcConstants::s_primaryCompsBegin;

        for(size_t phaseID : ids)
        {
            stripped.add(iPrimaryElement, phaseID);
        }
    }

    return stripped;
}

AnyElementUnbuffered::AnyElementUnbuffered(const PhaseIdPerElement& phaseIdPerElement)
    : m_phaseIdPerElement(stripToPrimary(phaseIdPerElement))
{

}

bool AnyElementUnbuffered::run(const Iterator begin, const Iterator end) const
{
    const size_t nElements = m_phaseIdPerElement.getNumberElements();
    for(size_t iElement = 0; iElement < nElements; ++iElement)
    {
        const bool isBuffered = isElementInAnyPhase(iElement, begin, end);

        if(isBuffered == false)
        {
            return true;
        }
    }

    return false;
}

bool AnyElementUnbuffered::isElementInAnyPhase(const size_t elementID,
                                               const Iterator begin,
                                               const Iterator end) const
{
    const PhaseIdPerElement::PhaseIds& phaseIDs = m_phaseIdPerElement.getPhaseIds(elementID);

    for(size_t phaseID : phaseIDs)
    {
        assert(begin + phaseID < end);
        Iterator it = begin;
        std::advance(it, phaseID);

        const double& amount = *it;

        if(amount > 0.0)
        {
            return true;
        }
    }

    return false;
}


}
