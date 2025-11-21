#ifndef ANYELEMENTUNBUFFERED_H
#define ANYELEMENTUNBUFFERED_H

#include <vector>

#include "phaseidperelement.h"

namespace IPP
{

class AnyElementUnbuffered
{
public:

    AnyElementUnbuffered(const PhaseIdPerElement& phaseIdPerElement);


    typedef std::vector<double>::const_iterator Iterator;

    bool run(const Iterator begin,
             const Iterator end) const;

private:
    bool isElementInAnyPhase(const size_t elementID,
                           const Iterator begin,
                           const Iterator end) const;

    const PhaseIdPerElement m_phaseIdPerElement;


};


}


#endif // ANYELEMENTUNBUFFERED_H
