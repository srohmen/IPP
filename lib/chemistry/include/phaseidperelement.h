#ifndef PHASEIDPERELEMENT_H
#define PHASEIDPERELEMENT_H

#include <vector>
#include <cstddef>

namespace IPP
{

class PhaseIdPerElement
{
public:
    PhaseIdPerElement();
    void init(const size_t nElements);

    void add(const size_t elementId, const size_t phaseId);

    size_t getNumberElements() const;

    typedef std::vector<size_t> PhaseIds;
    const PhaseIds& getPhaseIds(const size_t elementId) const;

private:
    typedef std::vector<PhaseIds> Lookup;
    Lookup lookup;

};

}

#endif // PHASEIDPERELEMENT_H
