#ifndef ABSTRACTDISSOLVEONLYCALC_H
#define ABSTRACTDISSOLVEONLYCALC_H

#include <vector>
#include <string>

#include "dissprecipbehaviour.h"
#include "phasenametoinfos.h"

namespace IPP
{

struct DissolveOnlyInfos;

class AbstractDissolveOnlyCalc
{
public:
    AbstractDissolveOnlyCalc() = default;

    virtual ~AbstractDissolveOnlyCalc() = default;

    virtual void init(const PhaseNameToInfos &phaseInfos,
                      const double& cellVol,
                      const double& cellArea,
                      const double& dt) = 0;

    virtual double getLowerPorosityThresh() const = 0;

    virtual bool needsNeighInfos() const = 0;
    virtual bool needsSaturationIndices() const = 0;
    virtual bool isNucleation() const = 0;

    virtual void findNucleationPhases(std::vector<std::string>& nucleationPhases,
                                      std::vector<std::string>& monomers) const = 0;
    virtual void findNucleationMonomers(std::vector<std::string>& nucleationMonomers) const = 0;

    virtual void setPrecipDissolveOnlyBehaviour(const std::vector<DissPrecipBehaviour>& behav) = 0;

    typedef std::vector<unsigned char>::iterator iterator;
    virtual bool evaluate(const DissolveOnlyInfos& data,
                          const iterator& begin,
                          const iterator& end) const = 0;

};

}

#endif // ABSTRACTDISSOLVEONLYCALC_H
