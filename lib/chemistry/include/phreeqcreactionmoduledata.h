#ifndef PHREEQCREACTIONMODULEDATA_H
#define PHREEQCREACTIONMODULEDATA_H

#include <cstddef>
#include <vector>

namespace IPP
{

class CellID_SI_Func;
class CellID_DissolveOnly_Func;
struct ReactionModuleOptimization;

class PhreeqcReactionModuleData
{
public:
    PhreeqcReactionModuleData();
    ~PhreeqcReactionModuleData();


    const CellID_SI_Func* siFunc;
    const CellID_DissolveOnly_Func* dissolveOnlyFunc;

    ReactionModuleOptimization* optim;

    std::vector<size_t> phaseToNuclPhaseID;
    std::vector<size_t> nuclPhaseToMonomerID;
    std::vector<double> monomerConc;

};

}

#endif // PHREEQCREACTIONMODULEDATA_H
