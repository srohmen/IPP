#ifndef REACTIONMODULEOPTIMIZATION_H
#define REACTIONMODULEOPTIMIZATION_H

#include <cstddef>

namespace IPP
{

class IsInsideVolume;
class SignificantChangePrediction;

struct ReactionModuleOptimization
{
    ReactionModuleOptimization();

    ~ReactionModuleOptimization();

    SignificantChangePrediction* changePredict;
    size_t forceRecalcFreq;

};

}

#endif // REACTIONMODULEOPTIMIZATION_H
