#include "reactionmoduleoptimization.h"

#include "isinsidevolume.h"
#include "significantchangeprediction.h"

namespace IPP
{

ReactionModuleOptimization::ReactionModuleOptimization()
    : changePredict(nullptr)
    , forceRecalcFreq(0)
{

}

ReactionModuleOptimization::~ReactionModuleOptimization()
{
    delete changePredict;
    changePredict = nullptr;
}

}

