#include "palabosboundaryconditions.h"

namespace IPP
{

namespace PalabosBoundaryConditions
{


typedef std::pair<BoundaryConditionDomain::BoundaryPosition,
BoundaryConditionDomain::BoundaryPosition> BoundaryPositionPair;
static BoundaryPositionPair getBCFromDim(const size_t dim)
{
    switch(dim)
    {
        case 0:
        {
            return std::make_pair(BoundaryConditionDomain::BP_left,
                                  BoundaryConditionDomain::BP_right);
            break;
        }

        case 1:
        {
            return std::make_pair(BoundaryConditionDomain::BP_bottom,
                                  BoundaryConditionDomain::BP_top);
            break;
        }

        case 2:
        {
            return std::make_pair(BoundaryConditionDomain::BP_back,
                                  BoundaryConditionDomain::BP_front);
            break;
        }

        default:
        {
            throw std::runtime_error("unknown dimension: " + std::to_string(dim));
        }
    }

}



static bool needsExtend(const BoundaryConditionType& type)
{
    switch(type)
    {
        case (BCT_Closed):
        {
            return false;
            break;
        }

        case (BCT_DensityDirichlet):
        case (BCT_VelocityDirichlet):
        case (BCT_DensityCauchy):
        case (BCT_VelocityOutflow):
        {

            return true;
            break;
        }

        default:
        {
            throw std::runtime_error("unknown boundary type: " + std::to_string(type));
            break;
        }
    }
}


BoundaryConditionSetup::MomentumToBCType
BoundaryConditionSetup::convertToPlbBoundaryType(const BoundaryConditionType& inputType)
{
    MomentumToBCType ret;

    switch(inputType)
    {
        case (BCT_Closed):
        {
            // impose no BC because non-wet bounce back is more suitable
            // ret[BCT_Density] = plb::boundary::neumann
            break;
        }

        case (BCT_DensityDirichlet):
        {
            ret[LBMT_Density] = plb::boundary::dirichlet;
            break;
        }

        case (BCT_VelocityDirichlet):
        {
            ret[LBMT_Velocity] = plb::boundary::dirichlet;
            break;
        }

        case (BCT_DensityCauchy):
        {
            ret[LBMT_Density] = plb::boundary::dirichlet;
            ret[LBMT_Velocity] = plb::boundary::dirichlet;
            break;
        }

        case (BCT_VelocityOutflow):
        {
            ret[LBMT_Velocity] = plb::boundary::outflow;
            break;
        }

        default:
        {
            throw std::runtime_error("unknown boundary type: " + std::to_string(inputType));
            break;
        }
    }

    return ret;
}

DomainInfos::BoundaryEnvelops BoundaryConditionSetup::calcNeededEnvelops(const std::vector<size_t>& forcedDims,
                                                                         const std::vector<BoundaryConditionData*>& bcs)
{
    DomainInfos::BoundaryEnvelops needsExtendSet;

    for(const size_t dim : forcedDims)
    {
        const BoundaryPositionPair pair = getBCFromDim(dim);
        needsExtendSet.insert(pair.first);
        needsExtendSet.insert(pair.second);
    }

    for(const BoundaryConditionData* bc : bcs)
    {
        const bool needsEvelop = needsExtend(bc->type);
        if(needsEvelop)
        {
            needsExtendSet.insert(bc->domain.position);
        }
    }

    return needsExtendSet;
}



}
}
