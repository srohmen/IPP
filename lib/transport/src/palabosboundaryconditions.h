#ifndef PALABOSBOUNDARYCONDITIONS_H
#define PALABOSBOUNDARYCONDITIONS_H

#include <vector>
#include <map>
#include <memory>
#include <unordered_set>


#include <palabos/boundaryCondition/boundaryCondition.h>


#include <palabos/dataProcessors/dataInitializerWrapper2D.h>

#include "plbtypededuction.h"
#include "boundaryconditiondata.h"
#include "palabosboundaryconditionutils.h"
#include "abstractboundaryconditions.h"
#include "palabosconversiontools.h"
#include "domaininfos.h"

namespace IPP
{

namespace PalabosBoundaryConditions
{

namespace BoundaryConditionSetup
{

    enum LBMomentType
    {
        LBMT_Density,
        LBMT_Velocity
    };


    typedef std::map<LBMomentType, plb::boundary::BcType> MomentumToBCType;
    MomentumToBCType convertToPlbBoundaryType(const BoundaryConditionType& inputType);


    DomainInfos::BoundaryEnvelops calcNeededEnvelops(const std::vector<size_t>& forcedDims,
                                                     const std::vector<BoundaryConditionData*>& bcs);


    template<typename BCSetup,
             typename T,
             template<typename U> class Descriptor,
             template<typename U, template<typename V> class Desc> class Lattice,
             typename TransFunc>
    static void defineBoundaryConditions(Lattice<T, Descriptor>& lattice,
                                         const std::vector<int>& periodicBC,
                                         const AbstractBoundaryConditions::BCVec& bcVec,
                                         const TransFunc& transFunc)
    {
        for(const AbstractBoundaryConditions::BoundaryConditionData& bc : bcVec)
        {
            const MomentumToBCType boundaryType = convertToPlbBoundaryType(bc.type);
            for(auto it = boundaryType.begin(); it != boundaryType.end(); ++it)
            {
                const LBMomentType moment = it->first;

                switch(moment)
                {

                using Box = typename PlbTypeDeduction::GetBoxXD<Descriptor<T>::d>::value;

                case (LBMT_Density):
                {
                    Box domain;
                    PalabosConversionTools::convertBox(bc.domain.range, domain);

                    BCSetup::setDensityBoundaryCondition(lattice, domain, periodicBC, it->second);

                    const auto value = transFunc.apply(bc.density);
                    plb::setBoundaryDensity(lattice, domain, (T)value);
                    break;
                }

                case (LBMT_Velocity):
                {
                    Box domain;
                    PalabosConversionTools::convertBox(bc.domain.range, domain);

                    BCSetup::setVelocityBoundaryCondition(lattice, domain, periodicBC, it->second);

                    using Vector= plb::Array<T, Descriptor<T>::d>;
                    Vector flowVec;
                    PalabosConversionTools::convertVector(bc.velocity, flowVec);
                    plb::setBoundaryVelocity(lattice, domain, flowVec);

                    break;
                }

                default:
                {
                    throw std::runtime_error("unknown moment type: " + std::to_string(moment));
                    break;
                }
                }
            }
        }
    }

}

}

}

#endif // PALABOSBOUNDARYCONDITIONS_H
