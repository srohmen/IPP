#ifndef SETUPBLOCKPERIODICITY_H
#define SETUPBLOCKPERIODICITY_H

#include <vector>

namespace IPP
{

namespace SetupBlockPeriodicity
{

template<typename Block>
void definePeriodicity(const std::vector<int>& periodicBC, Block& block)
{
    block.periodicity().toggleAll(false);

    for(const size_t d : periodicBC)
    {
        block.periodicity().toggle(d, true);
    }
}


}

}

#endif // SETUPBLOCKPERIODICITY_H
