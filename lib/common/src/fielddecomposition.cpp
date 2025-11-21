#include "fielddecomposition.h"

#include <cassert>

namespace IPP
{

FieldDecomposition::FieldDecomposition()
{

}

FieldDecomposition::~FieldDecomposition()
{

}

void FieldDecomposition::correctByTensorDim(const size_t tensorDim)
{
    const size_t nProcs = displs.size();
    assert(nProcs == sndCounts.size());

    int prevDispl = 0;
    for(size_t iProc = 0; iProc < nProcs; ++iProc)
    {
        displs[iProc] = prevDispl;
        sndCounts[iProc] *= tensorDim;
        prevDispl += sndCounts[iProc];
    }
}

const IPPVector3DLong &FieldDecomposition::getGlobalSize() const
{
    return globalSize;
}

size_t FieldDecomposition::getGlobalNumberCells() const
{
    return globalSize[0] * globalSize[1] * globalSize[2];
}

size_t FieldDecomposition::getLocalNumberOfCells() const
{
    return localNumberCells;
}

const IPPBox3DLong &FieldDecomposition::getOwnDomain() const
{
    return ownDomain;
}

const std::vector<IPPBox3DLong>& FieldDecomposition::getDomains() const
{
    return domains;
}

const std::vector<int> &FieldDecomposition::getDispls() const
{
    return displs;
}

const std::vector<int> &FieldDecomposition::getSndCounts() const
{
    return sndCounts;
}


}
