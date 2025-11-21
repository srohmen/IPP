#ifndef FIELDDECOMPOSITION_H
#define FIELDDECOMPOSITION_H

#include <vector>
#include <cstddef>

#include "ippbox.h"

namespace boost
{
namespace serialization
{
class access;
}
}

namespace IPP
{

class FieldDecomposition
{
public:
    FieldDecomposition();

    virtual ~FieldDecomposition();

    bool operator==(const FieldDecomposition& other) const
    {
        const bool result
                = globalSize == other.globalSize
                && localNumberCells == other.localNumberCells
                && ownDomain == other.ownDomain
                && domains == other.domains
                && displs == other.displs
                && sndCounts == other.sndCounts;
        return result;
    }

    bool operator!=(const FieldDecomposition& other) const
    {
        return !(*this == other);
    }

    void correctByTensorDim(const size_t tensorDim);

    const IPPVector3DLong& getGlobalSize() const;
    size_t getGlobalNumberCells() const;

    size_t getLocalNumberOfCells() const;
    const IPPBox3DLong& getOwnDomain() const;

    const std::vector<IPPBox3DLong>& getDomains() const;
    const std::vector<int>& getDispls() const;
    const std::vector<int>& getSndCounts() const;


protected:
    IPPVector3DLong globalSize;
    size_t localNumberCells;
    IPPBox3DLong ownDomain;
    std::vector<IPPBox3DLong> domains;
    std::vector<int> displs;
    std::vector<int> sndCounts;

private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & globalSize;
        ar & localNumberCells;
        ar & ownDomain;
        ar & domains;
        ar & displs;
        ar & sndCounts;
    }

};


}

#endif // FIELDDECOMPOSITION_H
