#ifndef CELLID_DISSOLVEONLY_FUNC_H
#define CELLID_DISSOLVEONLY_FUNC_H

#include <vector>
#include <string>
#include <cstddef>

namespace IPP
{

class CellID_DissolveOnly_Func
{
public:
    typedef unsigned char return_type;
    typedef std::vector<return_type>::iterator iterator;

    CellID_DissolveOnly_Func();
    virtual ~CellID_DissolveOnly_Func();

    virtual double getLowerPorosityThresh() const = 0;

    virtual bool needsNeighInfos() const = 0;
    virtual bool needsSaturationIndices() const = 0;

    virtual bool evaluate(const size_t iCell, const iterator& begin, const iterator& end) const = 0;
};

}

#endif // CELLID_DISSOLVEONLY_FUNC_H
