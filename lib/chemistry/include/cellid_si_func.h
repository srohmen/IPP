#ifndef CELLID_SI_FUNC_H
#define CELLID_SI_FUNC_H

#include <cstddef>
#include <vector>

namespace IPP
{

class CellID_SI_Func
{
public:
    typedef double return_type;

    CellID_SI_Func()
    {

    }
    virtual ~CellID_SI_Func()
    {

    }

//    virtual size_t nCells() const = 0;

    typedef std::vector<return_type>::iterator iterator;
    virtual void evaluate(const size_t iCell, const iterator& begin, const iterator& end) const = 0;

};

}

#endif // CELLID_SI_FUNC_H
