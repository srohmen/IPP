#ifndef CELLIDPASSINGFUNCWRAPPER_H
#define CELLIDPASSINGFUNCWRAPPER_H

#include <vector>
#include <cassert>

namespace IPP
{

template<typename Base, typename Func>
class CellIDPassingFuncWrapper : public Base
{
public:
    CellIDPassingFuncWrapper(const Func &func)
        : m_func(func)
    {

    }

    virtual size_t nCells() const
    {
        assert(false);
        return 0;
    }

    typedef typename Base::return_type return_type;
    virtual void evaluate(const size_t iCell, std::vector<return_type>& values) const
    {
        m_func.evaluate(iCell, values);
    }

private:
    const Func& m_func;
};

}

#endif // CELLIDPASSINGFUNCWRAPPER_H
