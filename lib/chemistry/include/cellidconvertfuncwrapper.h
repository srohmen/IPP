#ifndef CELLIDCONVERTFUNCWRAPPER_H
#define CELLIDCONVERTFUNCWRAPPER_H


#include <vector>

#include "arraydimensionconvert.h"

namespace IPP
{

template<typename Base, typename Func>
class CellIDConvertFuncWrapper : public Base
{
public:
    CellIDConvertFuncWrapper(const ArrayDimensionConvert* indexConv,
                             const Func &func)
        : m_indexConv(indexConv),
          m_func(func)
    {

    }

    virtual ~CellIDConvertFuncWrapper()
    {
        delete m_indexConv;
        m_indexConv = nullptr;
    }

    virtual size_t nCells() const
    {
        return m_indexConv->getNxyz();
    }

    typedef typename Base::return_type return_type;

    virtual void evaluate(const size_t iCell,
                          const typename Base::iterator& begin,
                          const typename Base::iterator& end) const
    {
        const IPPVector3DInt coord = m_indexConv->calcCoordinate(iCell);
        m_func.evaluate(coord, begin, end);
    }

//    virtual void evaluate(const size_t iCell,
//                          std::vector<return_type>& values) const
//    {
//        const IPPVector3DInt coord = m_indexConv->calcCoordinate(iCell);
//        m_func.evaluate(coord, values);
//    }



private:
    const ArrayDimensionConvert* m_indexConv;
    const Func& m_func;
};


}

#endif // CELLIDCONVERTFUNCWRAPPER_H
