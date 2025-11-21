#ifndef MPIDELEGATE_H
#define MPIDELEGATE_H

#include <PhreeqcRM.h>

namespace IPP
{

class MPIDelegate
{
public:
    MPIDelegate(PhreeqcRM& phreeqc)
        : m_phreeqc(phreeqc)
    {
        if(m_phreeqc.GetMpiMyself() != 0)
        {
            m_phreeqc.MpiWorker();
        }
    }

    ~MPIDelegate()
    {
        if(m_phreeqc.GetMpiMyself() == 0)
        {
            m_phreeqc.MpiWorkerBreak();
        }
    }

private:
    PhreeqcRM& m_phreeqc;
};

template<typename Func>
auto phreeqcMPIDispatch(PhreeqcRM& phreeqc, const Func& func ) -> decltype(func())
{
    const MPIDelegate delegate(phreeqc);

    if(phreeqc.GetMpiMyself() == 0)
    {
        return func();
    }
    else
    {
        typedef decltype(func()) RetType;
        return RetType();
    }
}

}

#endif // MPIDELEGATE_H
