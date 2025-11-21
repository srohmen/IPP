#ifndef ISENABLEDINCONTAINER_H
#define ISENABLEDINCONTAINER_H

#include <assert.h>
#include <cstddef>

namespace IPP
{

template<typename Container>
class IsEnabledInContainer
{
public:
    IsEnabledInContainer(const Container& isEnabled)
        : m_isEnabled(isEnabled)
    {

    }

    bool operator()(const size_t index) const
    {
        assert(index < m_isEnabled.size());
        return m_isEnabled[index];
    }

private:
    const Container& m_isEnabled;
};

template<typename Func>
class WrapIndexFunc
{
public:
    WrapIndexFunc(const Func& func, const size_t n)
        : m_func(func),
          m_n(n)
    {

    }

    bool operator()(const size_t index) const
    {
        const size_t wrapped = index % m_n;
        return m_func(wrapped);
    }

private:
    const Func& m_func;
    const size_t m_n;
};

}


#endif // ISENABLEDINCONTAINER_H
