#ifndef APPLYFUNCONTUPLEELEMENT_H
#define APPLYFUNCONTUPLEELEMENT_H

namespace IPP
{

template<typename Func, size_t pos>
class ApplyFuncOnTupleElement
{
public:
    ApplyFuncOnTupleElement(const Func& func)
        : m_func(func)
    {

    }

    template<typename Foo>
    bool operator()(const Foo& foo)
    {
        return m_func(std::get<pos>(foo));
    }

private:
    const Func& m_func;
};

}

#endif // APPLYFUNCONTUPLEELEMENT_H
