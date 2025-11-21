#ifndef LAZYDATAWRAPPER_H
#define LAZYDATAWRAPPER_H

namespace IPP
{

template<typename LazyContainer, typename SyncFunc>
class LazyDataWrapper
{
public:

    using Sync = SyncFunc;

    LazyDataWrapper(LazyContainer& data, const SyncFunc& func)
        : m_func(func),
          m_isDirty(false),
          m_data(data)
    {

    }

    LazyContainer& get()
    {
        return this->getImpl();
    }

    const LazyContainer& get() const
    {
        return this->getImpl();
    }

    void setDirty()
    {
        m_isDirty = true;
    }

    bool isDirty() const
    {
        return m_isDirty;
    }

    void sync()
    {
        if(m_isDirty)
        {
            m_func(m_data);
            m_isDirty = false;
        }
    }

private:

    LazyContainer& getImpl()
    {
        if(m_isDirty)
        {
            m_func(m_data);
            m_isDirty = false;
        }

        return m_data;
    }

    const LazyContainer& getImpl() const
    {
        if(m_isDirty)
        {
            m_func(m_data);
            m_isDirty = false;
        }

        return m_data;
    }

    const SyncFunc m_func;
    mutable bool m_isDirty;
    LazyContainer& m_data;
};


template<template <typename T> class SyncFuncT, typename LazyContainer, typename Container>
LazyDataWrapper<LazyContainer, SyncFuncT<Container> >
makeLazyDataWrapper(LazyContainer& lazyContainer, Container& originContainer)
{
    using SyncFunc = SyncFuncT<LazyContainer>;
    return LazyDataWrapper<LazyContainer, SyncFunc>(lazyContainer, SyncFunc(originContainer));
}

}

#endif // LAZYDATAWRAPPER_H
