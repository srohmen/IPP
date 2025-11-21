#include "arrayviewiterator.h"

#include "arrayviewiteratorimpl.h"

namespace IPP
{

ArrayViewIterator::ArrayViewIterator()
    : m_pImpl(nullptr)
{

}

ArrayViewIterator::ArrayViewIterator(const ArrayViewIterator& other)
    : m_pImpl(new ArrayViewIteratorImpl(*other.m_pImpl))
{

}


ArrayViewIterator::ArrayViewIterator(ArrayViewIteratorImpl *pImpl)
    : m_pImpl(pImpl)
{

}

ArrayViewIterator::~ArrayViewIterator()
{
    delete m_pImpl;
    m_pImpl = nullptr;
}

ArrayViewIterator& ArrayViewIterator::operator =(const ArrayViewIterator &other)
{
    m_pImpl = new ArrayViewIteratorImpl(*other.m_pImpl);
    return *this;
}

ArrayViewIterator &ArrayViewIterator::operator ++()
{
    ++(*m_pImpl);
    return *this;
}

bool ArrayViewIterator::operator ==(const ArrayViewIterator &other)
{
    return *m_pImpl == *other.m_pImpl;
}

bool ArrayViewIterator::operator !=(const ArrayViewIterator &other)
{
    return *m_pImpl != *other.m_pImpl;
}

const double &ArrayViewIterator::operator *()
{
    return *(*m_pImpl);
}

}

