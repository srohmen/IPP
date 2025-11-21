#ifndef ARRAYVIEWITERATOR_H
#define ARRAYVIEWITERATOR_H

namespace IPP
{


class ArrayViewIteratorImpl;

class ArrayViewIterator
{
public:
    ArrayViewIterator();
    ArrayViewIterator(const ArrayViewIterator &other);
    ArrayViewIterator(ArrayViewIteratorImpl* pImpl);
    ~ArrayViewIterator();

    ArrayViewIterator &operator =(const ArrayViewIterator &other);
    ArrayViewIterator& operator ++();
    bool operator ==(const ArrayViewIterator &other);
    bool operator !=(const ArrayViewIterator &other);
    const double& operator *();



private:

    ArrayViewIteratorImpl* m_pImpl;
};

}

#endif // ARRAYVIEWITERATOR_H
