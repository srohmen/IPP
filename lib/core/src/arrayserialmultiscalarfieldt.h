#ifndef ARRAYSERIALMULTISCALARFIELDT_H
#define ARRAYSERIALMULTISCALARFIELDT_H

#include "serialmultiscalar.h"

namespace IPP
{

template<typename T, typename MultiField, typename SerialField, size_t dim>
class ArraySerialMultiScalarField : public SerialMultiScalarT<T, MultiField, SerialField, dim>
{
public:
    using BaseClass = SerialMultiScalarT<T, MultiField, SerialField, dim>;
    using ArrayType = std::vector<T>;

    ArraySerialMultiScalarField(MultiField&& multiData, SerialField&& serialData, ArrayType& rawArr)
        : BaseClass(std::move(multiData), std::move(serialData)),
          m_arr(rawArr)
    {

    }

    const ArrayType& getRawDataReadOnly() const
    {
        this->getSerialDataReadOnly();
        return m_arr;
    }

    ArrayType& getRawDataMutable()
    {
        this->getSerialDataMutable();
        return m_arr;
    }


private:
     ArrayType& m_arr;

};

}

#endif // ARRAYSERIALMULTISCALARFIELDT_H
