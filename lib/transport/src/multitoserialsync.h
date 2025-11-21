#ifndef MULTITOSERIALSYNC_H
#define MULTITOSERIALSYNC_H

#include "multitoserialconversion.h"
#include "serialtomulticonversion.h"

namespace IPP
{

template<typename MultiField>
class SyncMultiToSerial
{
public:
    SyncMultiToSerial(MultiField& multiData)
        : m_multiData(multiData)
    {

    }

    template<typename SerialField>
    void operator()(SerialField& serialData) const
    {
        using T = typename SerialField::value_type;
        constexpr size_t dim = PlbTypeDeduction::GetFieldDim<T, MultiField>::value;
        using Sync = MultiToSerialConversion::ScalarFieldConvertT<T, SerialField, dim>;
        plb::applyProcessingFunctional(new Sync(serialData),
                                       m_multiData.getBoundingBox(), m_multiData);
    }

private:
    MultiField& m_multiData;

};

template<typename MultiField>
SyncMultiToSerial<MultiField> makeSyncMultiToSerial(MultiField& data)
{
    return SyncMultiToSerial<MultiField>(data);
}

template<typename SerialField>
class SyncSerialToMulti
{
public:
    SyncSerialToMulti(const SerialField& serialData)
        : m_serialData(serialData)
    {

    }

    template<typename MultiField>
    void operator()(MultiField& multiData) const
    {
        using T = typename SerialField::value_type;
        constexpr size_t dim = PlbTypeDeduction::GetFieldDim<T, MultiField>::value;
        using Sync = SerialToMultiConversion::ScalarFieldConvertT<T, SerialField, dim>;
        plb::applyProcessingFunctional(new Sync(m_serialData),
                                       multiData.getBoundingBox(), multiData);
    }

private:
    const SerialField& m_serialData;
};


}

#endif // MULTITOSERIALSYNC_H
