#ifndef SERIALMULTISCALAR_H
#define SERIALMULTISCALAR_H

#include "lazydatawrapper.h"
#include "multitoserialsync.h"

namespace IPP
{

template<typename T, typename MultiData, typename SerialData, size_t dim>
class SerialMultiScalarT
{
public:

    using MultiDataField = MultiData;
    using SerialDataField = SerialData;

    SerialMultiScalarT()
        : m_multiDataProxy(0, 0),
          m_serialDataProxy()
    {
    }

    SerialMultiScalarT(MultiDataField&& multiData, SerialDataField&& serialData)
        : m_multiData(multiData),
          m_serialData(serialData),
          m_multiDataProxy(m_multiData, typename LazyMultiDataWrapper::Sync(m_serialData)),
          m_serialDataProxy(m_serialData, typename LazySerialDataWrapper::Sync(m_multiData))
    {
        // initially the serial data are marked as dirty
        m_serialDataProxy.setDirty();
    }

    void sync()
    {
        assert(m_serialDataProxy.isDirty() == true || m_multiDataProxy.isDirty() == true);

        if(m_serialDataProxy.isDirty())
        {
            m_serialDataProxy.sync();
        }
        else if(m_multiDataProxy.isDirty())
        {
            m_multiDataProxy.sync();
        }

        assert(m_serialDataProxy.isDirty() == false && m_multiDataProxy.isDirty() == false);
    }

    SerialDataField& getSerialDataReadOnly()
    {
        return m_serialDataProxy.get();
    }

    const SerialDataField& getSerialDataReadOnly() const
    {
        return m_serialDataProxy.get();
    }


    SerialDataField& getSerialDataMutable()
    {
        m_multiDataProxy.setDirty();
        return m_serialDataProxy.get();
    }

    MultiDataField& getMultiDataReadOnly()
    {
        return m_multiDataProxy.get();
    }

    const MultiDataField& getMultiDataReadOnly() const
    {
        return m_multiDataProxy.get();
    }

    MultiDataField& getMultiDataMutable()
    {
        m_serialDataProxy.setDirty();
        return m_multiDataProxy.get();
    }

private:
    MultiDataField m_multiData;
    SerialDataField m_serialData;

    using LazyMultiDataWrapper  = LazyDataWrapper<MultiDataField,   SyncSerialToMulti<SerialDataField>  >;
    using LazySerialDataWrapper = LazyDataWrapper<SerialDataField,  SyncMultiToSerial<MultiDataField>   >;

    LazyMultiDataWrapper m_multiDataProxy;
    LazySerialDataWrapper m_serialDataProxy;



};


}

#endif // SERIALMULTISCALAR_H
