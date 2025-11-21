#ifndef SCALARTRANSFORMATIONS_H
#define SCALARTRANSFORMATIONS_H

#include <memory>
#include <vector>
#include <palabos/dataProcessors/dataAnalysisWrapper2D.h>
#include <palabos/dataProcessors/dataAnalysisWrapper3D.h>

namespace IPP
{

namespace ScalarTransformations
{

struct NoTransform
{
    template<typename DataField>
    void apply(DataField&) const
    {
        // do nothing
    }

    template<typename T>
    T apply(const T& val) const
    {
        return val;
    }
};

class NoTransformationForwardBackward
{
public:
    NoTransformationForwardBackward()
        : forward(),
          backward()
    {

    }

    ~NoTransformationForwardBackward()
    {
    }

    const NoTransform forward;
    const NoTransform backward;

};

template<typename ... T>
using NoTransFwdBackPtr = std::shared_ptr<NoTransformationForwardBackward>;

class TransBase
{
protected:
    static const double m_offset;
};

class ForwardTrans : public TransBase
{
public:
    ForwardTrans()
    {

    }

    template<typename DataField>
    void apply(DataField& field) const
    {
        plb::addInPlace(field, this->m_offset);
    }

    template<typename T>
    T apply(const size_t x, const size_t y, const size_t z, const T& val) const
    {
        T result = val;

        result += this->m_offset;

        return result;
    }

    template<typename T>
    T apply(const T& val) const
    {
        T result = val + this->m_offset;
        return result;
    }

};

class BackwardTrans : public TransBase
{
public:
    BackwardTrans()
    {

    }

    template<typename DataField>
    void apply(DataField& field) const
    {
        plb::subtractInPlace(field, this->m_offset);
    }

    template<typename T>
    T apply(const size_t x, const size_t y, const size_t z, const T& val) const
    {
        T result = val;
        result -= this->m_offset;
        return result;
    }

    template<typename T>
    T apply(const T& val) const
    {
        T result = val - this->m_offset;
        return result;
    }

};


class TransformationForwardBackward
{
public:
    TransformationForwardBackward()
    {

    }

    ~TransformationForwardBackward()
    {
    }

    const ForwardTrans forward;
    const BackwardTrans backward;
};


using TransFwdBackPtr = std::shared_ptr<TransformationForwardBackward>;

//template<typename Scalar>
//class NoDimensionTransform : public TransformBase<Scalar>
//{
//public:

//    virtual Scalar operator()(const size_t, const size_t, const Scalar& val) const
//    {
//        return (*this)(val);
//    }

//    virtual Scalar operator()(const size_t, const size_t, const size_t, const Scalar& val) const
//    {
//        return (*this)(val);
//    }

//    virtual Scalar operator ()(const Scalar& val) const = 0;


//private:

//};

//template<typename Scalar>
//class CompositeTransformFwdBwd : public TransformBase<Scalar>
//{
//public:
//    CompositeTransformFwdBwd()
//    {

//    }

//    virtual ~CompositeTransformFwdBwd()
//    {
//        for(size_t i = 0; i < m_trans.size(); ++i)
//        {
//            delete m_trans[i];
//        }
//        m_trans.clear();
//    }

//    void push_back(const TransformBase<Scalar>* trans)
//    {
//        m_trans.push_back(trans);
//    }

//    virtual Scalar operator()(const size_t x, const size_t y, const Scalar& val) const
//    {
//        Scalar result = val;
//        for(size_t i = 0; i < m_trans.size(); ++i)
//        {
//            const TransformBase<Scalar>* trans = m_trans[i];
//            result = (*trans)(x, y, result);
//        }

//        return result;
//    }

//    virtual Scalar operator()(const size_t x, const size_t y, const size_t z, const Scalar& val) const
//    {
//        Scalar result = val;
//        for(size_t i = 0; i < m_trans.size(); ++i)
//        {
//            const TransformBase<Scalar>* trans = m_trans[i];
//            result = (*trans)(x, y, z, result);
//        }

//        return result;
//    }

//private:



//    std::vector<TransformBase<Scalar>*> m_trans;
//};



//template<typename Scalar>
//class CompositeTransform : public TransformationForwardBackward<Scalar>
//{
//public:
//    typedef TransformationForwardBackward<Scalar> BaseClass;

//    CompositeTransform()
//        : concreteFwd(new CompositeTransformFwdBwd<Scalar>),
//          concreteBwd(new CompositeTransformFwdBwd<Scalar>),
//          BaseClass(concreteFwd, concreteBwd)
//    {

//    }

//    virtual ~CompositeTransform()
//    {
//        // base class deletes the objects
//        concreteFwd = nullptr;
//        concreteBwd = nullptr;
//    }

//    void add(const TransformBase<Scalar>* fwdTrans,
//             const TransformBase<Scalar>* backTrans)
//    {
//        concreteFwd->push_back(fwdTrans);
//        concreteBwd->push_back(backTrans);
//    }

//private:
//    CompositeTransformFwdBwd<Scalar>* concreteFwd;
//    CompositeTransformFwdBwd<Scalar>* concreteBwd;

//};

//template<typename T>
//using CompositeTransformPtr = std::shared_ptr<CompositeTransform<T>>;

//template<typename T>
//inline CompositeTransformPtr<T> makeCompositeTransform()
//{
//    return std::make_shared<CompositeTransformPtr<T>>;
//}





//template<typename Scalar>
//class NoTransformForward : public NoDimensionTransform<Scalar>
//{
//public:
//    virtual Scalar operator()(const Scalar& val) const
//    {
//        return val;
//    }
//};

//template<typename Scalar>
//using NoTransformBackward = NoTransformForward<Scalar>;

//template<typename Scalar>
//class NoTransform : public TransformationForwardBackward<Scalar>
//{
//public:
//    typedef TransformationForwardBackward<Scalar> BaseClass;
//    NoTransform()
//        : BaseClass(new NoTransformForward<Scalar>(),
//                    new NoTransformBackward<Scalar>())
//    {

//    }
//};

//template<typename T>
//using NoTransformPtr = std::shared_ptr<NoTransform<T>>;

//template<typename T>
//inline NoTransformPtr<T> makeNoTransform()
//{
//    return std::make_shared<NoTransform<T>>();
//}





//template<typename Scalar>
//class AddTransform : public NoDimensionTransform<Scalar>
//{
//public:
//    AddTransform(const Scalar& offset)
//        : offset(offset)
//    {

//    }

//    virtual Scalar operator()(const Scalar& val) const
//    {
//        return val + offset;
//    }

//private:
//    const Scalar offset;
//};

//template<typename Scalar>
//class SubstractTransform : public NoDimensionTransform<Scalar>
//{
//public:
//    SubstractTransform(const Scalar& offset)
//        : offset(offset)
//    {

//    }

//    virtual Scalar operator()(const Scalar& val) const
//    {
//        return val - offset;
//    }

//private:
//    const Scalar offset;
//};

//template<typename Scalar>
//class OffsetTransform : public TransformationForwardBackward<Scalar>
//{
//public:
//    typedef TransformationForwardBackward<Scalar> BaseClass;

//    OffsetTransform(const Scalar& offset)
//        : BaseClass(new AddTransform<Scalar>(offset),
//                    new SubstractTransform<Scalar>(offset))
//    {

//    }

//};

//template<typename T>
//using OffsetTransformPtr = std::shared_ptr<OffsetTransform<T>>;

//template<typename T>
//inline OffsetTransformPtr<T> makeOffsetTransform(const T& offset)
//{
//    return std::make_shared<OffsetTransform<T>>(offset);
//}



//template<typename Scalar>
//struct LinearTransformData
//{
//    LinearTransformData(const Scalar& srcMin, const Scalar& dstMin, const Scalar& factor)
//        : srcMin(srcMin),
//          dstMin(dstMin),
//          factor(factor)
//    {

//    }

//    const Scalar srcMin;
//    const Scalar dstMin;
//    const Scalar factor;
//};

//template<typename Scalar>
//class LinearTransformForward : public NoDimensionTransform<Scalar>
//{
//public:
//    LinearTransformForward(const LinearTransformData<Scalar>& data)
//        : data(data)
//    {

//    }

//    virtual Scalar operator()(const Scalar& val) const
//    {
//        const Scalar result = (val - data.srcMin) * data.factor + data.dstMin;
//        return result;
//    }

//private:
//    const LinearTransformData<Scalar>& data;
//};


//template<typename Scalar>
//class LinearTransformBackward : public NoDimensionTransform<Scalar>
//{
//public:
//    LinearTransformBackward(const LinearTransformData<Scalar>& data)
//        : data(data)
//    {

//    }

//    virtual Scalar operator()(const Scalar& val) const
//    {
//        const Scalar result = (val - data.dstMin) / data.factor + data.srcMin;
//        return result;
//    }


//private:
//    const LinearTransformData<Scalar>& data;
//};

//template<typename Scalar>
//class LinearTransform : public TransformationForwardBackward<Scalar>
//{
//public:
//    typedef TransformationForwardBackward<Scalar> BaseClass;

//    // TODO: check init order
//    LinearTransform(const Scalar& srcMin, const Scalar& dstMin, const Scalar& factor)
//        : BaseClass(new LinearTransformForward<Scalar>(data),
//                    new LinearTransformBackward<Scalar>(data)),
//          data(srcMin, dstMin, factor)
//    {

//    }


//private:
//    const LinearTransformData<Scalar> data;
//};

//template<typename T>
//using LinearTransformPtr = std::shared_ptr<LinearTransform<T>>;



//template<typename T>
//inline LinearTransformPtr<T> makeLinearTransform(const T& srcMin, const T& srcMax, const T& dstMin, const T& dstMax)
//{
//    assert(srcMin < srcMax);
//    const T srcDiff = srcMax - srcMin;

//    assert(dstMin < dstMax);
//    const T dstDiff = dstMax - dstMin;

//    const T factor = dstDiff / srcDiff;

//    return std::make_shared<LinearTransform<T>>(srcMin, dstMin, factor);
//}





//template<typename Scalar, typename Array>
//class MultiplyTransform : public TransformBase<Scalar>
//{
//public:
//    MultiplyTransform(const Array& arr)
//        : m_arr(arr)
//    {

//    }

//    virtual Scalar operator()(const size_t x, const size_t y, const Scalar& val) const
//    {
//        return val * m_arr.get(x, y, 0);
//    }

//    virtual Scalar operator()(const size_t x, const size_t y, const size_t z, const Scalar& val) const
//    {
//        return val * m_arr.get(x, y, z);
//    }

//private:
//    const Array& m_arr;
//};

//template<typename Scalar, typename Array>
//class DivideTransform : public TransformBase<Scalar>
//{
//public:
//    DivideTransform(const Array& arr)
//        : m_arr(arr)
//    {

//    }

//    virtual Scalar operator()(const size_t x, const size_t y, const Scalar& val) const
//    {
//        return val / m_arr.get(x, y, 0);
//    }

//    virtual Scalar operator()(const size_t x, const size_t y, const size_t z, const Scalar& val) const
//    {
//        return val / m_arr.get(x, y, z);
//    }


//private:
//    const Array& m_arr;
//};

} // end of namespace ScalarTransformations

} // end of namespace LBGeoChem

#endif // SCALARTRANSFORMATIONS_H
