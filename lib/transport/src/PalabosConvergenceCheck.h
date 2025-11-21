#ifndef PALABOSCONVERGENCECHECK_H
#define PALABOSCONVERGENCECHECK_H

#include <vector>
#include <palabos/dataProcessors/dataAnalysisWrapper2D.hh>
#include <palabos/dataProcessors/dataAnalysisWrapper3D.hh>
#include "plbtypededuction.h"
#include "domainiterator.h"
#include "ippstream.h"

namespace IPP
{

template<typename T, size_t dim>
class ComputeMaxDiff : public PlbTypeDeduction::
        GetReductiveBoxProcessingFunctionalXD_SS<T,T,dim>::value
{
public:
    ComputeMaxDiff()
        : maxScalarId(this->getStatistics().subscribeMax())
    {

    }

    virtual plb::BlockDomain::DomainT appliesTo() const override
    {
        return plb::BlockDomain::bulk;
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        modified[0] = plb::modif::staticVariables;
        modified[1] = plb::modif::staticVariables;
    }

    virtual ComputeMaxDiff<T,dim>* clone() const override
    {
        return new ComputeMaxDiff<T,dim>(*this);
    }


private:
    using BaseClass =
    typename PlbTypeDeduction::GetReductiveBoxProcessingFunctionalXD_SS<T,T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using Field1 = typename Traits::template argument<2>::type;
    using Field2 = typename Traits::template argument<3>::type;

public:

    virtual void process(Box domain, Field1& old, Field2& curr) override
    {
        plb::BlockStatistics& statistics = this->getStatistics();

        const auto offsetCurr = plb::computeRelativeDisplacement(old, curr);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& pos = *it;
            const T& oldVal = DataAccess::get(old, pos);
            const T& currVal = DataAccess::get(curr, pos + offsetCurr);
            const T diff = currVal - oldVal;
            const T relErr = diff / oldVal;
            const T absErr = std::abs(relErr);
            statistics.gatherMax(maxScalarId, absErr);
        }
    }

    T getMax() const
    {
        const plb::BlockStatistics& statistics = this->getStatistics();
        const T doubleMax = statistics.getMax(maxScalarId);
        if (std::numeric_limits<T>::is_integer)
        {
            return (T) plb::util::roundToInt(doubleMax);
        }
        return (T) doubleMax;
    }

private:
    plb::plint maxScalarId;
};

template<typename T, size_t dim>
class PalabosConvergenceCheck
{
public:
    explicit PalabosConvergenceCheck(const T& tolerance, bool& isConverged)
        : m_tolerance(tolerance)
        , m_isConverged(isConverged)
        , m_oldTime(0.0)
    {

    }

    void init(const size_t nFields)
    {
        m_oldData.resize(nFields, 0.0);
        m_oldTime = 0.0;
    }

    void setOldData(const size_t i, const T& data, const double& time)
    {
        assert(m_oldData.size() > i);
        m_oldData[i] = data;
        m_oldTime = time;
    }

    bool isConverged(const size_t i, const T& curr, const double& currTime)
    {
        assert(m_oldData.size() > i);

        const T& old = m_oldData[i];
        const T diff = curr - old;
        const T relErr = diff / old;
        const T maxErr = std::abs(relErr);

        const double dt = currTime - m_oldTime;
        const T errPerTime = maxErr / dt;

        IPP::pcout << "convergenceErr: " << errPerTime << std::endl;

        if(std::abs(errPerTime) <= m_tolerance)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    void setConvergence(const bool isConverged)
    {
        m_isConverged = isConverged;
    }

private:
    const double& m_tolerance;
    bool& m_isConverged;
    std::vector<T> m_oldData;
    double m_oldTime;
};

}

#endif // PALABOSCONVERGENCECHECK_H
