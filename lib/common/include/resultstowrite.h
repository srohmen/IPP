#ifndef RESULTSTOWRITE_H
#define RESULTSTOWRITE_H

#include <unordered_set>

namespace IPP
{

enum ResultType
{
    RT_ComponentConcentration,
    RT_ComponentFlux,
    RT_FluidDensityAndFlux,
    RT_Phase,
    RT_Porosity,
    RT_CapillaryPorosity,
    RT_InertFrac,
    RT_DiffusionCoef,
    RT_DiffusiveTransportScalar,
    RT_HydrodynamicTransportScalar,
    RT_DistanceField,
    RT_EnabledCells,

    // LBM debugging output
    RT_LBM_PostTransportConc,
    RT_LBM_Source,
    RT_LBM_Velocity,
    RT_LBM_Gradient,
    RT_LBM_TransPoros,

    RT_Budget,
    RT_PermeabilityInfo,

    RT_MaxRT
};



class ResultsToWrite
{
public:
    ResultsToWrite()
    {

    }

    void add(const ResultType& key)
    {
        m_data.insert(key);
    }

    void remove(const ResultType& key)
    {
        m_data.erase(key);
    }

    bool contains(const ResultType& key) const
    {
        return m_data.find(key) != m_data.end();
    }

private:
    struct ResultTypeHash
    {
        size_t operator()(const ResultType& key) const
        {
            return static_cast<size_t>(key);
        }
    };

    std::unordered_set<ResultType, ResultTypeHash> m_data;
};


}

#endif // RESULTSTOWRITE_H
