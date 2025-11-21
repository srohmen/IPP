#ifndef GETNEIGHINFO_H
#define GETNEIGHINFO_H

#include "cellneighborinfo.h"

#include <palabos/core/geometry2D.h>
#include <palabos/core/geometry3D.h>

#include "ippconstants.h"

namespace IPP
{

namespace GetNeighInfo
{

template<typename T>
static T calcTriangleSurface(const T& a, const T& b)
{
    const T res = a * b * 0.5;
    return res;
}


template<typename T, typename ScalarField, typename FlagField, size_t dim>
class GetNeighInfo;

template<typename T, typename ScalarField, typename FlagField>
class GetNeighInfo<T, ScalarField, FlagField, 2>
{
public:
    GetNeighInfo(const ScalarField& porosityField,
                 const T& thresh,
                 const FlagField& flagField)
        : m_porosityField(porosityField)
        , m_thresh(thresh)
        , m_flagField(flagField)
    {

    }

    CellNeighborInfo get(const plb::Dot2D& pos) const
    {
        // checking 4 combinations

        CellNeighborInfo result;

        if(this->isSolidOwn(pos))
        {
            result.location = CellNeighborInfo::L_Solid;
            result.shape = CellNeighborInfo::S_Filled;
            result.convexVolume = 0.0;
        }
        else
        {
            checkSurface(pos, result, -1,  0);
            checkSurface(pos, result,  1,  0);
            checkSurface(pos, result,  0, -1);
            checkSurface(pos, result,  0,  1);

            // skip edge check if position is already on surface
            if(result.location > CellNeighborInfo::L_Surface)
            {
                // not flat surface. check for edge and corner case

                checkEdge(pos, result, -1, -1);
                checkEdge(pos, result, -1, 1);
                checkEdge(pos, result, 1, -1);
                checkEdge(pos, result, 1, 1);
            }
        }

        return result;
    }

private:

    bool isSolidOwn(const plb::Dot2D& pos) const
    {
        const T& porosity = m_porosityField.get(pos.x, pos.y);
        return porosity < m_thresh;
    }

    bool isSolidNext(const plb::Dot2D& pos,
                     const plb::plint dx,
                     const plb::plint dy) const
    {
        const T& porosityNext = m_porosityField.get(pos.x + dx, pos.y + dy);
        const T& flagNext = m_flagField.get(pos.x + dx, pos.y + dy);
        const bool result = porosityNext < m_thresh || flagNext == IPPConstants::s_isSteadyFlag;
        return result;
    }

    void checkSurface(const plb::Dot2D& pos,
                      CellNeighborInfo& result,
                      const plb::plint dx,
                      const plb::plint dy) const
    {
        if(this->isSolidNext(pos, dx, dy))
        {
            if(result.location == CellNeighborInfo::L_Surface)
            {
                result.shape = CellNeighborInfo::S_Concave;
            }

            result.location = CellNeighborInfo::L_Surface;
            result.convexVolume = 0.0;
        }
    }

    void checkEdge(const plb::Dot2D& pos,
                   CellNeighborInfo& result,
                   const plb::plint dx,
                   const plb::plint dy) const
    {
        if(this->isSolidNext(pos, dx, dy))
        {
            result.location = std::min(result.location, CellNeighborInfo::L_Edge);
            result.shape = CellNeighborInfo::S_Convex;

            const T& a = m_porosityField.get(pos.x + dx, pos.y);
            const T& b = m_porosityField.get(pos.x, pos.y + dy);
            const T newVol = calcTriangleSurface<T>(1.0 - a, 1.0 - b);
            result.convexVolume -= static_cast<decltype(result.convexVolume)>(newVol);
        }

    }

    const ScalarField& m_porosityField;
    const T& m_thresh;
    const FlagField& m_flagField;
};



template<typename ScalarField, size_t dimToIgnore>
struct CalcEdgeNeighVol;

template<typename ScalarField>
struct CalcEdgeNeighVol<ScalarField, 0>
{
    static auto calcVol(const ScalarField& porosityField, const plb::Dot3D& pos,
                        const plb::plint /*dx*/, const plb::plint dy, const plb::plint dz)
    {
        const auto& a = porosityField.get(pos.x, pos.y,        pos.z + dz);
        const auto& b = porosityField.get(pos.x, pos.y + dy,   pos.z);
        using T = typename std::decay<decltype(a)>::type;
        const auto result = calcTriangleSurface<T>(1.0 - a, 1.0 - b);
        return result;
    }
};


template<typename ScalarField>
struct CalcEdgeNeighVol<ScalarField, 1>
{
    static auto calcVol(const ScalarField& porosityField, const plb::Dot3D& pos,
                        const plb::plint dx, const plb::plint /*dy*/, const plb::plint dz)
    {
        const auto& a = porosityField.get(pos.x,       pos.y, pos.z + dz);
        const auto& b = porosityField.get(pos.x + dx,  pos.y, pos.z);
        using T = typename std::decay<decltype(a)>::type;
        const auto result = calcTriangleSurface<T>(1.0 - a, 1.0 - b);
        return result;
    }
};


template<typename ScalarField>
struct CalcEdgeNeighVol<ScalarField, 2>
{
    static auto calcVol(const ScalarField& porosityField, const plb::Dot3D& pos,
                        const plb::plint dx, const plb::plint dy, const plb::plint /*dz*/)
    {
        const auto& a = porosityField.get(pos.x,       pos.y + dy, pos.z);
        const auto& b = porosityField.get(pos.x + dx,  pos.y,      pos.z);
        using T = typename std::decay<decltype(a)>::type;
        const auto result = calcTriangleSurface<T>(1.0 - a, 1.0 - b);
        return result;
    }
};


template<typename T, typename ScalarField, typename FlagField>
class GetNeighInfo<T, ScalarField, FlagField, 3>
{
public:
    GetNeighInfo(const ScalarField& porosityField,
                 const T& thresh,
                 const FlagField& flagField)
        : m_porosityField(porosityField)
        , m_thresh(thresh)
        , m_flagField(flagField)
    {

    }

    CellNeighborInfo get(const plb::Dot3D& pos) const
    {
        // checking 8 combinations

        CellNeighborInfo result;

        if(this->isSolidOwn(pos))
        {
            result.location = CellNeighborInfo::L_Solid;
            result.shape = CellNeighborInfo::S_Filled;
            result.convexVolume = 0.0;
        }
        else
        {
            checkSurface(pos, result, -1,  0,  0);
            checkSurface(pos, result,  1,  0,  0);
            checkSurface(pos, result,  0, -1,  0);
            checkSurface(pos, result,  0,  1,  0);
            checkSurface(pos, result,  0,  0, -1);
            checkSurface(pos, result,  0,  0,  1);

            // skip edge check if position is already on surface
            if(result.location > CellNeighborInfo::L_Surface)
            {
                // not flat surface. check for edge and corner case

                using S = ScalarField;
                checkEdge<CalcEdgeNeighVol<S,0>>(pos, result,  0, -1, -1);
                checkEdge<CalcEdgeNeighVol<S,0>>(pos, result,  0, -1,  1);
                checkEdge<CalcEdgeNeighVol<S,0>>(pos, result,  0,  1, -1);
                checkEdge<CalcEdgeNeighVol<S,0>>(pos, result,  0,  1,  1);

                checkEdge<CalcEdgeNeighVol<S,1>>(pos, result, -1,  0, -1);
                checkEdge<CalcEdgeNeighVol<S,1>>(pos, result, -1,  0,  1);
                checkEdge<CalcEdgeNeighVol<S,1>>(pos, result,  1,  0, -1);
                checkEdge<CalcEdgeNeighVol<S,1>>(pos, result,  1,  0,  1);

                checkEdge<CalcEdgeNeighVol<S,2>>(pos, result, -1, -1,  0);
                checkEdge<CalcEdgeNeighVol<S,2>>(pos, result, -1,  1,  0);
                checkEdge<CalcEdgeNeighVol<S,2>>(pos, result,  1, -1,  0);
                checkEdge<CalcEdgeNeighVol<S,2>>(pos, result,  1,  1,  0);

                // lower
                checkCorner(pos, result, -1, -1, -1);
                checkCorner(pos, result, -1,  1, -1);
                checkCorner(pos, result,  1, -1, -1);
                checkCorner(pos, result,  1,  1, -1);

                // upper
                checkCorner(pos, result, -1, -1,  1);
                checkCorner(pos, result, -1,  1,  1);
                checkCorner(pos, result,  1, -1,  1);
                checkCorner(pos, result,  1,  1,  1);
            }
        }

        return result;

    }


private:

    bool isSolidOwn(const plb::Dot3D& pos) const
    {
        const T& porosity = m_porosityField.get(pos.x, pos.y, pos.z);
        return porosity < m_thresh;
    }

    bool isSolidNext(const plb::Dot3D& pos,
                     const plb::plint dx,
                     const plb::plint dy,
                     const plb::plint dz) const
    {
        const T& porosityNext = m_porosityField.get(pos.x + dx, pos.y + dy, pos.z + dz);
        const T& flagNext = m_flagField.get(pos.x + dx, pos.y + dy, pos.z + dz);
        return porosityNext < m_thresh || flagNext == IPPConstants::s_isSteadyFlag;
    }

    void checkSurface(const plb::Dot3D& pos,
                      CellNeighborInfo& result,
                      const plb::plint dx, const plb::plint dy, const plb::plint dz) const
    {
        if(this->isSolidNext(pos, dx, dy, dz))
        {
            result.location = CellNeighborInfo::L_Surface;

            if(result.location == CellNeighborInfo::L_Surface)
            {
                result.shape = CellNeighborInfo::S_Concave;
            }

            result.convexVolume = 0.0;
        }
    }

    template<typename VolCalc>
    void checkEdge(const plb::Dot3D& pos, CellNeighborInfo& result,
                   const plb::plint dx, const plb::plint dy, const plb::plint dz) const
    {
        if(this->isSolidNext(pos, dx, dy, dz))
        {
            result.location = std::max(result.location, CellNeighborInfo::L_Edge);
            result.shape = CellNeighborInfo::S_Convex;
            const T newVol = VolCalc::calcVol(m_porosityField, pos, dx, dy, dz);
            result.convexVolume -= static_cast<decltype(result.convexVolume)>(newVol);
        }
    }

    void checkCorner(const plb::Dot3D& pos, CellNeighborInfo& result,
                     const plb::plint dx, const plb::plint dy, const plb::plint dz) const
    {
        if(this->isSolidNext(pos, dx, dy, dz))
        {
            result.location = std::min(result.location, CellNeighborInfo::L_Corner);
            result.shape = CellNeighborInfo::S_Convex;

            const T& a = m_porosityField.get(pos.x + dx,  pos.y,      pos.z);
            const T& b = m_porosityField.get(pos.x,       pos.y + dy, pos.z);
            const T& c = m_porosityField.get(pos.x,       pos.y,      pos.z + dz);
            const T newVol = calcPyramidVol(1.0 - a, 1.0 - b, 1.0 - c);
            result.convexVolume -= static_cast<decltype(result.convexVolume)>(newVol);
        }
    }

    static T calcPyramidVol(const T& x, const T& y, const T& z)
    {
        const T c =  1.0 / 6.0;
        return x * y * z * c;
    }



    const ScalarField& m_porosityField;
    const T& m_thresh;
    const FlagField& m_flagField;
};



}


} // end of namespace IPP


#endif // GETNEIGHINFO_H
