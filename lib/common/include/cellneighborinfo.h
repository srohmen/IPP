#ifndef CELLNEIGHBORINFO_H
#define CELLNEIGHBORINFO_H

#include <vector>

namespace IPP
{

struct CellNeighborInfo
{

    enum Location
    {
        L_Solid     = 0,
        L_Surface   = 1,
        L_Edge      = 2, // corresponding to "Corner" in 2D
        L_Corner    = 3,
        L_Bulk      = 4
    };

    enum Shape
    {
        // Shape values are shifted by 4 bytes because must fit into same byte as Location
        // no need of power 2 complements since mutual exclusive
        S_Filled    = (1 << 4),
        S_Flat      = (2 << 4),
        S_Concave   = (3 << 4),
        S_Convex    = (4 << 4)
    };


    CellNeighborInfo()
        : location(L_Bulk)
        , convexVolume(1.0)
        , shape(S_Flat)
        , noPrecip(false)
    {

    }


    Location location;
    double convexVolume;
    Shape shape;

    bool noPrecip;

    struct VolumeFractions
    {
        double self;
        double rest;
    };

    using NuclPhaseVolumeFractions = std::vector<VolumeFractions>;
    NuclPhaseVolumeFractions phaseVolFrac;


};

}

#endif // CELLNEIGHBORINFO_H
