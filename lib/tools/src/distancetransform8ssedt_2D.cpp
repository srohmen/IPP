#include "distancetransform8ssedt_2D.h"


#include <iostream>


#include "abstractscalarfield.h"
#include "deriv2d.h"
#include "rawscalarfield.h"
#include <cmath>
#include <assert.h>
#include <limits>

#include "bench_tools.h"

namespace IPP
{

namespace
{

typedef long long IntegerType;
using Deriv = Deriv2D<IntegerType>;
typedef RawScalarField2D<Deriv> ScalarFieldDeriv;

}

inline Deriv findHigher(const ScalarFieldDeriv &grid,
                        const Deriv &point,
                        int x, int y,
                        int offsetX, int offsetY )
{
    Deriv other = grid.getClamped(x+offsetX, y+offsetY);
    other.dx += offsetX;
    other.dy += offsetY;

    if (other.length2() < point.length2())
    {
        return other;
    }
    else
    {
        return point;
    }

}

inline void checkMask1(const size_t x, const size_t y, ScalarFieldDeriv& grid)
{
    Deriv point = grid.getClamped(x, y );
    point = findHigher( grid, point, x, y, -1,  0 );
    point = findHigher( grid, point, x, y,  0, -1 );
    point = findHigher( grid, point, x, y, -1, -1 );
    point = findHigher( grid, point, x, y,  1, -1 );
    grid[x][y] = point;
}

inline void checkMask2(const size_t x, const size_t y, ScalarFieldDeriv& grid)
{
    Deriv p = grid.getClamped(x, y);
    p = findHigher( grid, p, x, y, 1, 0 );
    grid[x][y] = p;
}

inline void passLowToHigh(ScalarFieldDeriv& grid)
{
    for (size_t y = 0; y < grid.getNy(); ++y)
    {
        for (size_t x = 0; x < grid.getNx(); ++x)
        {
            checkMask1(x, y, grid);
        }

        for (size_t x = grid.getNx(); x --> 0 ;)
        {
            checkMask2(x, y, grid);
        }
    }
}

inline void checkMask3(const size_t x, const size_t y, ScalarFieldDeriv& grid)
{
    Deriv point = grid.getClamped(x, y );
    point = findHigher( grid, point, x, y,  1,  0 );
    point = findHigher( grid, point, x, y,  0,  1 );
    point = findHigher( grid, point, x, y, -1,  1 );
    point = findHigher( grid, point, x, y,  1,  1 );
    grid[x][y] = point;
}

inline void checkMask4(size_t x, size_t y, ScalarFieldDeriv& grid)
{
    Deriv p = grid.getClamped(x, y);
    findHigher( grid, p, x, y, -1, 0 );
    grid[x][y] = p;
}

inline void passHighToLow(ScalarFieldDeriv& grid)
{
    for (size_t y = grid.getNy(); y --> 0 ;)
    {
        for (size_t x = grid.getNx(); x --> 0 ;)
        {
            checkMask3(x, y, grid);
        }

        for (size_t x = 0; x < grid.getNx(); x++)
        {
            checkMask4(x, y, grid);
        }
    }
}

inline void generateSDF(ScalarFieldDeriv& grid)
{
    passLowToHigh(grid);
    passHighToLow(grid);
}

void runTransform(const AbstractScalarField& input, const double& threshold,
                  AbstractScalarField& output)
{
    const size_t nx = input.getNx();
    const size_t ny = input.getNy();

    ScalarFieldDeriv grid1(nx, ny);
    ScalarFieldDeriv grid2(nx, ny);

    // Initialize the grids
    for( size_t y = 0; y < ny; ++y )
    {
        for ( size_t x = 0; x < nx; ++x )
        {
            // z always zero for 2D field
            const size_t z = 0;
            const double val = input.get(x, y, z);

            // Points inside get marked with a dx/dy of zero.
            // Points outside get marked with an infinitely large distance.
            const bool isInside = val < threshold;
            if(isInside)
            {
                grid1[x][y] = Deriv(0, 0);
                grid2[x][y] = Deriv();
            }
            else
            {
                grid2[x][y] = Deriv(0, 0);
                grid1[x][y] = Deriv();
            }
        }

    }

    // Generate the SDF.
    generateSDF( grid1 );
    generateSDF( grid2 );

    for( size_t y = 0; y < ny; ++y )
    {
        for ( size_t x = 0; x < nx; ++x )
        {
            const size_t z = 0;

            // Calculate the actual distance from the dx/dy
            const Deriv& p1 = grid1.getClamped(x, y);
            const Deriv& p2 = grid2.getClamped(x, y);
            double dist1 = std::sqrt( (double)p1.length2() );
            double dist2 = std::sqrt( (double)p2.length2() );
            double dist = dist1 - dist2;

            output.set(x, y, z, dist);
        }
    }
}

DistanceTransform8SSEDT_2D::DistanceTransform8SSEDT_2D()
{

}

void DistanceTransform8SSEDT_2D::calc(const AbstractScalarField& input, const double& threshold,
                                      AbstractScalarField& output)
{
    BENCH(runTransform(input, threshold, output));
}


} // end of namespace
