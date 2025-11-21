#include "distancetransform8ssedt_3D.h"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "deriv3d.h"
#include "abstractscalarfield.h"
#include "rawscalarfield.h"


#include "bench_tools.h"


namespace IPP
{


namespace
{



typedef int IntegerType;
using Deriv = Deriv3D<IntegerType>;
typedef RawScalarField3D<Deriv> ScalarFieldDeriv;


inline void findHigher(const ScalarFieldDeriv &grid,
                       Deriv &point,
                       int x, int y, int z,
                       int offsetX, int offsetY, int offsetZ)
{
    Deriv other = grid.getClamped(x+offsetX, y+offsetY, z+offsetZ);

    const Deriv offset(offsetX, offsetY, offsetZ);
    other += offset;

    if (other.length2() < point.length2())
    {
        point = other;
    }

}

// for mask check Ragnemalm (1993)

inline void checkMask1(const size_t x, const size_t y, const size_t z,
                       ScalarFieldDeriv& grid)
{
    Deriv point = grid.getClamped(x, y, z);

    findHigher( grid, point, x, y, z, -1, -1, -1);
    findHigher( grid, point, x, y, z, -1, -1,  0);
    findHigher( grid, point, x, y, z, -1,  0, -1);
    findHigher( grid, point, x, y, z, -1,  0,  0);
    findHigher( grid, point, x, y, z, -1,  1, -1);
    findHigher( grid, point, x, y, z, -1,  1,  0);

    findHigher( grid, point, x, y, z,  0, -1, -1);
    findHigher( grid, point, x, y, z,  0,  0, -1);
    findHigher( grid, point, x, y, z,  0,  1, -1);

    findHigher( grid, point, x, y, z,  1, -1, -1);
    findHigher( grid, point, x, y, z,  1,  0, -1);
    findHigher( grid, point, x, y, z,  1,  1, -1);

    grid[x][y][z] = point;
}

inline void checkMask2(const size_t x, const size_t y, const size_t z,
                       ScalarFieldDeriv& grid)
{
    Deriv point = grid.getClamped(x, y, z);

    findHigher( grid, point, x, y, z, -1, -1, -1);
    findHigher( grid, point, x, y, z, -1, -1,  0);
    findHigher( grid, point, x, y, z, -1, -1,  1);

    findHigher( grid, point, x, y, z,  0, -1, -1);
    findHigher( grid, point, x, y, z,  0, -1,  0);
    findHigher( grid, point, x, y, z,  0, -1,  1);

    findHigher( grid, point, x, y, z,  1, -1, -1);
    findHigher( grid, point, x, y, z,  1, -1,  0);
    findHigher( grid, point, x, y, z,  1, -1,  1);

    findHigher( grid, point, x, y, z,  1,  0, -1);
    findHigher( grid, point, x, y, z,  1,  0,  0);
    findHigher( grid, point, x, y, z,  1,  0,  1);

    grid[x][y][z] = point;
}

inline void checkMask3(const size_t x, const size_t y, const size_t z,
                       ScalarFieldDeriv& grid)
{
    Deriv point = grid.getClamped(x, y, z);

    findHigher( grid, point, x, y, z, -1, -1,  0);
    findHigher( grid, point, x, y, z, -1, -1,  1);
    findHigher( grid, point, x, y, z, -1,  0,  0);
    findHigher( grid, point, x, y, z, -1,  0,  1);
    findHigher( grid, point, x, y, z, -1,  1,  0);
    findHigher( grid, point, x, y, z, -1,  1,  1);

    findHigher( grid, point, x, y, z,  0, -1,  1);
    findHigher( grid, point, x, y, z,  0,  0,  1);
    findHigher( grid, point, x, y, z,  0,  1,  1);

    findHigher( grid, point, x, y, z,  1, -1,  1);
    findHigher( grid, point, x, y, z,  1,  0,  1);
    findHigher( grid, point, x, y, z,  1,  1,  1);

    grid[x][y][z] = point;
}

inline void checkMask4(const size_t x, const size_t y, const size_t z,
                       ScalarFieldDeriv& grid)
{
    Deriv point = grid.getClamped(x, y, z);

    findHigher( grid, point, x, y, z, -1,  1, -1);
    findHigher( grid, point, x, y, z, -1,  1,  0);
    findHigher( grid, point, x, y, z, -1,  1,  1);

    findHigher( grid, point, x, y, z,  0,  1, -1);
    findHigher( grid, point, x, y, z,  0,  1,  0);
    findHigher( grid, point, x, y, z,  0,  1,  1);

    findHigher( grid, point, x, y, z,  1,  0, -1);
    findHigher( grid, point, x, y, z,  1,  0,  0);
    findHigher( grid, point, x, y, z,  1,  0,  1);

    findHigher( grid, point, x, y, z,  1,  1, -1);
    findHigher( grid, point, x, y, z,  1,  1,  0);
    findHigher( grid, point, x, y, z,  1,  1,  1);

    grid[x][y][z] = point;
}

inline void passPosZ(ScalarFieldDeriv& grid)
{
    for (size_t y = 0; y < grid.getNy(); ++y)
    {
        for (size_t z = 0; z < grid.getNz(); ++z)
        {
            for (size_t x = 0; x < grid.getNx(); ++x)
            {
                checkMask1(x, y, z, grid);
            }
        }
    }
}

inline void passPosY(ScalarFieldDeriv& grid)
{
    for (size_t z = 0; z < grid.getNz(); ++z)
    {
        for (size_t y = 0; y < grid.getNy(); ++y)
        {
            for (size_t x = grid.getNx(); x --> 0 ;)
            {
                checkMask2(x, y, z, grid);
            }
        }
    }
}

inline void passNegZ(ScalarFieldDeriv& grid)
{
    for (size_t y = 0; y < grid.getNy(); ++y)
    {
        for (size_t z = grid.getNz(); z --> 0 ;)
        {
            for (size_t x = 0; x < grid.getNx(); ++x)
            {
                checkMask3(x, y, z, grid);
            }
        }
    }
}

inline void passNegY(ScalarFieldDeriv& grid)
{
    for (size_t z = 0; z < grid.getNz(); ++z)
    {
        for (size_t y = grid.getNy(); y --> 0 ;)
        {
            for (size_t x = grid.getNx(); x --> 0 ;)
            {
                checkMask4(x, y, z, grid);
            }
        }
    }
}

inline void generateSDF(ScalarFieldDeriv& grid)
{
    passPosZ(grid);
    passPosY(grid);
    passNegZ(grid);
    passNegY(grid);
}

template<typename Field>
void print(const Field& outputWrapper)
{
    for( size_t z = 0; z < outputWrapper.getNz(); ++z )
    {
        for( size_t y = 0; y < outputWrapper.getNy(); ++y )
        {
            for ( size_t x = 0; x < outputWrapper.getNx(); ++x )
            {
                const Deriv& p = outputWrapper.get(x, y, z);
                std::cout << "(" << p.dx << "," << p.dy << "," << p.dz << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}


static void runTransform(const AbstractScalarField& input,
                         const double& threshold,
                         AbstractScalarField& output)
{
    const size_t nx = input.getNx();
    const size_t ny = input.getNy();
    const size_t nz = input.getNz();


    ScalarFieldDeriv grid1(nx, ny, nz);
    ScalarFieldDeriv grid2(nx, ny, nz);

    // Initialize the grids
    for ( size_t x = 0; x < nx; ++x )
    {
        for( size_t y = 0; y < ny; ++y )
        {
            for ( size_t z = 0; z < nz; ++z )
            {
                const double val = input.get(x, y, z);

                // Points inside get marked with a dx/dy of zero.
                // Points outside get marked with an infinitely large distance.
                const bool isInside = val < threshold;
                if(isInside)
                {
                    grid1[x][y][z] = Deriv(0, 0, 0);
                    grid2[x][y][z] = Deriv();
                }
                else
                {
                    grid2[x][y][z] = Deriv(0, 0, 0);
                    grid1[x][y][z] = Deriv();
                }
            }
        }
    }


    generateSDF( grid1 );
    generateSDF( grid2 );


    for ( size_t x = 0; x < nx; ++x )
    {
        for( size_t y = 0; y < ny; ++y )
        {
            for ( size_t z = 0; z < nz; ++z )
            {
                const Deriv& p1 = grid1.getClamped(x, y, z);
                const Deriv& p2 = grid2.getClamped(x, y, z);
                const IntegerType l1 = p1.length2();
                const IntegerType l2 = p2.length2();
                double dist1 = std::sqrt( (double)l1 );
                double dist2 = std::sqrt( (double)l2 );
                double dist = dist1 - dist2;
                output.set(x, y, z, dist);
            }
        }
    }
}

}

DistanceTransform8SSEDT_3D::DistanceTransform8SSEDT_3D()
{

}

void DistanceTransform8SSEDT_3D::calc(const AbstractScalarField &input,
                                      const double &threshold,
                                      AbstractScalarField &output)
{
    BENCH(runTransform(input, threshold, output));
}


}

