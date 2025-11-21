#ifndef RAWSCALARFIELD_H
#define RAWSCALARFIELD_H

#include <algorithm>
#include <tuple>
#include <cstddef>
#include <assert.h>

namespace IPP
{


template<typename T>
struct RawScalarField2D
{
    RawScalarField2D()
        : nx(0),
          ny(0),
          field(nullptr),
          rawData(nullptr)
    {

    }

    RawScalarField2D(const RawScalarField2D<T>& other)
        : nx(other.nx),
          ny(other.ny)
    {
        allocate();
        std::copy(other.rawData, other.rawData + nx*ny, rawData);
    }

    RawScalarField2D& operator=(const RawScalarField2D& other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;

        // reuse storage when possible
        if(nx != other.nx || ny != other.ny)
        {
            deallocate();
            nx = other.nx;
            ny = other.ny;
            allocate();
        }
        std::copy(other.rawData, other.rawData + nx*ny, rawData);
        return *this;
    }

    RawScalarField2D(const size_t nx, const size_t ny)
        : nx(nx),
          ny(ny)
    {
        allocate();
    }

    RawScalarField2D(const size_t nx, const size_t ny, const T& ini)
        : nx(nx),
          ny(ny)
    {
        allocate();
        std::fill(rawData, rawData + nx*ny, ini);
    }

    ~RawScalarField2D()
    {
        deallocate();
    }

    inline size_t getNx() const
    {
        return nx;
    }

    inline size_t getNy() const
    {
        return ny;
    }

    inline T& get(const size_t x, const size_t y)
    {
        assert(x < nx);
        assert(y < ny);
        return field[x][y];
    }

    inline const T& get(const size_t x, const size_t y) const
    {
        assert(x < nx);
        assert(y < ny);
        return field[x][y];
    }

    typedef std::tuple<size_t, size_t> Coords;
    inline T& get(const Coords& coords)
    {
        return this->get(std::get<0>(coords), std::get<1>(coords));
    }

    inline const T& get(const Coords& coords) const
    {
        return this->get(std::get<0>(coords), std::get<1>(coords));
    }

    T& operator[](const Coords& coords)
    {
        return this->get(std::get<0>(coords), std::get<1>(coords));
    }

    const T& operator[](const Coords& coords) const
    {
        return this->get(std::get<0>(coords), std::get<1>(coords));
    }


    inline T* operator[](const size_t x)
    {
        assert(x < nx);
        return field[x];
    }

    inline const T* operator[](const size_t x) const
    {
        assert(x < nx);
        return field[x];
    }

    inline const T& getClamped(const int x, const int y) const
    {
        if ( x >= 0 && y >= 0 &&
             x < (int)nx && y < (int)ny)
            return field[x][y];
        else
            return empty;
    }

    inline T* data()
    {
        return rawData;
    }

    inline const T* data() const
    {
        return rawData;
    }

private:
    void allocate()
    {
        rawData = new T[nx * ny];
        field = new T*[nx];
        for (size_t x = 0; x < nx; ++x)
        {
            field[x] = rawData + ny * x;
        }
    }

    void deallocate()
    {
        delete[] field;
        field = nullptr;

        delete[] rawData;
        rawData = nullptr;
    }

    size_t nx, ny;
    T* rawData;
    T** field;

private:
    static const T empty;
};

template <typename T>
const T RawScalarField2D<T>::empty = T();

typedef RawScalarField2D<double> RawScalarField2DDouble;


template<typename T>
struct RawScalarField3D
{
    RawScalarField3D()
        : nx(0),
          ny(0),
          nz(0),
          field(nullptr),
          rawData(nullptr)
    {

    }

    RawScalarField3D(const RawScalarField3D<T>& other)
        : nx(other.nx),
          ny(other.ny),
          nz(other.nz)
    {
        allocate();
        std::copy(other.rawData, other.rawData + nx*ny*nz, rawData);
    }

    RawScalarField3D& operator=(const RawScalarField3D& other)
    {
        // check for self-assignment
        if(&other == this)
            return *this;

        // reuse storage when possible
        if(nx != other.nx || ny != other.ny || nz != other.nz)
        {
            deallocate();
            nx = other.nx;
            ny = other.ny;
            nz = other.nz;
            allocate();
        }
        std::copy(other.rawData, other.rawData + nx*ny*nz, rawData);
        return *this;
    }

    RawScalarField3D(const size_t nx, const size_t ny, const size_t nz)
        : nx(nx),
          ny(ny),
          nz(nz)
    {
        allocate();
    }

    RawScalarField3D(const size_t nx, const size_t ny, const size_t nz, const T& ini)
        : nx(nx),
          ny(ny),
          nz(nz)
    {
        allocate();
        std::fill(rawData, rawData + nx*ny*nz, ini);
    }

    ~RawScalarField3D()
    {
        deallocate();
    }

    inline size_t getNx() const
    {
        return nx;
    }

    inline size_t getNy() const
    {
        return ny;
    }

    inline size_t getNz() const
    {
        return nz;
    }


    inline T& get(const size_t x, const size_t y, const size_t z)
    {
        assert(x < nx);
        assert(y < ny);
        assert(z < nz);
        return field[x][y][z];
    }

    inline const T& get(const size_t x, const size_t y, const size_t z) const
    {
        assert(x < nx);
        assert(y < ny);
        assert(z < nz);
        return field[x][y][z];
    }

    typedef std::tuple<size_t, size_t, size_t> Coords;
    inline T& get(const Coords& coords)
    {
        return this->get(std::get<0>(coords), std::get<1>(coords), std::get<2>(coords));
    }

    inline const T& get(const Coords& coords) const
    {
        return this->get(std::get<0>(coords), std::get<1>(coords), std::get<2>(coords));
    }

    T& operator[](const Coords& coords)
    {
        return this->get(std::get<0>(coords), std::get<1>(coords), std::get<2>(coords));
    }

    const T& operator[](const Coords& coords) const
    {
        return this->get(std::get<0>(coords), std::get<1>(coords), std::get<2>(coords));
    }


    inline T** operator[](const size_t x)
    {
        assert(x < nx);
        return field[x];
    }

    inline const T** operator[](const size_t x) const
    {
        assert(x < nx);
        return field[x];
    }

    inline const T& getClamped(const int x, const int y, const int z) const
    {
        if ( x >= 0 && y >= 0 && z >= 0 &&
             x < (int)nx && y < (int)ny && z < (int)nz)
            return field[x][y][z];
        else
            return empty;
    }

    inline T* data()
    {
        return rawData;
    }

    inline const T* data() const
    {
        return rawData;
    }

private:
    void allocate()
    {
        rawData = new T[nx * ny * nz];
        field = new T**[nx];
        for (size_t x = 0; x < nx; ++x)
        {
            field[x] = new T*[ny];
            for (size_t y = 0; y < ny; ++y)
            {
                field[x][y] = rawData + nz * (y + ny * x);
            }
        }
    }

    void deallocate()
    {
        for (size_t x = 0; x < nx; ++x)
        {
            delete [] field[x];
        }

        delete[] field;
        field = nullptr;


        delete[] rawData;
        rawData = nullptr;
    }

    size_t nx, ny, nz;
    T* rawData;
    T*** field;

private:
    static const T empty;
};

template <typename T>
const T RawScalarField3D<T>::empty = T();

typedef RawScalarField3D<double> RawScalarField3DDouble;


} // end of namespace


#endif // RAWSCALARFIELD_H
