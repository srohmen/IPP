#ifndef GETPLBSIZE_H
#define GETPLBSIZE_H

#include <array>

namespace IPP
{

template<size_t dim>
struct GetPlbSize;

template<>
struct GetPlbSize<2>
{
    template<typename Container>
    static std::array<size_t, 3> get(const Container& container)
    {
        return { {(size_t)container.getNx(), (size_t)container.getNy(), 1} };
    }
};

template<>
struct GetPlbSize<3>
{
    template<typename Container>
    static std::array<size_t, 3> get(const Container& container)
    {
        return { {(size_t)container.getNx(), (size_t)container.getNy(), (size_t)container.getNz()} };
    }
};

}

#endif // GETPLBSIZE_H
