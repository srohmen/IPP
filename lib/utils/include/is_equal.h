#ifndef IS_EQUAL_H
#define IS_EQUAL_H

#include <type_traits>


template < std::size_t M, std::size_t N >
struct is_equal : std::false_type
{
  typedef std::false_type type;
};

template < std::size_t N >
struct is_equal < N, N > : std::true_type
{
  typedef std::true_type type;
};

#endif // IS_EQUAL_H
