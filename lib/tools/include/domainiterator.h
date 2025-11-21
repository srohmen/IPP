#ifndef DOMAINITERATOR_H
#define DOMAINITERATOR_H

#include <cstddef>

namespace IPP
{

namespace DataAccess
{

template<size_t dim>
class DomainIterator;

//template<typename Box>
//auto begin(const Box& domain);

//template<typename Box>
//auto end(const Box& domain);

}
}


#include "domainiterator2d.h"
#include "domainiterator3d.h"

#endif // DOMAINITERATOR_H
