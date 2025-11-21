#ifndef TRANSPORTTYPES_H
#define TRANSPORTTYPES_H

#ifdef IPP_MULTI_PRECISION
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

namespace IPP
{

#ifdef IPP_MULTI_PRECISION
using TransportFloatingPointType = boost::multiprecision::cpp_dec_float_50;
#else
using TransportFloatingPointType = double;
#endif

}

#endif // TRANSPORTTYPES_H
