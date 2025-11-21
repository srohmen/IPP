#ifndef PHREEQCDEFAULTCONVERGENCEVALUES_H
#define PHREEQCDEFAULTCONVERGENCEVALUES_H

#include <cstddef>

namespace PhreeqcDefaultConvergenceValues
{

static const size_t iterations = 100;
static const double tolerance = 1.0E-15;
static const double convergence_tolerance = 1.0E-8;
static const bool diagonalScaling = false;

}


#endif // PHREEQCDEFAULTCONVERGENCEVALUES_H
