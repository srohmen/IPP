#ifndef GINZBURGTRT_H
#define GINZBURGTRT_H

#include <palabos/core/array.h>

namespace IPP
{

namespace GinzburgTRT
{

template<typename T, typename Descriptor>
T equilibrium(const plb::plint iPop,
              const T& cPhi,
              const T& rhoBar,
              const plb::Array<T,Descriptor::d>& j,
              const T& jSqr,
              const T& poros)
{
    assert(iPop > 0);

    T jd = 0;
    for(plb::plint d = 0; d < Descriptor::d; ++d)
    {
        jd += Descriptor::c[iPop][d] * j[d];
    }

    const T rhoBarPhi = rhoBar * cPhi;
    const T w = 0.5;

    // TODO: check "In"-factor for lattice type
    // const T eq = w * (rhoBarPhi + poros * jd);
    const T eq = w * (rhoBarPhi + poros * jd + poros * jSqr);

    return eq;

}

template<typename T, typename Descriptor>
void equilibrium(const T& cPhi,
                 const T& rhoBar,
                 const plb::Array<T,Descriptor::d>& j,
                 const T& jSqr,
                 const T& poros,
                 plb::Array<T,Descriptor::q>& eq)
{
    T fSum = 0;

    for(plb::plint iPop = 1; iPop < Descriptor::q; ++iPop)
    {
        eq[iPop] = equilibrium<T,Descriptor>(iPop, cPhi, rhoBar, j, jSqr, poros);
        fSum += eq[iPop];
    }
    eq[0] = poros * rhoBar - fSum;
}

}

}

#endif // GINZBURGTRT_H
