#ifndef PLBCHECKDATA_H
#define PLBCHECKDATA_H

#include "checknan.h"
#include "checkneg.h"

namespace IPP
{


template<typename Scalar, typename Field>
inline bool hasNoNeg(Field& field)
{
    bool hasNeg = false;
    CheckNeg<Scalar> check(hasNeg);
    plb::apply(check, field);

    return hasNeg == false;
}

template<typename Scalar, typename Field>
inline bool hasNoNan(Field& field)
{
    bool hasNaN = false;
    CheckNaN<Scalar> check(hasNaN);
    plb::apply(check, field);

    return hasNaN == false;
}


}

#endif // PLBCHECKDATA_H
