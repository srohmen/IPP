#ifndef CELLIDCONVERTTOSIFUNC_H
#define CELLIDCONVERTTOSIFUNC_H

#include "cellidconvertfuncwrapper.h"
#include "cellid_si_func.h"

#include "abstractgeometricsicalc.h"

namespace IPP
{

class CellIDConvertToSIFunc : public CellIDConvertFuncWrapper<CellID_SI_Func, AbstractGeometricSICalc>
{
public:
    CellIDConvertToSIFunc(const ArrayDimensionConvert* indexConv, const AbstractGeometricSICalc& siCalc);
};

} // end of namespace

#endif // CELLIDCONVERTTOSIFUNC_H
