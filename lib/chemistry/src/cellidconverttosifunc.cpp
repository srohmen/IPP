#include "cellidconverttosifunc.h"
#include "abstractgeometricsicalc.h"

namespace IPP
{

CellIDConvertToSIFunc::CellIDConvertToSIFunc(const ArrayDimensionConvert *indexConv, const AbstractGeometricSICalc &siCalc)
    :  CellIDConvertFuncWrapper<CellID_SI_Func, AbstractGeometricSICalc>(indexConv, siCalc)
{

}


} // end of namespace LBGeoChem
