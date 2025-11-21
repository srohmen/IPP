#ifndef CSHGELDIFFUSIONCOEFCALC_H
#define CSHGELDIFFUSIONCOEFCALC_H

namespace IPP
{

class CSHGelDiffusionCoefCalc
{
public:
    CSHGelDiffusionCoefCalc();

    static double calc(const double &Dgel, const double &gelPoros);
};

}

#endif // CSHGELDIFFUSIONCOEFCALC_H
