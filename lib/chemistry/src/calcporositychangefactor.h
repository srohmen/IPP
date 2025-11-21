#ifndef CALCPOROSITYCHANGEFACTOR_H
#define CALCPOROSITYCHANGEFACTOR_H


namespace CalcPorosityChangeFactor
{

inline double calc(const double& oldPorosity,
                   const double& currPorosity)
{
    if(currPorosity == 0.0)
    {
        return 0.0;
    }
    else if(oldPorosity == currPorosity)
    {
        return 1.0;
    }
    else
    {
        return oldPorosity / currPorosity;
    }
}

}

#endif // CALCPOROSITYCHANGEFACTOR_H
