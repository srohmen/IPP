#ifndef GOLDENSECTIONSEARCH_H
#define GOLDENSECTIONSEARCH_H

#include <vector>
#include <cstddef>
#include <assert.h>
#include <cmath>

namespace IPP
{

namespace GoldenSectionSearch
{

// scalar version
bool isConverged(const double& x0,
                 const double& x1,
                 const double& tolerance);

template<typename Func>
void minimum( Func& func, const double& tolerance,
              double& a, double& fa,
              double& b, double& fb)
{
    assert(tolerance > 0.0);

    static const double sqrt5 = std::sqrt(5); // 2.236067977499789696;
    static const double lambda = 0.5 * (sqrt5 - 1.0);
    static const double mu = 0.5 * (3.0 - sqrt5);  // = 1 - lambda
    double x1;
    double x2;
    double fx1;
    double fx2;


    // finding  first two internal points
    x1 = b - lambda * (b - a);
    x2 = a + lambda * (b - a);
    fx1 = func(x1);
    fx2 = func(x2);



    // finfing minumum in interval and reducing interval size gradually
    while ( isConverged( a, b, tolerance) == false)
    {
        if (fx1 > fx2)
        {
            a = x1;
            fa = fx1;

            if ( isConverged( a, b, tolerance) )
            {
                break;
            }

            x1 = x2;
            fx1 = fx2;
            x2 = b - mu * (b - a);
            fx2= func(x2);
        }
        else
        {
            b = x2;
            fb = fx2;

            if ( isConverged( a, b, tolerance) )
            {
                break;
            }

            x2 = x1;
            fx2 = fx1;
            x1 = a + mu * (b - a);
            fx1 = func(x1);
        }
    }
}


// vector version
bool isConverged(const std::vector<double>& x0Vec,
                 const std::vector<double>& x1Vec,
                 const double& tolerance);

bool isConverged(std::vector<size_t>& toCheck,
                 const std::vector<double>& x0Vec,
                 const std::vector<double>& x1Vec,
                 const double& tolerance);

template<typename Func>
void minimum( Func& func, const double& tolerance,
              std::vector<double>& a, std::vector<double>& fa,
              std::vector<double>& b, std::vector<double>& fb,
              std::vector<size_t> toOptimize)
{
    assert(tolerance > 0.0);

    const size_t n = a.size();

    static const double sqrt5 = std::sqrt(5); // 2.236067977499789696;
    static const double lambda = 0.5 * (sqrt5 - 1.0);
    static const double mu = 0.5 * (3.0 - sqrt5);  // = 1 - lambda
    std::vector<double> x1 = a;
    std::vector<double> x2 = b;
    std::vector<double> fx1 = fa;
    std::vector<double> fx2 = fb;

    func.setEnabledCells(toOptimize);


    // finding  first two internal points
    for(size_t i = 0; i < toOptimize.size(); ++i)
    {
        const size_t iCell = toOptimize[i];
        x1[iCell] = b[iCell] - lambda * (b[iCell] - a[iCell]);
        x2[iCell] = a[iCell] + lambda * (b[iCell] - a[iCell]);
    }
    func(x1, fx1);
    func(x2, fx2);



    // finding minumum in interval and reducing interval size gradually
    while ( isConverged(toOptimize, a, b, tolerance) == false)
    {
        if (fx1 > fx2)
        {
            a = x1;
            fa = fx1;

            if ( isConverged(toOptimize, a, b, tolerance) )
            {
                break;
            }

            x1 = x2;
            fx1 = fx2;
            for(size_t i = 0; i < n; ++i)
            {
                x2[i] = b[i] - mu * (b[i] - a[i]);
            }
            func(x2, fx2);
        }
        else
        {
            b = x2;
            fb = fx2;

            if ( isConverged(toOptimize, a, b, tolerance) )
            {
                break;
            }

            x2 = x1;
            fx2 = fx1;
            for(size_t i = 0; i < n; ++i)
            {
                x1[i] = a[i] + mu * (b[i] - a[i]);
            }
            func(x1, fx1);
        }
    }
}

}

}

#endif // GOLDENSECTIONSEARCH_H
