#include "findzeroporosityroot.h"

#include <unordered_map>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "porositycalcutils.h"
#include "retrievephaseamounts.h"
#include "localaccessphreeqcrm.h"
#include "indexhelper.h"
#include "geometrytools.h"
#include "ippexception.h"

namespace IPP
{

namespace FindZeroPorosityRoot
{

SetTotalsAndRun::SetTotalsAndRun(LocalAccessPhreeqcRM& phreeqc,
                                 const std::vector<double>& totalDiff,
                                 const double& porosity, const size_t iCell)
    : m_phreeqc(phreeqc)
    , m_totalDiff(totalDiff)
    , m_porosity(porosity)
    , m_iCell(iCell)
    , m_currInterp(1.0)
{

}

void SetTotalsAndRun::setAndRun(const double& t)
{
    const size_t nComps = m_totalDiff.size();
    const size_t nCells = m_phreeqc.getNumberCells();

    std::vector<double> concVec;
    m_phreeqc.getConcentrationsLocal(concVec);


    // t = 1 -> moleDiff = (1 - t) * -totalDiff = 0
    // t = 0 -> moleDiff = (1 - t) * -totalDiff = -totalDiff
    const double interpDiff = t - m_currInterp;
    m_currInterp = t;

    for(size_t iComp = 0; iComp < nComps; ++iComp)
    {
        const double& totalDiff = m_totalDiff[iComp];
        const double moleDiff = interpDiff * totalDiff;

        const size_t index = IndexHelper::toPhreeqcRawIndex(nCells, m_iCell, iComp);

        const double concDiffNormalized = moleDiff / m_porosity;

        assert(index < concVec.size());
        double& conc = concVec[index];
        assert(iComp == 3 || conc >= 0.0);


        // new conc could be negative if e.g. a large
        // amount of H was added in last iteration
        conc += concDiffNormalized;
        // assert(iComp == 3 || newConc >= 0.0);
    }


    m_phreeqc.setConcentrationsLocal(concVec);
    m_phreeqc.runCellsLocal();


}

class FindPorosZero
{

public:
    FindPorosZero(SetTotalsAndRun& setTotals,
                  const RetrievePhaseAmounts& retrievePhases,
                  const PorosityCalc& porosCalc,
                  const double& porosCacheA,
                  const double& porosCacheB,
                  const size_t iCell,
                  const double& porosOff)
        : m_setTotals(setTotals)
        , m_retrievePhases(retrievePhases)
        , m_porosCalc(porosCalc)
        , m_iCell(iCell)
        , m_porosOffset(porosOff)
    {
        // TODO: enable optim with poros result caching
        // m_cache[0.0] = m_porosCacheA;
        // m_cache[1.0] = m_porosCacheB;
    }

    double operator()(const double& t)
    {
        const double poros = calcCached(t);
        return poros - m_porosOffset;
    }

    double calcCached(const double& t)
    {
        const auto it = m_cache.find(t);
        if(it == m_cache.end())
        {
            const double poros = runCalcPorosity(t);
            m_cache[t] = poros;
            return poros;
        }
        else
        {
            return it->second;
        }
    }

    double getCached(const double& t) const
    {
        const auto it = m_cache.find(t);
        if(it == m_cache.end())
        {
            throw std::runtime_error("error in poros caching");
            return -1.0;
        }
        else
        {
            return it->second;
        }
    }

    double runCalcPorosity(const double& t) const
    {
        m_setTotals.setAndRun(t);

        const size_t nCells = m_porosCalc.getnDomains();
        const size_t nPhases = m_porosCalc.getnPhases();

        std::vector<double> precipMol(nCells * nPhases);
        m_retrievePhases.retrieve(precipMol);

        GeometryTools::transpose(precipMol, nCells, nPhases);

        const PorosityCalc::ConstCellPhasesView precipMolView(precipMol, { (ptrdiff_t)nCells,
                                                                           (ptrdiff_t)nPhases });
        std::vector<double> porosity;
        m_porosCalc.calcPorosity(precipMolView, porosity);

        const double& poros = porosity[m_iCell];

        std::cout << "expensive: " << t << "\t" << poros << std::endl;

        return poros;
    }

    SetTotalsAndRun& m_setTotals;
    const RetrievePhaseAmounts& m_retrievePhases;
    const PorosityCalc& m_porosCalc;
    const size_t m_iCell;
    const double m_porosOffset;

    std::unordered_map<double, double> m_cache;

};


static double funcWrapper (double x, void *funcRaw)
{
    FindPorosZero *func = static_cast<FindPorosZero*>(funcRaw);
    return (*func)(x);
}

static double findRoot(const FindPorosZero& porosFunc, gsl_root_fsolver* solver)
{
    int status;
    double xResult;
    do
    {
        status = gsl_root_fsolver_iterate(solver);
        xResult = gsl_root_fsolver_root(solver);
        const double x0 = gsl_root_fsolver_x_lower(solver);
        const double x1 = gsl_root_fsolver_x_upper(solver);


        status = gsl_root_test_interval(x0, x1, 1.0E-9, 1.E-3);

        if(status == GSL_CONTINUE)
        {
            const double porosX0 = porosFunc.getCached(x0);
            const double porosX1 = porosFunc.getCached(x1);

            if(porosX1 < porosX0)
            {
                if(porosX1 > 0.0 && porosX1 <= 1.0E-6)
                {
                    xResult = x1;
                    status = GSL_SUCCESS;
                }
            }
            else
            {
                if(porosX0 > 0.0 && porosX0 <= 1.0E-6)
                {
                    xResult = x0;
                    status = GSL_SUCCESS;
                }
            }
        }
        else if(status == GSL_SUCCESS)
        {
            // special case: negative porosity is very close to zero / within tolerance
            // take the positive value

            const double porosX0 = porosFunc.getCached(x0);
            const double porosX1 = porosFunc.getCached(x1);


            if(porosX0 < 0.0)
            {
                IPPCheck::assertCheck(porosX1 >= 0.0);
                xResult = x1;
            }
            else if(porosX1 < 0.0)
            {
                IPPCheck::assertCheck(porosX0 >= 0.0);
                xResult = x0;
            }
            else
            {
                if(porosX0 >= porosX1)
                {
                    xResult = x0;
                }
                else if(porosX1 > porosX0)
                {
                    xResult = x1;
                }
                else
                {
                    throw std::runtime_error("error in converging zero porosity: "
                                             + std::to_string(porosX0) + " <-> "
                                             + std::to_string(porosX1));
                }
            }


        }
    }
    while(status == GSL_CONTINUE);


    return xResult;
}

Result findRoot(SetTotalsAndRun& setAndRun,
                const RetrievePhaseAmounts& retrievePhases,
                const PorosityCalc& porosCalc,
                const double& oldPoros,
                const double& currPoros,
                const size_t iCell)
{
    FindPorosZero porosFunc(setAndRun, retrievePhases, porosCalc,
                            oldPoros, currPoros, iCell, 0.0);

#ifndef NDEBUG
    const double one = porosFunc(1.0);
    const double zero = porosFunc(0.0);

    assert(std::abs(one - currPoros) < 1.0E-2);
    assert(std::abs(zero - oldPoros) < 1.0E-2);
#endif


    gsl_function func;
    func.function = &funcWrapper;
    func.params = &porosFunc;

    const gsl_root_fsolver_type* solverType = gsl_root_fsolver_brent;
    gsl_root_fsolver* solver = gsl_root_fsolver_alloc (solverType);

    double x0 = 0.0;
    double x1 = 1.0;
    gsl_root_fsolver_set (solver, &func, x0, x1);


    const double xResult = findRoot(porosFunc, solver);


    gsl_root_fsolver_free (solver);


    // expensive recalc because the internal state of phreeqc must be set to result
    const double resultPoros = porosFunc.runCalcPorosity(xResult);
    std::cout << "result: " << xResult << "\t" << resultPoros << std::endl;

    IPPCheck::assertCheck(resultPoros >= 0.0);


    return { xResult, resultPoros } ;
}




}

}


