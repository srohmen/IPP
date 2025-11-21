#ifndef UPDATEDISTANCEFIELD_H
#define UPDATEDISTANCEFIELD_H

#include "mpitools.h"
#include "simulationexchangedata.h"
#include "multitoserialdatasync.h"
#include "serialtomultidatasync.h"
#include "array_view.h"
#include "calcdistance.h"

namespace IPP
{

namespace UpdateDistanceField
{

template<size_t dim>
static void execute(SimulationExchangeData& simData, const double& porosThreshold)
{
    MPI_CHECK_SYNC;

    // currently the dist field algorithm is currently only working in global mode

    const std::vector<double>& porosLocal = simData.getPorosity();
    const FieldDecomposition& decomp = *simData.getDecomp();

    std::vector<double> porosGlobal;
    MultiToSerialDataSync<double> localToGlobal(porosGlobal, decomp);
    localToGlobal.pullToRoot(porosLocal);

    MPI_CHECK_SYNC;


    std::vector<double> distFieldGlobal;

    if(plb::global::mpi().isMainProcessor())
    {
        const IPPVector3DLong size = decomp.getGlobalSize();
        const av::bounds<dim> bounds = av::MakeBounds<dim>::get(size[0], size[1], size[2]);

        using ArrView = av::array_view<double, dim>;
        const ArrView porosityView(porosGlobal, bounds);

        distFieldGlobal.resize(decomp.getGlobalNumberCells());
        ArrView resultView(distFieldGlobal, bounds);
        CalcDistance::calcDistanceField<dim>(porosityView, porosThreshold, resultView);
    }

    MPI_CHECK_SYNC;

    std::vector<double>& distFieldLocal = simData.getDistField();

    SerialToMultiDataSync<double> globalToLocal(distFieldLocal, decomp);
    globalToLocal.pushFromRoot(distFieldGlobal);


    MPI_CHECK_SYNC;
}

}

}
#endif // UPDATEDISTANCEFIELD_H
