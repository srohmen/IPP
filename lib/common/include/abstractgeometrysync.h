#ifndef ABSTRACTGEOMETRYSYNC_H
#define ABSTRACTGEOMETRYSYNC_H

#include <vector>

#include "ippvector.h"
#include "cellneighborinfo.h"


namespace IPP
{

class AbstractGeometrySync
{
public:
    virtual ~AbstractGeometrySync() = default;


    typedef std::vector<IPPVector3DLong> NodeList;

    virtual void run(const NodeList& newPermNodes,
                     const NodeList& newInterfaceNodes,
                     const NodeList& newNonPermNonInterfaceNodes,
                     std::vector<double>& conc) = 0;

    virtual bool needDistField() const = 0;
    virtual bool needNeighMinPoros() const = 0;

    virtual void setPrintDebug(const bool printDebug) = 0;


    struct NeighborInfos
    {
        NeighborInfos(const double porosityThresh,
                      std::vector<CellNeighborInfo>& results)
            : porosityThresh(porosityThresh)
            , results(results)
        {

        }

        const double porosityThresh;
        std::vector<CellNeighborInfo>& results;
    };

    virtual NeighborInfos& getNeighInfos() = 0;
};


}

#endif // ABSTRACTGEOMETRYSYNC_H
