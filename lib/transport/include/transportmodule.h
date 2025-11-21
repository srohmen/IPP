#ifndef TRANSPORTMODULE_H
#define TRANSPORTMODULE_H

#include <string>
#include <vector>

#include <ipprendererfwd.h>
// #include "abstractscalarfieldfwd.h"
#include "domaininfos.h"

namespace IPP
{

class FieldDecomposition;
class NonLocalOperations;
class AbstractGeometrySync;
class ResultsToWrite;

class TransportModule
{
public:
    TransportModule();

    virtual ~TransportModule();

    virtual void init() = 0;

    virtual NonLocalOperations& getNonLocalOperations() = 0;
    virtual const DomainInfos& getDomainInfos() const = 0;
    virtual const FieldDecomposition* getDecomposition() const = 0;
    virtual void setGeomSync(AbstractGeometrySync* geomSync) = 0;

    // ordinary loop
    virtual void updatePostReactionState() = 0;
    virtual void run() = 0;
    virtual void collectData() = 0;
    //////////

    virtual void writeDebugData(const ResultsToWrite &resultsToWrite) const = 0;
    virtual void writeResults(const ResultsToWrite& resultsToWrite) = 0;
    virtual void writeOutletFlux() = 0;

    virtual void saveCheckpoint() const = 0;
    virtual void loadCheckpoint(const size_t iteration) = 0;

    virtual void updateRenderer(IPPRendererPtr& renderer) = 0;

    virtual void updateTransportScalarDref(const double& newPorosRef) = 0;

};

}

#endif // TRANSPORTMODULE_H
