#ifndef IPPIMPL_H
#define IPPIMPL_H


#include "ippconfig.h"
#include "transreacstate.h"

namespace IPP
{

class IPPimpl
{
public:
    IPPimpl();
    ~IPPimpl();

    void init(const IPPConfigPtr conf);
    void loadCheckpoint(const size_t iteration);

    void execute();

    void setRenderingModule(IPPRendererPtr renderer);

    
private:
    void initialStep();
    void runStep();

    void writeResult() const;
    void saveCheckpoint() const;

    inline bool mustWriteResults(const size_t it) const;

    IPPConfigPtr m_conf;
    TransReacState* m_state;


    // hardcoded optim
    size_t m_lastPorsRefChangeIt;
    double m_oldPorosRef;
    bool updatePorosRef(const size_t it);
    void setPorosRef(const double& newPorosRef);
};

}

#endif // IPPIMPL_H
