#ifndef IPP_H
#define IPP_H

#include "ippconfig.h"

namespace IPP
{

class IPPimpl;

class IPPRenderer;
typedef std::shared_ptr<IPPRenderer> IPPRendererPtr;

class IPP
{
public:
    IPP();
    ~IPP();

    void init(const IPPConfigPtr conf);
    void loadCheckpoint(const size_t iteration);

    void execute();

    void setRenderingModule(IPPRendererPtr renderer);
    
private:
    IPPimpl* m_pImpl;

};

} // end of namespace

#endif // IPP_H
