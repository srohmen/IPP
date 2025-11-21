#include "ipp.h"

#include "ippimpl.h"

namespace IPP
{

IPP::IPP()
    : m_pImpl(new IPPimpl())
{

}

IPP::~IPP()
{
    delete m_pImpl;
}

void IPP::init(const IPPConfigPtr conf)
{
    m_pImpl->init(conf);
}

void IPP::loadCheckpoint(const size_t iteration)
{
    m_pImpl->loadCheckpoint(iteration);
}

void IPP::execute()
{
    m_pImpl->execute();
}

void IPP::setRenderingModule(IPPRendererPtr renderer)
{
    m_pImpl->setRenderingModule(renderer);
}


} // end of namespace LBGeoChem
