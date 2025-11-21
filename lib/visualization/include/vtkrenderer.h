#ifndef VTKRENDERER_H
#define VTKRENDERER_H

#include <memory>
#include <ipprenderer.h>

namespace IPP
{

class VTKRenderer;
typedef std::shared_ptr<VTKRenderer> VTKRendererPtr;

namespace VTKRendererFactory
{

VTKRendererPtr create();

}



class VTKRenderer : public IPPRenderer
{
public:
    virtual void init(const size_t winX, const size_t winY,
                      const size_t nx, const size_t ny, const size_t nz,
                      const double& dx,
                      const double& min, const double& max, const size_t nContours = 0) = 0;
};

typedef std::shared_ptr<VTKRenderer> VTKRendererPtr;

}
#endif // VTKRENDERER_H

