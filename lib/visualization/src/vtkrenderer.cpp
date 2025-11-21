#include "vtkrenderer.h"


#include <memory>


#ifndef NO_RENDERING

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkImageMapper3D.h>
#include <vtkImageActor.h> // Note: this is a 3D actor (c.f. vtkImageMapper which is 2D)
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>


#include <vtkActor2D.h>
#include <vtkImageMapper.h>
#include <vtkCamera.h>
#include <vtkImageProperty.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkImageReslice.h>
#include <vtkTransform.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>



#include <vtkInteractorStyleImage.h>
#include <vtkImageResliceMapper.h>
#include <abstractscalarfield.h>
#endif

namespace IPP
{

struct IPP_VTKRendererImpl : public VTKRenderer
{

public:
    IPP_VTKRendererImpl();

    virtual void init(const size_t winX, const size_t winY,
                      const size_t nx, const size_t ny, const size_t nz,
                      const double& dx,
                      const double& min, const double& max, const size_t nContours = 0);

    virtual void updateImageData(const AbstractScalarField &imageSrc);


    virtual void render();


private:
#ifndef NO_RENDERING
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkImageData> imageData;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
#endif
};

typedef std::shared_ptr<IPP_VTKRendererImpl> LBGeoChemVTKRendererImplPtr;


IPP_VTKRendererImpl::IPP_VTKRendererImpl()
{

}



#ifdef NO_RENDERING
// replace by empty implementation
void IPP_VTKRendererImpl::init(const size_t winX, const size_t winY,
                             const size_t nx, const size_t ny, const size_t nz,
                             const double& dx,
                             const double& min, const double& max, const size_t nContours)
{

}

void IPP_VTKRendererImpl::updateImageData(const AbstractScalarField& imageSrc)
{

}

void IPP_VTKRendererImpl::render()
{

}

#else
static void createColorImage(vtkImageData* imageData)
{
    int* dims = imageData->GetDimensions();

    for (int z = 0; z < dims[2]; z++)
    {
        for (int y = 0; y < dims[1]; y++)
        {
            for (int x = 0; x < dims[0]; x++)
            {
                double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
                *pixel = (double)x / (double)dims[0];
            }
        }
    }

    imageData->Modified();
}




void IPP_VTKRendererImpl::init(const size_t winX, const size_t winY,
                             const size_t nx, const size_t ny, const size_t nz,
                             const double& dx,
                             const double& min, const double& max, const size_t nContours)
{
    imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(nx, ny, nz);
    imageData->AllocateScalars(VTK_DOUBLE, 1);
    imageData->SetSpacing(dx, dx, dx);
    //    imageData->SetOrigin(offset[0], offset[1], offset[2]);
    createColorImage(imageData);


    // Create a lookup table to map cell data to colors
    vtkSmartPointer<vtkLookupTable> lookupTable =
      vtkSmartPointer<vtkLookupTable>::New();
    lookupTable->SetTableRange(min, max);
    lookupTable->SetHueRange(0.667, 0.0);
//    lookupTable->SetBelowRangeColor(0, 0, 0, 1);
//    lookupTable->SetAboveRangeColor(1, 1, 1, 1);
    lookupTable->Build();

    vtkSmartPointer<vtkImageReslice> reslice =
            vtkSmartPointer<vtkImageReslice>::New();
    reslice->SetInputData(imageData);
    reslice->SetInterpolationModeToCubic();

    vtkSmartPointer<vtkImageMapToColors> colorMapping =
      vtkSmartPointer<vtkImageMapToColors>::New();
    colorMapping->SetLookupTable(lookupTable);
    colorMapping->SetInputConnection(reslice->GetOutputPort());




    vtkSmartPointer<vtkImageResliceMapper> mapper = vtkSmartPointer<vtkImageResliceMapper>::New();
    mapper->SetInputConnection(colorMapping->GetOutputPort());

    vtkSmartPointer<vtkImageSlice> imageSlice = vtkSmartPointer<vtkImageSlice>::New();
    imageSlice->SetMapper(mapper);
    imageSlice->GetProperty()->SetInterpolationTypeToCubic();


    vtkSmartPointer<vtkRenderer> geoRenderer = vtkSmartPointer<vtkRenderer>::New();
    geoRenderer->SetLayer(0);
    geoRenderer->SetBackground(.1,.2,.3); // Background color dark blue
    geoRenderer->AddViewProp(imageSlice);
    geoRenderer->ResetCamera();

    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetNumberOfLayers(2);
    renderWindow->AddRenderer(geoRenderer);

    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

//    vtkSmartPointer<vtkInteractorStyleImage> style =
//      vtkSmartPointer<vtkInteractorStyleImage>::New();
//    renderWindowInteractor->SetInteractorStyle(style);

//    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
//      vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview

//    renderWindowInteractor->SetInteractorStyle( style );


    if(nContours > 0)
    {
        // Create an isosurface
        vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
        contourFilter->SetInputConnection(reslice->GetOutputPort());
        contourFilter->GenerateValues(nContours, min, max); // (numContours, rangeStart, rangeEnd)
//        contourFilter->ComputeScalarsOn();


        vtkSmartPointer<vtkLookupTable> isoLUT =
          vtkSmartPointer<vtkLookupTable>::New();
        isoLUT->SetTableRange(min, max);
//        isoLUT->SetNumberOfTableValues(1);
//        isoLUT->SetTableValue(0, 1,1,1,1);
//        isoLUT->SetHueRange(-0.6667, 0.0);
        isoLUT->SetHueRange(lookupTable->GetHueRange());
        const double hueVal = 0.5;
        isoLUT->SetValueRange(hueVal, hueVal);

        isoLUT->Build();


        // Map the contours to graphical primitives
        vtkSmartPointer<vtkPolyDataMapper> contourMapper =
                vtkSmartPointer<vtkPolyDataMapper>::New();
        contourMapper->SetInputConnection(contourFilter->GetOutputPort());
        contourMapper->SetLookupTable(isoLUT);
        contourMapper->SetScalarRange(min, max);
        contourMapper->Update();

        // Create an actor for the contours
        vtkSmartPointer<vtkActor> contourActor =
                vtkSmartPointer<vtkActor>::New();
        contourActor->SetMapper(contourMapper);


        vtkSmartPointer<vtkRenderer> contourRenderer = vtkSmartPointer<vtkRenderer>::New();
        contourRenderer->SetLayer(1);
        renderWindow->AddRenderer(contourRenderer);
        contourRenderer->AddActor(contourActor);
        contourRenderer->SetActiveCamera(geoRenderer->GetActiveCamera());
//        contourRenderer->ResetCamera();

    }

    renderWindow->SetSize(winX, winY);
    renderWindow->SetFullScreen(false);

//    renderWindowInteractor->Initialize();
//    renderWindowInteractor->Start();
    renderWindow->Render();

    //    imageSlice->SetVisibility(false);
}

void IPP_VTKRendererImpl::updateImageData(const AbstractScalarField& imageSrc)
{
    int* dims = imageData->GetDimensions();

    for (int z = 0; z < dims[2]; z++)
    {
        for (int y = 0; y < dims[1]; y++)
        {
            for (int x = 0; x < dims[0]; x++)
            {
                double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
                const double value = imageSrc.get(x, y, z);
                pixel[0] = value;
            }
        }
    }

    imageData->Modified();

}

void IPP_VTKRendererImpl::render()
{
    renderWindow->Render();
//    renderWindowInteractor->Start();

}
#endif

VTKRendererPtr VTKRendererFactory::create()
{
    VTKRendererPtr renderer = std::make_shared<IPP_VTKRendererImpl>();
    return renderer;
}


}
