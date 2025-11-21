#ifndef PALABOSRESULTWRITER_HH
#define PALABOSRESULTWRITER_HH

#include "palabosresultwriter.h"

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <palabosfieldwrapper.h>
#include <arrayviewiterator.h>

#include "palaboslatticevalueaccess.h"
#include "arraydimensionconvert.h"
#include "mpi_types.h"
#include "palabostransportdata.h"
#include "simulationexchangedata.h"
#include "resultcalc.h"
#include "plbtypededuction.h"
#include "palabosconversiontools.h"
#include "createfilename.h"
#include "concarrayview.h"
#include "resultstowrite.h"
#include "ippstream.h"
#include "permeabiltyflag.h"

#include "arrayviewiteratorimpl.h"
#include "fielddecomposition.h"
#include "multitoserialdatasync.h"
#include "array_view_strm.h"

#include "mpitools.h"

namespace IPP
{


template<typename T, int nDim>
class Tensor_A_dividedBy_Scalar_B_inplace_functional2D :
        public plb::BoxProcessingFunctional2D_ST<T,T,nDim>
{
public:
    virtual void process(plb::Box2D domain, plb::ScalarField2D<T>& B,
                         plb::TensorField2D<T,nDim>& A)
    {
        plb::Dot2D offset = plb::computeRelativeDisplacement(B,A);
        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                A.get(iX+offset.x,iY+offset.y) /= B.get(iX, iY);
            }
        }
    }

    virtual Tensor_A_dividedBy_Scalar_B_inplace_functional2D<T,nDim>* clone() const
    {
        return new Tensor_A_dividedBy_Scalar_B_inplace_functional2D<T,nDim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::staticVariables;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return plb::BlockDomain::bulk;
    }
};

template<typename T, int tensorDim, size_t dim>
struct GetTensorDividedByScalarInPlaceFunction;

template<typename T, int tensorDim>
struct GetTensorDividedByScalarInPlaceFunction<T,tensorDim,2>
{
    using value = Tensor_A_dividedBy_Scalar_B_inplace_functional2D<T,tensorDim>;
};

template<typename T, int tensorDim>
struct GetTensorDividedByScalarInPlaceFunction<T,tensorDim,3>
{
    using value = plb::Tensor_A_dividedBy_Scalar_B_inplace_functional3D<T,tensorDim>;
};





template<typename T, int tensorDim, int dim>
class Tensor_A_Times_Scalar_B_InplaceFunctional :
        public PlbTypeDeduction::GetBoxProcessingFunctionalXD_ST<T,T,tensorDim,dim>::value
{
public:


    using BaseClass = typename PlbTypeDeduction::
    GetBoxProcessingFunctionalXD_ST<T, T, dim, dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;

    using Box = typename Traits::template argument<1>::type;
    using ScalarField = typename Traits::template argument<2>::type;
    using TensorField = typename Traits::template argument<3>::type;

    virtual void process(Box domain, ScalarField& B, TensorField& A)
    {
        const auto offset = plb::computeRelativeDisplacement(B,A);

        using namespace DataAccess;
        const DomainIterator<dim> end = DataAccess::end(domain);
        for(DomainIterator<dim> it = DataAccess::begin(domain);
            it < end; ++it)
        {
            const auto& posScalar = *it;
            const auto posTensor = posScalar + offset;

            plb::Array<T,dim>& vec = DataAccess::get(A, posTensor);
            const T& factor = DataAccess::get(B, posScalar);
            vec *= factor;
        }
    }

    virtual Tensor_A_Times_Scalar_B_InplaceFunctional<T,tensorDim,dim>* clone() const
    {
        return new Tensor_A_Times_Scalar_B_InplaceFunctional<T,tensorDim,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::staticVariables;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return plb::BlockDomain::bulk;
    }
};



namespace ResultWriter
{


template<typename TransportTraits, bool isTRT>
struct GetFluxCorr;


template<typename TransportTraits>
struct GetFluxCorr<TransportTraits, true>
{
    using Scalar = typename TransportTraits::Scalar;
    using DiffusionLattice = typename TransportTraits::DiffusionLattice;
    using ScalarFieldPtr = typename TransportTraits::ScalarFieldPtr;
    using DiffusionDescriptor = typename TransportTraits::DiffusionDescriptor;

    static ScalarFieldPtr get(DiffusionLattice& diffLattice)
    {
        using EF = typename DiffusionDescriptor::ExternalField;

        // for TRT the omega(Minus) are at external field cell variable
        static_assert(PlbTypeDeduction::HasTransportScalar<DiffusionDescriptor>::value, "need transportScalar");
        ScalarFieldPtr omegaMinus = plb::computeExternalScalar(diffLattice, EF::transportScalarBeginsAt);

        static_assert(PlbTypeDeduction::HasPorosity<DiffusionDescriptor>::value, "need porosity");
        ScalarFieldPtr porosity = plb::computeExternalScalar(diffLattice, EF::porosityBeginsAt);

        auto fluxFactorCalc = [](const Scalar& omega){ return (1.0 - 0.5 * omega); };
        plb::apply(fluxFactorCalc, *omegaMinus);

        ScalarFieldPtr fluxCorr = plb::multiply(*porosity, *omegaMinus);

        return fluxCorr;
    }
};

template<typename TransportTraits>
struct GetFluxCorr<TransportTraits, false>
{
    using Scalar = typename TransportTraits::Scalar;
    using DiffusionLattice = typename TransportTraits::DiffusionLattice;
    using ScalarFieldPtr = typename TransportTraits::ScalarFieldPtr;

    static ScalarFieldPtr get(DiffusionLattice& diffLattice)
    {
        // for SRT the omega is saved at the dynamic
        ScalarFieldPtr fluxCorr = plb::computeOmega(diffLattice);

        auto fluxFactorCalc = [](const Scalar& omega){ return (1.0 - 0.5 * omega); };
        plb::apply(fluxFactorCalc, *fluxCorr);

        return fluxCorr;
    }
};


template<typename TransportTraits>
struct PalabosResultWriterDetail
{
    // TODO: remove types which are not needed in detail
    using Scalar = typename TransportTraits::Scalar;

    using ScalarField = typename TransportTraits::ScalarField;
    using ScalarFieldPtr = typename TransportTraits::ScalarFieldPtr;
    using TensorField = typename TransportTraits::TensorField;
    using TensorFieldPtr = typename TransportTraits::TensorFieldPtr;
    using DiffusionDescriptor = typename TransportTraits::DiffusionDescriptor;
    using DiffusionLattice = typename TransportTraits::DiffusionLattice;
    using DiffusionLatticePtr = typename TransportTraits::DiffusionLatticePtr;
    using Vector = typename TransportTraits::Vector;
    using Box = typename TransportTraits::Box;

    static constexpr size_t dim = TransportTraits::dim;
    using PlbFactory = PalabosObjectFactory<dim>;


    using PlbTransData = PalabosTransportData<TransportTraits>;
    using FlowDiffResults = IPPResults::FlowDiffResults<TransportTraits>;

    using VtkOut = typename TransportTraits::VtkOut;


    static void writeVTKSingleFile(const Scalar& dx, const Vector& offset,
                                   const IPPResults::DensityVelocityResults<TransportTraits>& results,
                                   plb::pluint iter, const std::string& name)
    {

        const std::string fileName = CreateFileName::create(name, iter);

        const double ddx = dx;
        VtkOut vtkOut(fileName, ddx, offset);
        vtkOut.template writeData<float>(*results.density, "density");

        if(results.velocity.get() != nullptr)
        {
            // TODO: fix velocity conversion factor for LBM -> physical units
            const double velConvFac = 1.0;
            vtkOut.template writeData<TransportTraits::dim, float>(*results.velocity, "velocity", velConvFac);
        }
    }

    static void writeVTK(const Scalar& dx, const Vector& offset,
                         typename TransportTraits::ScalarField& scalarField,
                         plb::pluint iter, const std::string& name)
    {
        const std::string fileName = CreateFileName::create(name, iter);

        VtkOut vtkOut(fileName, dx, offset);
        vtkOut.template writeData<float>(scalarField, "scalar");
    }

    static inline double getViewData(const av::array_view<const double, 3>& view,
                                     const size_t x, const size_t y, const size_t z)
    {
        return view[x][y][z];
    }

    static inline double getViewData(const av::array_view<const double, 2>& view,
                                     const size_t x, const size_t y, const size_t /*z*/)
    {
        return view[x][y];
    }

    template<typename ArrView>
    static void writeVTKZonal(const Scalar& dx, const std::array<size_t, 3>& dims,
                              const plb::Array<Scalar,3>& offset,
                              const plb::pluint iter,
                              const ArrView& view,
                              const std::string& filePrefix)
    {
        vtkSmartPointer<vtkImageData> imageData = createVTKCellImage(dims, offset, dx);
        const size_t totalSize = imageData->GetNumberOfCells();

        vtkSmartPointer<vtkDoubleArray> arr =
                vtkSmartPointer<vtkDoubleArray>::New();

        arr->SetName("scalar");
        arr->SetNumberOfComponents(1);
        arr->SetNumberOfValues(totalSize);
        arr->SetNumberOfTuples(totalSize);

        double* dstPtr = arr->GetPointer(0);



        for (size_t z = 0; z < dims[2]; z++)
        {
            for (size_t y = 0; y < dims[1]; y++)
            {
                for (size_t x = 0; x < dims[0]; x++)
                {
                    *dstPtr = getViewData(view, x, y, z);
                    ++dstPtr;
                }
            }
        }

        imageData->GetCellData()->AddArray(arr);

        vtkSmartPointer<vtkXMLImageDataWriter> writer =
                vtkSmartPointer<vtkXMLImageDataWriter>::New();

        std::string filenName = plb::global::directories().getOutputDir()
                                + CreateFileName::create(filePrefix, iter) + ".vti";
        writer->SetFileName(filenName.c_str());
        writer->SetInputData(imageData);
        writer->Write();

    }

    static void writeVTK(const double& dx, const std::array<size_t, 3>& dims,
                         const plb::Array<double,3>& offset,
                         const plb::pluint iter,
                         const av::array_view<const double, 3>& view,
                         const std::string& filePrefix)
    {
        vtkSmartPointer<vtkImageData> imageData =
                vtkSmartPointer<vtkImageData>::New();
        imageData->SetDimensions(dims[0], dims[1], dims[2]);
        imageData->AllocateScalars(VTK_DOUBLE, 1);
        imageData->SetSpacing(dx, dx, dx);
        imageData->SetOrigin(offset[0], offset[1], offset[2]);


        for (size_t z = 0; z < dims[2]; z++)
        {
            for (size_t y = 0; y < dims[1]; y++)
            {
                for (size_t x = 0; x < dims[0]; x++)
                {
                    double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
                    pixel[0] = view[x][y][z];
                }
            }
        }

        vtkSmartPointer<vtkXMLImageDataWriter> writer =
                vtkSmartPointer<vtkXMLImageDataWriter>::New();

        std::string filenName = plb::global::directories().getOutputDir()
                                + CreateFileName::create(filePrefix, iter) + ".vti";
        writer->SetFileName(filenName.c_str());
        writer->SetInputData(imageData);
        writer->Write();

    }


    static void writeVTKZonal(const Scalar& dx, const std::array<size_t, 3>& dims,
                              const plb::Array<Scalar,3>& offset,
                              const plb::pluint iter,
                              const std::vector<double>& rawData,
                              const std::string& filePrefix)
    {
        const av::bounds<3> bounds = { (ptrdiff_t)dims[0], (ptrdiff_t)dims[1], (ptrdiff_t)dims[2] };
        const av::array_view<const double, 3> view(rawData, bounds);
        writeVTKZonal(dx, dims, offset, iter, view, filePrefix);
    }

    static void writeVTK(const Scalar& dx, const ArrayDimensionConvert& indexConv,
                         const plb::Array<Scalar,3>& offset,
                         const plb::pluint iter,
                         const std::vector<double>::const_iterator& begin,
                         const std::string& filePrefix)
    {
        size_t nx, ny, nz;
        indexConv.getSize(nx, ny, nz);

        vtkSmartPointer<vtkImageData> imageData =
                vtkSmartPointer<vtkImageData>::New();
        imageData->SetDimensions(nx, ny, nz);
        imageData->AllocateScalars(VTK_DOUBLE, 1);
        imageData->SetSpacing(dx, dx, dx);
        imageData->SetOrigin(offset[0], offset[1], offset[2]);

        int* dims = imageData->GetDimensions();

        for (int z = 0; z < dims[2]; z++)
        {
            for (int y = 0; y < dims[1]; y++)
            {
                for (int x = 0; x < dims[0]; x++)
                {
                    double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
                    const size_t index = indexConv.calcIndex(x, y, z);
                    const std::vector<double>::const_iterator currIt = begin + index;
                    pixel[0] = *currIt;
                }
            }
        }

        vtkSmartPointer<vtkXMLImageDataWriter> writer =
                vtkSmartPointer<vtkXMLImageDataWriter>::New();

        std::string filenName = plb::global::directories().getOutputDir()
                                + CreateFileName::create(filePrefix, iter) + ".vti";
        writer->SetFileName(filenName.c_str());
        writer->SetInputData(imageData);
        writer->Write();

    }

    static void writeVTKFromRaw(const Scalar& dx, const ArrayDimensionConvert& indexConv,
                                const plb::Array<Scalar,3>& offset,
                                const plb::pluint iter,
                                const std::vector<double>& rawData,
                                const std::vector<std::string>& dataSetNames,
                                const std::string& filePrefix)
    {
        size_t nx, ny, nz;
        indexConv.getSize(nx, ny, nz);
        const std::array<size_t, 3> dims = {{nx, ny, nz}};
        vtkSmartPointer<vtkImageData> imageData = createVTKCellImage(dims, offset, dx);
        const size_t totalSize = imageData->GetNumberOfCells();

        for(size_t i = 0; i < dataSetNames.size(); ++i)
        {
            const size_t index = totalSize * i;

            vtkSmartPointer<vtkDoubleArray> arr =
                    vtkSmartPointer<vtkDoubleArray>::New();

            const std::string& name = dataSetNames[i];
            arr->SetName(name.c_str());
            arr->SetNumberOfComponents(1);
            arr->SetNumberOfValues(totalSize);
            arr->SetNumberOfTuples(totalSize);

            double* dstPtr = arr->GetPointer(0);
            const std::vector<double>::const_iterator begin = rawData.begin() + index;
            const std::vector<double>::const_iterator end = begin + totalSize;

            std::copy(begin, end, dstPtr);


            imageData->GetCellData()->AddArray(arr);
        }

        vtkSmartPointer<vtkXMLImageDataWriter> writer =
                vtkSmartPointer<vtkXMLImageDataWriter>::New();

        std::string filenName = plb::global::directories().getOutputDir()
                                + CreateFileName::create(filePrefix, iter) + ".vti";
        writer->SetFileName(filenName.c_str());
        writer->SetInputData(imageData);

        writer->Write();

    }


    static vtkSmartPointer<vtkImageData> createVTKCellImage(const std::array<size_t, 3>& dims,
                                                            const plb::Array<Scalar, 3>& offset,
                                                            const Scalar& dx)
    {
        vtkSmartPointer<vtkImageData> imageData =
                vtkSmartPointer<vtkImageData>::New();
        imageData->SetSpacing(dx, dx, dx);

        const double dx_2 = dx * 0.5;
        if(dims[2] == 1)
        {
            imageData->SetDimensions(dims[0]+1, dims[1]+1, 1);
            imageData->SetOrigin(offset[0] - dx_2, offset[1] - dx_2, offset[2]);
        }
        else
        {
            imageData->SetDimensions(dims[0]+1, dims[1]+1, dims[2] + 1);
            imageData->SetOrigin(offset[0] - dx_2, offset[1] - dx_2, offset[2] - dx_2);
        }

        return imageData;
    }

    struct Nodal
    {
        static vtkSmartPointer<vtkImageData> createImage(const std::array<size_t, 3>& dims,
                                                         const plb::Array<double, 3>& offset,
                                                         const double& dx)
        {
            vtkSmartPointer<vtkImageData> image =
                    vtkSmartPointer<vtkImageData>::New();
            image->SetDimensions(dims[0], dims[1], dims[2]);
            // image->AllocateScalars(VTK_DOUBLE, 1);
            image->SetSpacing(dx, dx, dx);
            image->SetOrigin(offset[0], offset[1], offset[2]);
            return image;
        }

        static size_t getSize(vtkImageData& image)
        {
            return image.GetNumberOfPoints();
        }

        template<typename Arr>
        static void addArray(vtkImageData& image, Arr& arr)
        {
            image.GetPointData()->AddArray(arr);
        }
    };

    struct Zonal
    {
        static vtkSmartPointer<vtkImageData> createImage(const std::array<size_t, 3>& dims,
                                                         const plb::Array<Scalar, 3>& offset,
                                                         const Scalar& dx)
        {
            vtkSmartPointer<vtkImageData> image = createVTKCellImage(dims, offset, dx);
            return image;
        }

        static size_t getSize(vtkImageData& image)
        {
            return image.GetNumberOfCells();
        }

        template<typename Arr>
        static void addArray(vtkImageData& image, Arr& arr)
        {
            image.GetCellData()->AddArray(arr);
        }
    };


    template<typename NodalZonalImageHandler>
    static void writeVTKFromRaw(const std::array<size_t, 3>& dims,
                                const plb::Array<Scalar,3>& offset,
                                const av::array_view<const double, 4> view,
                                const std::vector<std::string>& dataSetNames,
                                const size_t nComps, const plb::pluint iter,
                                const std::string& filePrefix, const double& dx)
    {
        vtkSmartPointer<vtkImageData> imageData = NodalZonalImageHandler::createImage(dims, offset, dx);
        const size_t totalSize = NodalZonalImageHandler::getSize(*imageData);

        for(size_t i = 0; i < nComps; ++i)
        {
            vtkSmartPointer<vtkDoubleArray> arr =
                    vtkSmartPointer<vtkDoubleArray>::New();

            const std::string& name = dataSetNames[i];
            arr->SetName(name.c_str());
            arr->SetNumberOfComponents(1);
            arr->SetNumberOfValues(totalSize);
            arr->SetNumberOfTuples(totalSize);

            double* dstPtr = arr->GetPointer(0);

            const av::array_view<const double, 3> slice = view[i];

            for (size_t z = 0; z < dims[2]; z++)
            {
                for (size_t y = 0; y < dims[1]; y++)
                {
                    for (size_t x = 0; x < dims[0]; x++)
                    {
                        *dstPtr = slice[x][y][z];
                        ++dstPtr;
                    }
                }
            }

            NodalZonalImageHandler::addArray(*imageData, arr);

        }

        vtkSmartPointer<vtkXMLImageDataWriter> writer =
                vtkSmartPointer<vtkXMLImageDataWriter>::New();

        std::string filenName = plb::global::directories().getOutputDir()
                                + CreateFileName::create(filePrefix, iter) + ".vti";
        writer->SetFileName(filenName.c_str());
        writer->SetInputData(imageData);

        writer->Write();
    }

    static void writeVTKFromRaw(const double& dx, const std::array<size_t, 3>& dims,
                                const plb::Array<Scalar,3>& offset,
                                const plb::pluint iter,
                                const std::vector<double>& rawData,
                                const std::string& filePrefix)
    {
        const av::bounds<3> bounds = { (ptrdiff_t)dims[0], (ptrdiff_t)dims[1], (ptrdiff_t)dims[2] };
        const av::array_view<const double, 3> view(rawData, bounds);
        writeVTK(dx, dims, offset, iter, view, filePrefix);
    }

    template<typename NodalZonalImageHandler>
    static void writeVTKFromRaw(const double& dx, const std::array<size_t, 3>& dims,
                                const plb::Array<Scalar,3>& offset,
                                const plb::pluint iter,
                                const std::vector<double>& rawData,
                                const std::vector<std::string>& dataSetNames,
                                const std::string& filePrefix)
    {
        const size_t nComps = dataSetNames.size();
        const auto view = ConcArrayViewTool::makeConcView<3>(nComps, dims[0], dims[1], dims[2], rawData);

        writeVTKFromRaw<NodalZonalImageHandler>(dims, offset, view, dataSetNames,
                                                nComps, iter, filePrefix, dx);

    }


    static void writeTransportScalar(const Scalar& dx, const plb::Array<Scalar,3>& offset, const size_t it,
                                     PalabosTransportData<TransportTraits>& palabosData)
    {
        if(palabosData.diffLattices.empty() == false)
        {
            typename TransportTraits::DiffusionLatticePtr& lattice = palabosData.diffLattices.front();
            typename TransportTraits::ScalarFieldPtr scalarField;
            PalabosLatticeValueAccess::computeTransportScalarIfPossible<typename TransportTraits::DiffusionDescriptor,
                    typename TransportTraits::DiffusionLattice>(*lattice, scalarField);

            if(scalarField.get() != nullptr)
            {
                writeVTK(dx, offset, *scalarField, it, "trans_scalar_");
            }

        }
    }


    static Box createBoundaryDomain(const size_t nx, const size_t ny,
                                    const size_t nz, const size_t direction)
    {
        const size_t offsetFromEnd = 1;
        switch(direction)
        {
            case 0:
            {
                return PalabosObjectFactoryUtils<Scalar, dim>::createBox(nx - offsetFromEnd, nx - offsetFromEnd,
                                                                         0, ny-1,
                                                                         0, nz-1);
            }

            case 1:
            {
                return PalabosObjectFactoryUtils<Scalar, dim>::createBox(0, nx-1,
                                                                         ny - offsetFromEnd, ny - offsetFromEnd,
                                                                         0, nz-1);
            }

            case 2:
            {
                return PalabosObjectFactoryUtils<Scalar, dim>::createBox(0, nx-1,
                                                                         0, ny-1,
                                                                         nz - offsetFromEnd, nz - offsetFromEnd);
            }

            default:
            {
                throw std::runtime_error("unknown dimension for diffusion coeff calculation: "
                                         + std::to_string(dim));
            }

        }
    }


    static Box createBoundaryDomainCornerExcluded(const size_t nx, const size_t ny,
                                                  const size_t nz, const size_t direction)
    {
        switch(direction)
        {
            case 0:
            {
                return PalabosObjectFactoryUtils<Scalar, dim>::createBox(nx-1, nx-1, 1, ny-2, 1, nz-2);
            }

            case 1:
            {
                return PalabosObjectFactoryUtils<Scalar, dim>::createBox(1, nx-2, ny-1, ny-1, 1, nz-2);
            }

            case 2:
            {
                return PalabosObjectFactoryUtils<Scalar, dim>::createBox(1, nx-2, 1, ny-2, nz-1, nz-1);
            }

            default:
            {
                throw std::runtime_error("unknown dimension for diffusion coeff calculation: "
                                         + std::to_string(dim));
            }

        }
    }

};

template<typename T, size_t dim>
class DecomposePermInfo : public PlbTypeDeduction::
        GetScalarFieldBoxProcessingFunctionalXD<T,dim>::value
{
public:
    DecomposePermInfo()
    {

    }

    virtual DecomposePermInfo<T,dim>* clone() const override
    {
        return new DecomposePermInfo<T,dim>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        assert(modified.size() == 4);
        for(plb::modif::ModifT& val : modified)
        {
            val = plb::modif::nothing;
        }
    }

private:
    using BaseClass =
    typename PlbTypeDeduction::GetScalarFieldBoxProcessingFunctionalXD<T,dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::process)>;
    using Box = typename Traits::template argument<1>::type;
    using InputFieldVec = typename Traits::template argument<2>::type;

    using ScalarField = typename PlbTypeDeduction::GetScalarField_XD<T, dim>::value;

public:
    virtual void process(Box domain, InputFieldVec fields) override
    {
        assert(fields.size() == 4);

        ScalarField& permMask = *fields[0];
        ScalarField& isPermeable = *fields[1];
        ScalarField& isInterface = *fields[2];
        ScalarField& permNeighCount = *fields[3];

        const auto offsetIsPerm = plb::computeRelativeDisplacement(permMask, isPermeable);
        const auto offsetIsInterface = plb::computeRelativeDisplacement(permMask, isInterface);
        const auto offsetPermNeighCount = plb::computeRelativeDisplacement(permMask, permNeighCount);

        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const auto& pos = *it;

            const T& maskVal = DataAccess::get(permMask, pos);

            T& isPermeableVal = DataAccess::get(isPermeable, pos + offsetIsPerm);
            T& isInterfaceVal = DataAccess::get(isInterface, pos + offsetIsInterface);
            T& permNeighCountVal = DataAccess::get(permNeighCount, pos + offsetPermNeighCount);

            isPermeableVal = (maskVal & PF_isPermeable) ? true : false;
            isInterfaceVal = (maskVal & PF_isInterface) ? true : false;
            permNeighCountVal = (unsigned int)maskVal >> 8;

        }
    }

private:

};


template<typename TransportTraits>
void PalabosResultWriter<TransportTraits>::writeDebugData(const ResultsToWrite& resultsToWrite,
                                                          const Scalar& dx,
                                                          const plb::Array<Scalar, 3>& offset3d,
                                                          const std::array<size_t, 3>& dims,
                                                          const SimulationExchangeData& simData,
                                                          const PlbTransData& palabosData,
                                                          const PlbFactory& factory)
{
    MPI_CHECK_SYNC;

    const size_t it = simData.getIteration();

    Vector offset;
    PalabosConversionTools::convertVector(offset3d, offset);


#ifdef IPP_DEBUG_ENABLED

    if(resultsToWrite.contains(RT_EnabledCells))
    {
        const DebugData::EnabledFieldVec& enabledFields = simData.debugData.enabledCellsSteps;

        for(size_t iEnab = 0; iEnab < enabledFields.size(); ++iEnab)
        {
            const std::vector<char>& localFieldRaw = enabledFields[iEnab];

            // TODO: optimize in order to prevent char -> double conversion
            std::vector<double> localField(localFieldRaw.size());
            for(size_t i = 0; i < localFieldRaw.size(); ++i)
            {
                localField[i] = localFieldRaw[i];
            }

            const FieldDecomposition& decomp = *simData.getDecomp();

            assert(decomp.getOwnDomain().getDiagonalSize() == (long)localField.size());

            std::vector<double> globalField;
            MultiToSerialDataSync<double> sync(globalField, decomp);
            sync.pullToRoot(localField);

            if(plb::global::mpi().isMainProcessor())
            {
                using PlbWriterDetail = PalabosResultWriterDetail<TransportTraits>;
                const std::string filePrefix = "enabled_" + std::to_string(iEnab) + "_";
                PlbWriterDetail::writeVTKZonal(dx, dims, offset3d, it, globalField, filePrefix);
            }
        }
    }

#endif


    if(resultsToWrite.contains(RT_PermeabilityInfo))
    {
        using MaskType = typename PlbTransData::MaskType;
        using MaskField = typename PlbTypeDeduction::GetMultiScalarField_XD<MaskType, dim>::value;
        using VtkOutMask = typename PlbTypeDeduction::GetVtkOutXD<MaskType, dim>::value;

        const std::string fileName = CreateFileName::create("perm_infos_", it);
        const double cdx = dx;
        VtkOutMask vtkOut(fileName, cdx, offset);

        MaskField permeable = factory.template createMultiScalarField<MaskType>(0);
        MaskField interface = factory.template createMultiScalarField<MaskType>(0);
        MaskField neigh = factory.template createMultiScalarField<MaskType>(0);

        std::vector<MaskField*> fields;
        fields.push_back(&(*palabosData.permMask));
        fields.push_back(&permeable);
        fields.push_back(&interface);
        fields.push_back(&neigh);

        plb::applyProcessingFunctional(new DecomposePermInfo<MaskType, dim>,
                                       permeable.getBoundingBox(), fields);

        vtkOut.template writeData<MaskType>(permeable, "permeable");
        vtkOut.template writeData<MaskType>(interface, "interface");
        vtkOut.template writeData<MaskType>(neigh, "nNeighbors");
    }
}

template<typename TransportTraits>
void PalabosResultWriter<TransportTraits>::writeResults(const ResultsToWrite& resultsToWrite,
                                                        const Scalar& dx,
                                                        const plb::Array<Scalar, 3>& offset3d,
                                                        const std::array<size_t, 3>& dims,
                                                        const FlowDiffResults& results,
                                                        const SimulationExchangeData& simData,
                                                        const PlbTransData& palabosData,
                                                        const PlbFactory& factory,
                                                        std::shared_ptr<AuxResultProcessing>& resultProcess,
                                                        ConvCheck &convergenceCheck)
{
    using PlbWriterDetail = PalabosResultWriterDetail<TransportTraits>;

    const size_t it = simData.getIteration();
    const double time = simData.getTime();
    Vector offset;
    PalabosConversionTools::convertVector(offset3d, offset);

    if(resultsToWrite.contains(RT_FluidDensityAndFlux))
    {
        const IPPResults::DensityVelocityResults<TransportTraits>& flowResults = results.flowResults;
        PlbWriterDetail::writeVTKSingleFile(dx, offset, flowResults, it, "fluidDens_");
    }

    const double ddx = dx;


    if(resultsToWrite.contains(RT_LBM_Velocity))
    {

        const typename PlbTransData::DiffusionLatticePtrVec& diffLattices = palabosData.diffLattices;

        const std::string fileName = CreateFileName::create("velocities_", it);
        VtkOut vtkOut(fileName, ddx, offset);

        for(size_t iSpecies = 0; iSpecies < diffLattices.size(); ++iSpecies)
        {
            const DiffusionLatticePtr diffLattice = diffLattices[iSpecies];

            ScalarField rhoBar = factory.template createMultiScalarField<Scalar>();
            TensorField velocity = factory.template createMultiTensorField<Scalar, dim>();
            plb::computeRhoBarJ(*diffLattice, rhoBar, velocity, diffLattice->getBoundingBox());

            const std::vector<std::string>& speciesNames = palabosData.componentNames;
            const std::string& arrName = speciesNames[iSpecies];

            vtkOut.template writeData<TransportTraits::dim, float>(velocity, arrName);
        }
    }

    if(resultsToWrite.contains(RT_LBM_Gradient))
    {
        const std::string fileName = CreateFileName::create("gradients_", it);
        VtkOut vtkOut(fileName, ddx, offset);

        const std::vector<std::string>& speciesNames = palabosData.componentNames;
        for(size_t iSpecies = 0; iSpecies < speciesNames.size(); ++iSpecies)
        {
            const ScalarField& conc = palabosData.postReactionConcentrations[iSpecies];
            const TensorFieldPtr grad = plb::computeGradient(const_cast<ScalarField&>(conc));
            const std::string& arrName = speciesNames[iSpecies];
            vtkOut.template writeData<TransportTraits::dim, float>(*grad, arrName);
        }
    }


    if(resultsToWrite.contains(RT_LBM_PostTransportConc))
    {

        const std::string fileName = CreateFileName::create("ptc_", it);
        VtkOut vtkOut(fileName, ddx, offset);

        const std::vector<std::string>& compNames = palabosData.componentNames;
        for(size_t iComp = 0; iComp < compNames.size(); ++iComp)
        {
            const ScalarField& conc = palabosData.postTransportConcentrations[iComp];
            const std::string &arrName = compNames[iComp];
            vtkOut.template writeData<float>(const_cast<ScalarField&>(conc), arrName);
        }
    }

    if(resultsToWrite.contains(RT_LBM_TransPoros))
    {
        const std::string fileName = CreateFileName::create("poros_trans_", it);
        VtkOut vtkOut(fileName, ddx, offset);

        const std::vector<std::string>& speciesNames = palabosData.componentNames;
        for(size_t iSpecies = 0; iSpecies < speciesNames.size(); ++iSpecies)
        {
            const DiffusionLatticePtr& diffLattice = palabosData.diffLattices[iSpecies];
            ScalarFieldPtr poros;
            PalabosLatticeValueAccess::computePorosityIfPossible<DiffusionDescriptor>(*diffLattice, poros);

            if(poros.get())
            {
                const std::string& arrName = speciesNames[iSpecies];
                vtkOut.template writeData<float>(*poros, arrName);
            }
        }

    }



    if(resultsToWrite.contains(RT_LBM_Source))
    {
        const std::string fileName = CreateFileName::create("source_", it);
        VtkOut vtkOut(fileName, ddx, offset);

        const std::vector<std::string>& compNames = palabosData.componentNames;
        for(size_t iComp = 0; iComp < compNames.size(); ++iComp)
        {
            const DiffusionLatticePtr& diffLattice = palabosData.diffLattices[iComp];
            const ScalarFieldPtr source =
                    plb::computeExternalScalar(*diffLattice,
                                               DiffusionDescriptor::ExternalField::scalarBeginsAt);

            const std::string& arrName = compNames[iComp];
            vtkOut.template writeData<float>(*source, arrName);
        }

    }


    using ArrayView = av::array_view<const double, dim+1>;

    const FieldDecomposition &decomp = *simData.getDecomp();
    const IPPBox3DLong& ownDomain = decomp.getOwnDomain();
    const IPPVector3DLong localSize = ownDomain.getSize();
    const IPPVector3DLong& location = ownDomain.lower;

    const std::vector<std::string>& compNames = simData.getCompNames();


    if(resultsToWrite.contains(RT_Budget))
    {
        std::vector<ScalarField> postReacTotals;
        std::vector<ScalarField> postTransTotals;

        {
            const ScalarField& porosity = palabosData.getPorosity();

            for(size_t iComp = 0; iComp < compNames.size(); ++iComp)
            {
                const ScalarField& postReacConc = palabosData.postReactionConcentrations[iComp];
                postReacTotals.push_back(postReacConc);
                ScalarField& postReacTotal = postReacTotals.back();
                plb::multiplyInPlace(postReacTotal, const_cast<ScalarField&>(porosity));

                const ScalarField& postTransConc = palabosData.postTransportConcentrations[iComp];
                postTransTotals.push_back(postTransConc);
                ScalarField& postTransTotal = postTransTotals.back();
                plb::multiplyInPlace(postTransTotal, const_cast<ScalarField&>(porosity));

            }
        }



        {

            const PhaseNameToInfos& phaseInfos = simData.getPhaseNameToInfos();
            const std::vector<double>& precipPhases = simData.getPrecipField();

            const av::bounds<dim+1> bounds =
                    av::MakeBounds<dim>::get(phaseInfos.size(), localSize[0], localSize[1], localSize[2]);
            const ArrayView phasesView(precipPhases, bounds);

            size_t iPhase = 0;

            for(PhaseNameToInfos::const_iterator phaseIt = phaseInfos.begin();
                phaseIt != phaseInfos.end(); ++phaseIt, ++iPhase)
            {
                assert(simData.getPhaseNames().at(iPhase) == phaseIt->first);

                ScalarField molePhase = factory.template createMultiScalarField<Scalar>();
                const auto phaseSlice = phasesView[iPhase];
                plb::setToFunction(molePhase, molePhase.getBoundingBox(),
                                   SetScalarsFromArray<const double, dim>(phaseSlice, location));

                const PhaseInfo& info = phaseIt->second;
                const PhaseInfo::ElemIdToStoich& stoich = info.stoich;

                assert(compNames.size() == stoich.size());

                for(size_t iComp = 0; iComp < stoich.size(); ++iComp)
                {
                    const Scalar& fac = stoich[iComp];
                    const ScalarFieldPtr moleElemInPhase = plb::multiply(molePhase, fac);

                    ScalarField& postReacTotal = postReacTotals[iComp];
                    plb::addInPlace(postReacTotal, *moleElemInPhase);

                    ScalarField& postTransTotal = postTransTotals[iComp];
                    plb::addInPlace(postTransTotal, *moleElemInPhase);
                }

            }
        }



        plb::plb_ofstream file;
        std::ios_base::openmode flag = std::ios_base::trunc;
        if(it > 0)
        {
            flag = std::ios_base::app;
        }

        const boost::filesystem::path resultDir(plb::global::directories().getOutputDir());
        const boost::filesystem::path resultOut = resultDir / "budget.csv";
        file.open(resultOut.string().c_str(), flag);

        // write header
        if(it == 0)
        {
            file << "iteration,time,";

            for(size_t iComp = 0; iComp < compNames.size(); ++iComp)
            {
                const std::string& name = compNames[iComp];
                file << name << "_t" << "," << name << "_r" << "," ;
            }

            file << std::endl;
        }

        file << it << "," << time << ",";


        {
            std::ios_base::fmtflags f( std::cout.flags() );

            for(size_t iComp = 0; iComp < postReacTotals.size(); ++iComp)
            {
                ScalarField& postTransTotal = postTransTotals[iComp];
                const Scalar postTransSum = plb::computeSum(postTransTotal);
                file << std::setprecision(std::numeric_limits<Scalar>::max_digits10) << postTransSum << ",";

                ScalarField& postReacTotal = postReacTotals[iComp];
                const Scalar postReacSum = plb::computeSum(postReacTotal);
                file << std::setprecision(std::numeric_limits<Scalar>::max_digits10) << postReacSum
                     << ",";
            }

            std::cout.flags( f );
        }

        file << std::endl;
    }


    if(resultsToWrite.contains(RT_Porosity))
    {
        const std::string fileName = CreateFileName::create("poros_", it);
        typename TransportTraits::VtkOut vtkOut(fileName, ddx, offset);
        const ScalarField& porosity = palabosData.getPorosity();
        vtkOut.template writeData<float>(const_cast<ScalarField&>(porosity), "porosity");


        if(plb::global::mpi().isMainProcessor())
        {
            if(resultProcess)
            {
                assert(false);
            }

        }

    }

    if(resultsToWrite.contains(RT_CapillaryPorosity))
    {
        const std::vector<double>& capPoros = simData.getCapillaryPorosity();
        assert(decomp.getOwnDomain().getDiagonalSize() == (long)capPoros.size());

        std::vector<double> capPorosGlobal;
        MultiToSerialDataSync<double> sync(capPorosGlobal, decomp);
        sync.pullToRoot(capPoros);

        if(plb::global::mpi().isMainProcessor())
        {
            using PlbWriterDetail = PalabosResultWriterDetail<TransportTraits>;
            const std::string filePrefix = "capillary_poros_";
            PlbWriterDetail::writeVTKFromRaw(ddx, dims, offset3d, it, capPorosGlobal, filePrefix);
        }

    }

    if(resultsToWrite.contains(RT_InertFrac))
    {
        const std::vector<double>& inertVolFrac = simData.getInertVolFrac();
        assert(decomp.getOwnDomain().getDiagonalSize() == (long)inertVolFrac.size());

        std::vector<double> valuesGlobal;
        MultiToSerialDataSync<double> sync(valuesGlobal, decomp);
        sync.pullToRoot(inertVolFrac);

        if(plb::global::mpi().isMainProcessor())
        {
            using PlbWriterDetail = PalabosResultWriterDetail<TransportTraits>;
            const std::string filePrefix = "inert_frac_";
            PlbWriterDetail::writeVTKFromRaw(ddx, dims, offset3d, it, valuesGlobal, filePrefix);
        }
    }


    if(resultsToWrite.contains(RT_DistanceField))
    {
        const std::vector<double>& distFieldLocal = simData.getDistField();
        assert(decomp.getOwnDomain().getDiagonalSize() == (long)distFieldLocal.size());

        std::vector<double> distFieldGlobal;
        MultiToSerialDataSync<double> sync(distFieldGlobal, decomp);
        sync.pullToRoot(distFieldLocal);

        if(plb::global::mpi().isMainProcessor())
        {
            using PlbWriterDetail = PalabosResultWriterDetail<TransportTraits>;
            const std::string filePrefix = "dist_field_";
            PlbWriterDetail::writeVTKFromRaw(ddx, dims, offset3d, it, distFieldGlobal, filePrefix);
        }
    }


    if(resultsToWrite.contains(RT_DiffusionCoef))
    {
        const std::vector<double>& diffCoefsLocal = simData.getDiffusionCoefs();
        assert(decomp.getOwnDomain().getDiagonalSize() == (long)diffCoefsLocal.size());

        std::vector<double> diffCoefsGlobal;
        MultiToSerialDataSync<double> sync(diffCoefsGlobal, decomp);
        sync.pullToRoot(diffCoefsLocal);

        if(plb::global::mpi().isMainProcessor())
        {
            using PlbWriterDetail = PalabosResultWriterDetail<TransportTraits>;
            const std::string filePrefix = "diff_coefs_";
            PlbWriterDetail::writeVTKFromRaw(ddx, dims, offset3d, it, diffCoefsGlobal, filePrefix);
        }
    }

    if(resultsToWrite.contains(RT_Phase))
    {
        const std::vector<std::string>& phaseNames = simData.getPhaseNames();
        const size_t nPhases = phaseNames.size();

        const std::vector<double>& localField = simData.getPrecipField();
        FieldDecomposition phaseDecomp = *simData.getDecomp();
        phaseDecomp.correctByTensorDim(nPhases);

        std::vector<double> globalPrecipField;
        MultiToSerialDataSync<double> sync(globalPrecipField, phaseDecomp, nPhases);
        sync.pullToRoot(localField);


        if(plb::global::mpi().isMainProcessor())
        {
            PlbWriterDetail::template writeVTKFromRaw<typename PlbWriterDetail::Nodal>
                    (ddx, dims, offset3d, it, globalPrecipField, phaseNames, "phases_");

        }

    }

    {
        const AuxDataVec& auxDataVec = simData.getAuxData();


        if(resultProcess && plb::global::mpi().isMainProcessor())
        {
            std::vector<std::string> names;
            for(size_t iData = 0; iData < auxDataVec.size(); ++iData)
            {
                names.push_back(auxDataVec[iData].name);
            }

            resultProcess->begin(it, time, names);
        }


        for(size_t iData = 0; iData < auxDataVec.size(); ++iData)
        {
            const AuxDataName& auxData = auxDataVec[iData];
            const std::string& name = auxData.name;
            const std::vector<double>& localData = auxData.data;

            std::vector<double> globalData;
            MultiToSerialDataSync<double> sync(globalData, decomp);
            sync.pullToRoot(localData);

            if(plb::global::mpi().isMainProcessor())
            {
                PlbWriterDetail::writeVTKFromRaw(ddx, dims, offset3d, it, globalData, name + "_");

                if(resultProcess)
                {
                    resultProcess->process(name, globalData);
                }
            }
        }

        if(resultProcess && plb::global::mpi().isMainProcessor())
        {
            resultProcess->end();
        }
    }


    if(resultsToWrite.contains(RT_ComponentConcentration))
    {
        const std::string fileName = CreateFileName::create("comp_", it);
        VtkOut vtkOut(fileName, ddx, offset);

        for(size_t iComp = 0; iComp < compNames.size(); ++iComp)
        {
            const ScalarField& conc = palabosData.postReactionConcentrations[iComp];
            const std::string& arrName = compNames[iComp];
            vtkOut.template writeData<float>(const_cast<ScalarField&>(conc), arrName);
        }
    }

    if(resultsToWrite.contains(RT_DiffusiveTransportScalar))
    {
        const std::string fileName = CreateFileName::create("diff_transport_scalar_", it);
        typename TransportTraits::VtkOut vtkOut(fileName, ddx, offset);
        for(size_t iSpecies = 0; iSpecies < palabosData.diffLattices.size(); ++iSpecies)
        {
            const DiffusionLatticePtr& diffLattice = palabosData.diffLattices[iSpecies];
            const std::string& name = palabosData.componentNames[iSpecies];
            ScalarFieldPtr result = plb::computeExternalScalar(*diffLattice,
                                                               DiffusionDescriptor::ExternalField::transportScalarBeginsAt);


            vtkOut.template writeData<float>(*result, name);
        }
    }

    writeOutletFlux(dx, offset3d, dims, simData, palabosData, factory, true, convergenceCheck);

    MPI_CHECK_SYNC;
}

template<typename TransportTraits>
void PalabosResultWriter<TransportTraits>::writeOutletFlux(const Scalar& dx,
                                                           const plb::Array<Scalar, 3>& offset3d,
                                                           const std::array<size_t, 3>& dims,
                                                           const SimulationExchangeData& simData,
                                                           const PlbTransData& palabosData,
                                                           const PlbFactory& factory,
                                                           const bool writeFluxImage,
                                                           ConvCheck &convergenceCheck)
{
    if(simData.getDimToDiffSpecies().empty() == false)
    {
        plb::plb_ofstream file;
        std::ios_base::openmode flag = std::ios_base::trunc;

        const size_t it = simData.getIteration();
        if(it > 0)
        {
            flag = std::ios_base::app;
        }

        const boost::filesystem::path resultDir(plb::global::directories().getOutputDir());
        const boost::filesystem::path resultOut = resultDir / "outlet_diff_flux";
        file.open(resultOut.string().c_str(), flag);


        std::unique_ptr<VtkOut> vtkOut;
        if(writeFluxImage)
        {
            Vector offset;
            PalabosConversionTools::convertVector(offset3d, offset);

            const std::string fileName = CreateFileName::create("diff_flux_", it);
            vtkOut.reset(new VtkOut(fileName, dx, offset));
        }

        bool allConverged = true;

        for(auto diffIt = simData.getDimToDiffSpecies().begin();
            diffIt != simData.getDimToDiffSpecies().end(); ++diffIt)
        {
            const size_t diffDim = diffIt->first;
            const size_t iSpecies = diffIt->second;

            const std::string& name = palabosData.componentNames[iSpecies];
            const DiffusionLatticePtr diffLattice = palabosData.diffLattices[iSpecies];

            using PlbWriterDetail = PalabosResultWriterDetail<TransportTraits>;
            const Box avgDomain =
                    PlbWriterDetail::createBoundaryDomain(dims[0], dims[1], dims[2], diffDim);

            const Box wholeLattice = diffLattice->getBoundingBox();

            ScalarField rhoBar = factory.template createMultiScalarField<Scalar>();
            TensorField j = factory.template createMultiTensorField<Scalar, dim>();
            plb::computeRhoBarJ(*diffLattice, rhoBar, j, wholeLattice);

            ScalarFieldPtr fluxCorr = GetFluxCorr<TransportTraits, TransportTraits::isTRT>::get(*diffLattice);

            const Scalar dt = simData.getDtCurr();
            const Scalar conv = dx/dt;
            plb::multiplyInPlace(*fluxCorr, conv); // convert to physical units

            using MultiplyFunc = Tensor_A_Times_Scalar_B_InplaceFunctional<Scalar, dim, dim>;
            plb::applyProcessingFunctional(new MultiplyFunc, j.getBoundingBox(), *fluxCorr, j);

            const ScalarFieldPtr flux = plb::extractComponent(j, diffDim);
            const Scalar avgFlux = plb::computeAverage(*flux, avgDomain);

            const size_t iField = diffIt - simData.getDimToDiffSpecies().begin();
            const double time = simData.getTime();
            const bool isConverged = convergenceCheck.isConverged(iField, avgFlux, time);
            if(isConverged == false)
            {
                allConverged = false;
            }
            convergenceCheck.setOldData(iField, avgFlux, time);

            file << it << "\t";
            file << std::setprecision (15) << time
                 << "\t" << name << "\t";
            file << std::setprecision (15) << avgFlux << std::endl;

            if(writeFluxImage)
            {
                vtkOut->template writeData<dim, float>(j, name + "_j");
            }
        }

        convergenceCheck.setConvergence(allConverged);
    }

}

} // end of namespace
} // end of namespace

#endif // PALABOSRESULTWRITER_HH
