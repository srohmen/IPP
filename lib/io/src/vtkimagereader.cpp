#include "vtkimagereader.h"

#include <boost/filesystem.hpp>

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include "ippexception.h"

#include "ippconfig.h"
#include "arraydimensionconvert.h"
#include "ippstream.h"
#include "geometrytools.h"

namespace IPP
{

static vtkSmartPointer<vtkImageData> readVTK(const std::string &vtkFileName)
{
    vtkSmartPointer<vtkXMLImageDataReader> reader =
            vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(vtkFileName.c_str());
    reader->Update();

    return reader->GetOutput();
}

void VTKImageReader::read(const std::string &vtkFileName, std::vector<size_t> &arr)
{
    vtkSmartPointer<vtkXMLImageDataReader> reader =
            vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(vtkFileName.c_str());
    reader->Update();

    vtkSmartPointer<vtkImageData> image = reader->GetOutput();

    const int* dims = image->GetDimensions();
    const size_t nx = dims[0];
    const size_t ny = dims[1];
    const size_t nz = dims[2];

    arr.resize(nx*ny*nz);


    for(size_t z = 0; z < nz; ++z)
    {
        for(size_t y = 0; y < ny; ++y)
        {
            for(size_t x = 0; x < nx; ++x)
            {
                const void* rawVal = image->GetScalarPointer(x,y,z);
                const size_t val = *static_cast<const unsigned char*>(rawVal);

                const size_t index = z * nx * ny + y*nx + x;
                arr[index] = val;

            }
        }
    }
}

void VTKImageReader::readDegrVTKs(const std::string &phasesFileName,
                                  const std::string &compFileName,
                                  const std::string &inertFracMap,
                                  Composition compTemplate,
                                  const IPPVector3DInt &origin,
                                  IPPConfig &config)
{
    IPPCheck::assertCheck(boost::filesystem::exists(boost::filesystem::path(phasesFileName)),
                          "phase file not found: " + phasesFileName);

    IPPCheck::assertCheck(boost::filesystem::exists(boost::filesystem::path(compFileName)),
                          "comp file not found: " + compFileName);

    if(inertFracMap.empty() == false)
    {
        IPPCheck::assertCheck(boost::filesystem::exists(boost::filesystem::path(inertFracMap)),
                              "inertFracMap file not found: " + inertFracMap);
    }

    pcout << "reading: " << phasesFileName << "\t" << compFileName << "\t" << inertFracMap << std::endl;


    vtkSmartPointer<vtkImageData> phasesImage = readVTK(phasesFileName);
    vtkSmartPointer<vtkImageData> compImage = readVTK(compFileName);


    const int* dimsPhases = phasesImage->GetDimensions();
    {
        // all input VTK files must have the same dimensions
        const int* dimsComp = compImage->GetDimensions();
        IPPCheck::assertCheck(dimsPhases[0] == dimsComp[0]);
        IPPCheck::assertCheck(dimsPhases[1] == dimsComp[1]);
        IPPCheck::assertCheck(dimsPhases[2] == dimsComp[2]);
    }


    vtkPointData* phaseData = phasesImage->GetPointData();
    const int nPhases = phaseData->GetNumberOfArrays();
    for(int iArr = 0; iArr < nPhases; ++iArr)
    {
        const std::string phaseName = phaseData->GetArrayName(iArr);
        compTemplate.phases.emplace_back(phaseName, 0.0);
    }



    vtkPointData* compData = compImage->GetPointData();

    const int nComps = compData->GetNumberOfArrays();
    std::vector<size_t> srcToDst(nComps, -1);

    for(int iArr = 0; iArr < nComps; ++iArr)
    {
        const std::string compName = compData->GetArrayName(iArr);

        // exclude water components
        if(compName == "H2O" || compName == "H" || compName == "O" || compName == "Charge")
        {
            continue;
        }

        compTemplate.elemConcVec.emplace_back(compName, 0.0);
        srcToDst[iArr] = compTemplate.elemConcVec.size() - 1;
    }


    const size_t nxyz = config.nx * config.ny * config.nz;
    assert(config.compositions.empty());
    config.compositions.resize(nxyz, compTemplate);

    for(size_t iCell = 0; iCell < nxyz; ++iCell)
    {
        Composition& comp = config.compositions[iCell];
        comp.name = "cell_" + std::to_string(iCell);

        config.domains.emplace_back(&comp);
        Domain& domain = config.domains.back();
        domain.cells.resize(1, iCell);
    }

    const ArrayDimensionConvert dstConv(config.nx, config.ny, config.nz);
    const IPPVector3DInt vtkBounds = {dimsPhases[0], dimsPhases[1], dimsPhases[2]};

    if(inertFracMap.empty() == false)
    {
        vtkSmartPointer<vtkImageData> inertImage = readVTK(inertFracMap);

        const int* dimsInert = inertImage->GetDimensions();

        {
            // all input VTK files must have the same dimensions
            const int* dimsPhases = phasesImage->GetDimensions();
            IPPCheck::assertCheck(dimsPhases[0] == dimsInert[0], "inert map x dimension mismatch: "
                    + std::to_string(dimsPhases[0]) + " vs. " + std::to_string(dimsInert[0]));
            IPPCheck::assertCheck(dimsPhases[1] == dimsInert[1], "inert map y dimension mismatch: "
                    + std::to_string(dimsPhases[1]) + " vs. " + std::to_string(dimsInert[1]));
            IPPCheck::assertCheck(dimsPhases[2] == dimsInert[2], "inert map z dimension mismatch: "
                    + std::to_string(dimsPhases[2]) + " vs. " + std::to_string(dimsInert[2]));
        }

        for(size_t z = 0; z < config.nz; ++z)
        {
            for(size_t y = 0; y < config.ny; ++y)
            {
                for(size_t x = 0; x < config.nx; ++x)
                {
                    const IPPVector3DInt absPos =  {(int)x, (int)y, (int)z};
                    const IPPVector3DInt localPos = absPos + origin;

                    // just crop if its outside of domain
                    if(GeometryTools::isWithin(localPos, vtkBounds))
                    {
                        const void* rawVal = inertImage->GetScalarPointer(
                                    localPos[0], localPos[1], localPos[2]);
                        const double val = *static_cast<const double*>(rawVal);

                        const size_t index = dstConv.calcIndex(x, y, z);
                        config.compositions[index].inertVolFraction = val;
                    }
                    else
                    {

                    }
                }
            }
        }
    }
    else
    {
        for(size_t iCell = 0; iCell < nxyz; ++iCell)
        {
            Composition& comp = config.compositions[iCell];
            comp.inertVolFraction = 0.0;
        }
    }


    for(int iArr = 0; iArr < nPhases; ++iArr)
    {
        const std::string name = phaseData->GetArrayName(iArr);
        pcout << "reading data for: " << name << std::endl;
        vtkDoubleArray* scalars = dynamic_cast<vtkDoubleArray*>(phaseData->GetArray(name.c_str()));

        const int* dims = phasesImage->GetDimensions();

        for(size_t z = 0; z < config.nz; ++z)
        {
            for(size_t y = 0; y < config.ny; ++y)
            {
                for(size_t x = 0; x < config.nx; ++x)
                {
                    const IPPVector3DInt absPos =  {(int)x, (int)y, (int)z};
                    const IPPVector3DInt localPos = absPos + origin;

                    // just crop if its outside of domain
                    if(GeometryTools::isWithin(localPos, vtkBounds))
                    {
                        const size_t srcIndex =
                                localPos[0]
                                + localPos[1] * dims[0]
                                + localPos[2] * dims[0] * dims[1];

                        const double val = scalars->GetValue(srcIndex);

                        const size_t dstIndex = dstConv.calcIndex(x, y, z);
                        Composition& comp = config.compositions[dstIndex];
                        comp.phases[iArr].amount = val;
                    }
                }
            }
        }

    }

    for(int iArr = 0; iArr < nComps; ++iArr)
    {
        const size_t dstComp = srcToDst[iArr];
        const std::string name = compData->GetArrayName(iArr);
        if(dstComp == (size_t)-1)
        {
            pcout << "excluding: " << name << std::endl;
            continue;
        }

        pcout << "reading data for: " << name << std::endl;

        // palabos writes as float...
        vtkFloatArray* scalars = dynamic_cast<vtkFloatArray*>(compData->GetArray(name.c_str()));

        const int* dims = compImage->GetDimensions();

        for(size_t z = 0; z < config.nz; ++z)
        {
            for(size_t y = 0; y < config.ny; ++y)
            {
                for(size_t x = 0; x < config.nx; ++x)
                {
                    const IPPVector3DInt absPos =  {(int)x, (int)y, (int)z};
                    const IPPVector3DInt localPos = absPos + origin;

                    // just crop if its outside of domain
                    if(GeometryTools::isWithin(localPos, vtkBounds))
                    {
                        const size_t srcIndex =
                                localPos[0]
                                + localPos[1] * dims[0]
                                + localPos[2] * dims[0] * dims[1];

                        const double val = scalars->GetValue(srcIndex);

                        const size_t dstIndex = dstConv.calcIndex(x, y, z);
                        Composition& comp = config.compositions[dstIndex];
                        comp.elemConcVec[dstComp].conc = val;
                    }
                }
            }
        }

    }
}

} // end of namespace IPP







