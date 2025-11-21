#ifndef VTKIMAGEREADER_H
#define VTKIMAGEREADER_H

#include <string>
#include <vector>

#include "ippvector.h"

namespace IPP
{

class IPPConfig;
class Composition;

namespace VTKImageReader
{

void read(const std::string& vtkFileName, std::vector<size_t>& arr);


void readDegrVTKs(const std::string& phasesFileName,
                  const std::string& compFileName,
                  const std::string &inerFracMap,
                  Composition compTemplate,
                  const IPPVector3DInt& origin,
                  IPPConfig& config);

}


}

#endif // VTKIMAGEREADER_H
