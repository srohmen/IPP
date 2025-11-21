#ifndef CEMHYD3DIMAGEREADER_H
#define CEMHYD3DIMAGEREADER_H

#include <string>
#include <vector>
#include <unordered_map>

#include "ippconfig.h"

namespace IPP
{

class ArrayDimensionConvert;

namespace Cemhyd3DImageReader
{

typedef std::unordered_map<size_t, IPP::Composition*> IDtoComposition;
typedef std::vector<IPP::Domain> DomainVec;

void read(const std::string& fileName, std::vector<size_t>& compIdArr);


void convert(const ArrayDimensionConvert& inputIndexConv,
             const IDtoComposition& idToComposition,
             const std::vector<size_t>& compIdArr,
             const IPPVector3DInt &origin,
             const IPPVector3DInt &bounds,
             DomainVec& domains);
}

} // end of namespace IPP

#endif // CEMHYD3DIMAGEREADER_H
