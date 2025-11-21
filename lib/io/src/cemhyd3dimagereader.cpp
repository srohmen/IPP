#include "cemhyd3dimagereader.h"

#include <fstream>
#include <stdexcept>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include <boost/tokenizer.hpp>
#include <boost/range.hpp>

#include "ippexception.h"
#include "arraydimensionconvert.h"
#include "geometrytools.h"

namespace IPP
{

//#define POROSITY 0
//#define C3S 1
//#define C2S 2
//#define C3A 3
//#define C4AF 4
//#define GYPSUM 5
//#define HEMIHYD 6
//#define ANHYDRITE 7
//#define POZZ 8
//#define INERT 9
//#define SLAG 10
//#define ASG 11
//#define CAS2 12
//#define CH 13
//#define CSH 14
//#define C3AH6 15
//#define ETTR 16
//#define ETTRC4AF 17
//#define AFM 18
//#define FH3 19
//#define POZZCSH 20
//#define SLAGCSH 21
//#define CACL2 22
//#define FREIDEL 23
//#define STRAT 24
//#define GYPSUMS 25
//#define CACO3 26
//#define AFMC 27

//static const size_t s_maxID = 50;

void Cemhyd3DImageReader::read(const std::string& fileName, std::vector<size_t>& compIdArr)
{
    std::ifstream is(fileName);

    if(is.is_open() == false)
    {
        throw std::runtime_error("Could not open file: " + fileName);
    }

    // get length of file:
    is.seekg (0, is.end);
    int length = is.tellg();
    is.seekg (0, is.beg);

    std::vector<char> buffer(length);


    // read data as a block:
    is.read(buffer.data(), length);
    if(!is)
    {
        throw std::runtime_error("Reading file failed: " + fileName);
    }
    is.close();


    boost::char_separator<char> sep("\n");
    boost::tokenizer<boost::char_separator<char>, std::vector<char>::const_iterator > tokens(buffer, sep);

    compIdArr.clear();
    for(auto beg = tokens.begin(); beg != tokens.end(); ++beg)
    {
        const size_t compID = std::stoi(*beg);
        compIdArr.push_back(compID);
    }

}


void Cemhyd3DImageReader::convert(const ArrayDimensionConvert& inputIndexConv,
                                  const IDtoComposition& idToComposition,
                                  const std::vector<size_t>& compIdArr,
                                  const IPPVector3DInt& origin,
                                  const IPPVector3DInt& bounds,
                                  DomainVec& domains)
{
    typedef std::unordered_multimap<size_t, size_t> CompToCellMultiMap;
    typedef std::unordered_multimap<size_t, size_t>::const_iterator CompToCellMultiMapIterator;
    CompToCellMultiMap compToCellMap;
    for(size_t i = 0; i < compIdArr.size(); ++i)
    {
        const size_t compID = compIdArr[i];
        compToCellMap.insert(std::make_pair(compID, i));
    }




    const ArrayDimensionConvert outputIndexConv(bounds);

    for(CompToCellMultiMapIterator cemHydIDDomain = compToCellMap.begin();
        cemHydIDDomain != compToCellMap.end(); ++cemHydIDDomain)
    {
        const size_t compID = cemHydIDDomain->first;

        const IDtoComposition::const_iterator compIt = idToComposition.find(compID);

        if(compIt == idToComposition.end())
        {
            throw std::runtime_error("Composition not defined. ID: " + std::to_string(compID));
        }

        const Composition* comp = compIt->second;

        domains.push_back(Domain(comp));
        Domain& domain = domains.back();

        //        std::cout << compID << "\t" << comp->name << std::endl;

        CompToCellMultiMapIterator innerIt = cemHydIDDomain;
        CompToCellMultiMapIterator prevIt = innerIt;
        while(innerIt != compToCellMap.end() && innerIt->first == compID)
        {
            const size_t colMajorCellID = innerIt->second;
            const IPPVector3DInt pos = inputIndexConv.colMajorToRowMajorCoord(colMajorCellID);
            const IPPVector3DInt absPos = pos + origin;

            if(GeometryTools::isWithin(absPos, bounds))
            {
                const size_t index = outputIndexConv.calcIndex(absPos);
                domain.cells.push_back(index);
            }
            else
            {
//                std::cout << "excluding: " << absPos << std::endl;
            }


            prevIt = innerIt;
            ++innerIt;
        }

        cemHydIDDomain = prevIt;
    }
}



} // end of namespace LBGeoChem
