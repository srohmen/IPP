#include "fixcellsidsfornewfielddecomposition.h"

#include <chrono>
#include <regex>

#include "fielddecomposition.h"
#include "ippexception.h"
#include "ippbox.h"
#include "ippstream.h"


namespace IPP
{
namespace FixCellsIDsForNewFieldDecomposition
{


static size_t findCoord(const std::vector<IPPBox3DLong>& domains,
                        const IPPVector3DLong& coord)

{
    static const size_t invalidId = -1;
    size_t result = invalidId;

    size_t id = 0;

    for(size_t iDomain = 0; result == invalidId && iDomain < domains.size(); ++iDomain)
    {
        const IPPBox3DLong& domain = domains[iDomain];

        if(domain.contains(coord))
        {
            for(long x = domain.lower[0]; result == invalidId && x <= domain.upper[0]; ++x)
            {
                for(long y = domain.lower[1]; result == invalidId && y <= domain.upper[1]; ++y)
                {
                    for(long z = domain.lower[2]; result == invalidId && z <= domain.upper[2]; ++z)
                    {
                        if(coord[0] == x && coord[1] == y && coord[2] == z)
                        {
                            result = id;
                        }

                        ++id;
                    }
                }
            }
        }

        id += domain.getDiagonalSize();
    }

    IPPCheck::assertCheck(result != invalidId);

    return result;
}

static void findOldToNewCellsIDs(const FieldDecomposition& oldDecomp,
                                 const FieldDecomposition& newDecomp,
                                 std::vector<size_t>& oldToNewCellsID)
{
    const std::vector<IPPBox3DLong>& oldDomains = oldDecomp.getDomains();
    const std::vector<IPPBox3DLong>& newDomains = newDecomp.getDomains();

    IPPCheck::assertCheck(oldDecomp.getGlobalNumberCells() == newDecomp.getGlobalNumberCells());
    oldToNewCellsID.resize(newDecomp.getGlobalNumberCells(), -1);

    size_t oldId = 0;

    for(size_t i = 0; i < oldDomains.size(); ++i)
    {
        const IPPBox3DLong& oldDomain = oldDomains[i];
        for(long x = oldDomain.lower[0]; x <= oldDomain.upper[0]; ++x)
        {
            for(long y = oldDomain.lower[1]; y <= oldDomain.upper[1]; ++y)
            {
                for(long z = oldDomain.lower[2]; z <= oldDomain.upper[2]; ++z)
                {
                    const IPPVector3DLong coord = {{ x,y,z }};
                    const size_t newId = findCoord(newDomains, coord);
                    oldToNewCellsID[oldId] = newId;
                    ++oldId;
                }
            }
        }
    }
}

static void fixCellIds(const std::vector<size_t>& oldToNewCellsID,
                       std::istream& input,
                       std::ostream& output)
{
    pcout << "fixing number cell IDs: " << oldToNewCellsID.size() << std::endl;

    auto last = std::chrono::steady_clock::now();
    size_t lineCounter = 0;

    const std::regex rawRx("_RAW\\s+\\d+");
    const std::regex digitRx("\\d+");


    for(std::string line; std::getline(input, line); )
    {
        std::smatch rawMatch;
        if(std::regex_search(line, rawMatch, rawRx))
        {
            // extract old ID
            std::smatch digitMatch;
            std::regex_search(line, digitMatch, digitRx);
            IPPCheck::assertCheck(digitMatch.size() == 1);

            const std::ssub_match base_sub_match = digitMatch[0];
            const std::string idStr = base_sub_match.str();
            const int oldId = std::stoi(idStr);
            const size_t newId = oldToNewCellsID.at(oldId);
            const std::string newIdStr = std::to_string(newId);

            // pcout << idStr << " -> " << newIdStr << std::endl;

            line = std::regex_replace(line, digitRx, newIdStr);

        }

        output << line << std::endl;


        const auto duration = std::chrono::duration_cast<std::chrono::seconds>
                              (std::chrono::steady_clock::now() - last);
        if(duration.count() >= 10)
        {
            std::cout << "processed lines: " << lineCounter << std::endl;
            last = std::chrono::steady_clock::now();
        }

        ++lineCounter;
    }
}

void fix(const IPP::FieldDecomposition& oldDecomp,
         const IPP::FieldDecomposition& newDecomp,
         std::istream& input,
         std::ostream& output)
{
    std::vector<size_t> oldToNewCellsID;
    findOldToNewCellsIDs(oldDecomp, newDecomp, oldToNewCellsID);
    fixCellIds(oldToNewCellsID, input, output);
}


}

}
