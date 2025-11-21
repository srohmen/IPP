#ifndef INDEXHELPER_H
#define INDEXHELPER_H

#include <cstddef>

namespace IPP
{

namespace IndexHelper
{

inline size_t toPhreeqcRawIndex(const size_t nCells, const size_t iCell, const size_t iComp)
{
    const size_t index = nCells * iComp + iCell;
    return index;
}

inline size_t toPhreeqcTransposedIndex(const size_t nComp, const size_t iCell, const size_t iComp)
{
    const size_t index = nComp * iCell + iComp;
    return index;
}


}

}

#endif // INDEXHELPER_H
