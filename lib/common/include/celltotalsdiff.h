#ifndef CELLTOTALSDIFF_H
#define CELLTOTALSDIFF_H

#include <vector>

namespace IPP
{

struct CellTotalsDiff
{
    std::size_t iCell;
    std::vector<double> totalsDiff;


    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & iCell;
        ar & totalsDiff;
    }

};

}

#endif // CELLTOTALSDIFF_H
