#ifndef PHASENAMETOINFOS_H
#define PHASENAMETOINFOS_H

#include <map>
#include <vector>
#include <string>

namespace IPP
{

struct PhaseInfo
{
    PhaseInfo()
        : molarVolume(-1.0)
        , stoich()
        , gfw(-1.0)
        , isPermeable(false)
    {

    }

    double molarVolume; // cm3/mol

    typedef std::vector<double> ElemIdToStoich;
    ElemIdToStoich stoich;

    double gfw; // g/mol
    bool isPermeable;



    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & molarVolume;
        ar & stoich;
        ar & gfw;
        ar & isPermeable;
    }
};

typedef std::map<std::string, PhaseInfo> PhaseNameToInfos;

}

#endif // PHASENAMETOINFOS_H
