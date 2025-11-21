#ifndef NODEINFOS_H
#define NODEINFOS_H

#include <vector>

namespace IPP
{

struct NodeInfos
{
    std::vector<char> interfaceCells;
    std::vector<char> nonPermNonInterfaceCells;
    std::size_t nInterfaceNodesGlobal;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & interfaceCells;
        ar & nonPermNonInterfaceCells;
        ar & nInterfaceNodesGlobal;
    }
};

}

#endif // NODEINFOS_H
