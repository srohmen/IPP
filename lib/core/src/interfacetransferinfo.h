#ifndef INTERFACETRANSFERINFO_H
#define INTERFACETRANSFERINFO_H

#include <vector>

namespace IPP
{

struct TransferInfo
{
    size_t rank;

};

struct InterfaceTransferInfo
{
    TransferInfo send;
    TransferInfo recieve;
};

typedef std::vector<InterfaceTransferInfo> InterfaceTransferInfo;

}

#endif // INTERFACETRANSFERINFO_H
