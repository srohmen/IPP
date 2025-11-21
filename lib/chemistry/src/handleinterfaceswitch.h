#ifndef HANDLEINTERFACESWITCH_H
#define HANDLEINTERFACESWITCH_H

#include <cstddef>
#include <vector>

namespace IPP
{

class HandleInterfaceSwitch
{
public:
    HandleInterfaceSwitch(std::vector<char>& isNonPermNonInterfaceCell,
                          std::vector<char>& isInterfaceCell);

    void setIsInterfaceRun(const bool isInterfaceRun);

    bool mustRun(const size_t iCell);

private:
    bool m_isInterfaceRun;
    std::vector<char>& m_isNonPermNonInterfaceCell;
    std::vector<char>& m_isInterfaceCell;
};


}


#endif // HANDLEINTERFACESWITCH_H
