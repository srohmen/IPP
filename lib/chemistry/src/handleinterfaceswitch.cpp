#include "handleinterfaceswitch.h"

#include <cassert>

namespace IPP
{

HandleInterfaceSwitch::HandleInterfaceSwitch(std::vector<char>& isNonPermNonInterfaceCell, std::vector<char>& isInterfaceCell)
    : m_isInterfaceRun(false)
    , m_isNonPermNonInterfaceCell(isNonPermNonInterfaceCell)
    , m_isInterfaceCell(isInterfaceCell)
{

}

void HandleInterfaceSwitch::setIsInterfaceRun(const bool isInterfaceRun)
{
    m_isInterfaceRun = isInterfaceRun;
}

bool HandleInterfaceSwitch::mustRun(const size_t iCell)
{
    assert(m_isNonPermNonInterfaceCell.size() > iCell);

    if(m_isNonPermNonInterfaceCell[iCell])
    {
        return false;
    }

    assert(m_isInterfaceCell.size() > iCell);
    if(m_isInterfaceRun == m_isInterfaceCell[iCell])
    {
        return true;
    }
    else
    {
        return false;
    }

}



}
