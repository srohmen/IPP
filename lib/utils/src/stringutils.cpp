#include "stringutils.h"

namespace IPP
{
namespace StringUtils
{

StringBeginsWith::StringBeginsWith(const std::string& beginStr)
    : m_beginStr(beginStr)
{

}

bool StringBeginsWith::operator()(const std::string& toCheck) const
{
    const std::string subStr = toCheck.substr(0, m_beginStr.size());
    if(subStr == m_beginStr)
    {
        return true;
    }
    else
    {
        return false;
    }
}

}
}
