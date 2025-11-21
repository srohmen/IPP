#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include <string>
#include <vector>
#include <sstream>


namespace IPP
{

namespace StringUtils
{

    class StringBeginsWith
    {
    public:
        StringBeginsWith(const std::string& beginStr);

        bool operator()(const std::string& toCheck) const;

    private:
        const std::string& m_beginStr;
    };


    template<typename T>
    std::string toString(const std::vector<T>& vec)
    {
        std::stringstream ss;
        ss << "[ ";
        for(size_t i = 0; i < vec.size(); ++i)
        {
            ss << vec[i];
            if(i+1 < vec.size())
            {
                ss << ", ";
            }
        }

        ss << " ]";
        return ss.str();
    }
}

}

#endif // STRINGUTILS_H
