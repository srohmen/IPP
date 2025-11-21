#ifndef ISINSIDEVOLUME_H
#define ISINSIDEVOLUME_H

#include <cassert>
#include <cstddef>
#include <vector>

namespace IPP
{

class IsInsideVolume
{
public:
    IsInsideVolume(const std::vector<double>& distField);

    inline bool evaluate(const size_t iCell) const
    {
        assert(iCell < m_distField.size());
        const double& latticeDistance = m_distField[iCell];

        if(latticeDistance < -1.0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

private:
    const std::vector<double>& m_distField;
};

}


#endif // ISINSIDEVOLUME_H
