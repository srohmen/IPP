#ifndef HANDLEDISABLEDREACTIONCELLS_H
#define HANDLEDISABLEDREACTIONCELLS_H

#include <assert.h>
#include <cstddef>

namespace IPP
{

namespace HandleDisabledReactionCells
{

    template<typename Container, typename Predicate>
    void updateEnabledValues(const Container& src, const Predicate& func, Container& dst)
    {
        assert(src.size() == dst.size());

        for(size_t i = 0; i < src.size(); ++i)
        {
            const bool isEnabled = func(i);
            if(isEnabled)
            {
                dst[i] = src[i];
            }
            else
            {
                // do not change old value
            }
        }
    }



}

}

#endif // HANDLEDISABLEDREACTIONCELLS_H
