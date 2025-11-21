#ifndef CONFIGDATAWRAPPER_H
#define CONFIGDATAWRAPPER_H

#include <json/value.h>

namespace IPP
{

struct ConfigDataWrapper
{
    ConfigDataWrapper(const Json::Value& node)
        : node(node)
    {

    }

    const Json::Value& node;
};

} // end of namespace IPP

#endif // CONFIGDATAWRAPPER_H
