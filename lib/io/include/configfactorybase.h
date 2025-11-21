#ifndef CONFIGFACTORYBASE_H
#define CONFIGFACTORYBASE_H

namespace IPP
{

struct ConfigDataWrapper;

class ConfigFactoryBase
{
public:
    virtual ~ConfigFactoryBase();

    virtual void setCookie(const ConfigDataWrapper* configData);

protected:
    const ConfigDataWrapper* m_configData;
};

}

#endif // CONFIGFACTORYBASE_H
