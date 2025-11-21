#include "configfactorybase.h"

namespace IPP
{

ConfigFactoryBase::~ConfigFactoryBase()
{

}

void ConfigFactoryBase::setCookie(const ConfigDataWrapper *configData)
{
    m_configData = configData;
}

}
