#ifndef JSONMULTISCALEDIFFUSIONFACTORY_H
#define JSONMULTISCALEDIFFUSIONFACTORY_H

#include "multiscalediffusionfactory.h"


namespace IPP
{

class JSONMultiScaleDiffusionFactory : public MultiScaleDiffusionFactory
{
public:
    JSONMultiScaleDiffusionFactory();    
    virtual ~JSONMultiScaleDiffusionFactory();

    virtual AbstractMultiScaleDiffusionCalc* generate() const;

private:

};

}

#endif // JSONMULTISCALEDIFFUSIONFACTORY_H
