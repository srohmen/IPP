#include "nonreturningsicalcfactory.h"


namespace IPP
{


NonReturningSICalcFactory::NonReturningSICalcFactory()
{

}

AbstractSICalc *IPP::NonReturningSICalcFactory::generate() const
{
    return nullptr;
}

}
