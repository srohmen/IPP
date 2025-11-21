#include "nonreturningporositymodelfactory.h"

namespace IPP
{


NonReturningPorosityModelFactory::NonReturningPorosityModelFactory()
{

}

AbstractPorosityCalc* NonReturningPorosityModelFactory::generate() const
{
    return nullptr;
}


}
