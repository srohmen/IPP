#ifndef PALABOSFWD_H
#define PALABOSFWD_H

#include <palabos/dataProcessors/dataInitializerFunctional2D.h>
#include <palabos/dataProcessors/dataInitializerFunctional3D.h>

#include <palabos/atomicBlock/dataProcessingFunctional2D.h>
#include <palabos/atomicBlock/dataProcessingFunctional3D.h>

#include <palabos/atomicBlock/reductiveDataProcessingFunctional2D.h>
#include <palabos/atomicBlock/reductiveDataProcessingFunctional3D.h>



namespace plb
{

template<typename T>
class VtkImageOutput2D;


template<typename T>
class VtkImageOutput3D;

}

#endif // PALABOSFWD_H
