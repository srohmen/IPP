#ifndef ABSTRACTSCALARFIELD_H
#define ABSTRACTSCALARFIELD_H

#include <cstddef>

#include "ippvector.h"

namespace IPP
{

class AbstractScalarField
{
public:
    AbstractScalarField()
    {

    }

    virtual ~AbstractScalarField()
    {

    }

    virtual size_t getNx() const = 0;
    virtual size_t getNy() const = 0;
    virtual size_t getNz() const = 0;

    double get(const IPPVector3DInt& pos) const
    {
        return this->get(pos[0], pos[1], pos[2]);
    }

    virtual double get(int x, int y, int z) const = 0;
    virtual void set(int x, int y, int z, const double& val) = 0;

};


class ConstAbstractScalarField
{
public:
    ConstAbstractScalarField()
    {

    }

    virtual ~ConstAbstractScalarField()
    {

    }

    virtual size_t getNx() const = 0;
    virtual size_t getNy() const = 0;
    virtual size_t getNz() const = 0;

    double get(const IPPVector3DInt& pos) const
    {
        return this->get(pos[0], pos[1], pos[2]);
    }

    virtual double get(int x, int y, int z) const = 0;

};


}


#endif // ABSTRACTSCALARFIELD_H
