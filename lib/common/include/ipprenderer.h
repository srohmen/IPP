#ifndef IPPRENDERER_H
#define IPPRENDERER_H

namespace IPP
{

class AbstractScalarField;

class IPPRenderer
{
public:
    IPPRenderer()
    {

    }

    virtual ~IPPRenderer()
    {

    }

    virtual void updateImageData(const AbstractScalarField &imageSrc) = 0;

    virtual void render() = 0;

};

}

#endif // IPPRENDERER_H
