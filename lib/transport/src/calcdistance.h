#ifndef CALCDISTANCE_H
#define CALCDISTANCE_H

#include <palabos/atomicBlock/dataProcessingFunctional2D.h>
#include <palabos/atomicBlock/atomicBlock2D.h>
#include <palabos/atomicBlock/dataProcessingFunctional3D.h>
#include <palabos/atomicBlock/atomicBlock3D.h>

#include <atomicBlock/dataField2D.h>


#include "distancetransform8ssedt_2D.h"
#include "distancetransform8ssedt_3D.h"

#include "scalarfieldwrapper.h"
#include "palabosfieldwrapper.h"
#include "rawscalarfield.h"
#include "virtualpalabosfieldwrapper.h"
#include "porosityaccess.h"
#include "arrayviewscalarfield.h"


#include "bench_tools.h"

namespace IPP
{

namespace CalcDistance
{

template<size_t dim>
struct GetDistanceTransform;

template<>
struct GetDistanceTransform<2>
{
    using type = DistanceTransform8SSEDT_2D;
};

template<>
struct GetDistanceTransform<3>
{
    using type = DistanceTransform8SSEDT_3D;
};


template<size_t dim, typename T>
void calcDistanceField(const av::array_view<T, dim>& input, const T& thresh, av::array_view<T, dim>& output)
{
    using DistTransform = typename GetDistanceTransform<dim>::type;
    DistTransform transform;

    const ArrayViewScalarField<av::array_view<T, dim>> in(input);
    ArrayViewScalarField<av::array_view<T, dim>> out(output);
    transform.calc(in, thresh, out);
}

template<typename T, size_t dim>
class CalcDistanceT;

template<typename T>
class CalcDistanceT<T, 2> : public plb::BoxProcessingFunctional2D_SS<T, T>
{
public:
    CalcDistanceT(const T& threshold)
        : m_thresh(threshold)
    {

    }

    virtual CalcDistanceT<T, 2>* clone() const
    {
        return new CalcDistanceT<T, 2>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::staticVariables;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return plb::BlockDomain::bulk;
    }



    virtual void process(plb::Box2D domain,
                         plb::ScalarField2D<T>& porosField,
                         plb::ScalarField2D<T>& distanceField)
    {
        // simpleNeighbor(porosField, domain, distanceField);
        BENCH(signedDistanceField(porosField, domain, distanceField));
    }

private:    
    void simpleNeighbor(plb::ScalarField2D<T>& porosField, const plb::Box2D& domain,
                        plb::ScalarField2D<T>& distanceField)
    {
        const plb::Dot2D offset = plb::computeRelativeDisplacement(porosField,distanceField);
        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                int numNeighbors = 0;
                for (plb::plint dx=-1; dx<=+1; ++dx)
                {
                    for (plb::plint dy=-1; dy<=+1; ++dy)
                    {
                        const plb::plint x = iX+dx;
                        const plb::plint y = iY+dy;

                        const T& porosity = porosField.get(x, y);
                        if(porosity <= m_thresh)
                        {
                            ++numNeighbors;
                        }

                    }
                }

                T& result = distanceField.get(iX+offset.x,iY+offset.y);
                if (numNeighbors > 0)
                {
                    result = 0.0;
                }
                else
                {
                    result = 2.0;
                }

            }
        }
    }

    void signedDistanceField(plb::ScalarField2D<T>& porosField, const plb::Box2D& domain,
                             plb::ScalarField2D<T>& distanceField)
    {

        VirtualPalabosFieldWrapper<T, plb::ScalarField2D<T>, 2> input(porosField);

        RawScalarField2D<T> tmpOutput(porosField.getNx(), porosField.getNy());
        RawScalarFieldWrapper2D<T> output(tmpOutput);

        DistanceTransform8SSEDT_2D transform;
        transform.calc(input, m_thresh, output);


        const plb::Dot2D offset = plb::computeRelativeDisplacement(porosField,distanceField);

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                T& result = distanceField.get(iX+offset.x,iY+offset.y);
                result = tmpOutput[iX][iY];
            }
        }

    }


    const T m_thresh;
};


template<typename T>
class CalcDistanceT<T, 3> : public plb::BoxProcessingFunctional3D_SS<T,T>
{
public:
    CalcDistanceT(const T& threshold)
        : m_thresh(threshold)
    {

    }

    virtual CalcDistanceT<T, 3>* clone() const
    {
        return new CalcDistanceT<T, 3>(*this);
    }

    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const
    {
        modified[0] = plb::modif::nothing;
        modified[1] = plb::modif::staticVariables;
    }

    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return plb::BlockDomain::bulk;
    }

    virtual void process(plb::Box3D domain,
                         plb::ScalarField3D<T>& porosField,
                         plb::ScalarField3D<T>& distanceField)
    {
        signedDistanceField(porosField, domain, distanceField);
    }

private:

    void signedDistanceField(plb::ScalarField3D<T>& porosField, const plb::Box3D& domain,
                             plb::ScalarField3D<T>& distanceField)
    {

        VirtualPalabosFieldWrapper<T, plb::ScalarField3D<T>, 3> input(porosField);

        RawScalarField3D<T> tmpOutput(porosField.getNx(), porosField.getNy());
        RawScalarFieldWrapper3D<T> output(tmpOutput);


        DistanceTransform8SSEDT_3D transform;
        transform.calc(input, m_thresh, output);


        const plb::Dot3D offset = plb::computeRelativeDisplacement(porosField,distanceField);

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                for (plb::plint iZ=domain.z0; iY<=domain.z1; ++iZ)
                {
                    T& result = distanceField.get(iX+offset.x, iY+offset.y, iZ + offset.z);
                    result = tmpOutput[iX][iY][iZ];
                }
            }
        }

    }


    const T m_thresh;
};

}


}


#endif // CALCDISTANCE_H
