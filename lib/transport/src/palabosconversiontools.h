#ifndef PALABOSCONVERSIONTOOLS_H
#define PALABOSCONVERSIONTOOLS_H

#include <palabos/core/geometry2D.h>
#include <palabos/core/geometry3D.h>

#include "ippvector.h"
#include "ippbox.h"
#include "array_view.h"

namespace IPP
{

namespace PalabosConversionTools
{

    template<size_t dim>
    struct ToPlb;

    template<>
    struct ToPlb<2>
    {
        template<typename Vector>
        static plb::Dot2D convertVector(const Vector& pos)
        {
            return plb::Dot2D(pos[0], pos[1]);
        }

        template<typename Box>
        static plb::Box2D convertBox(const Box& box)
        {
            return plb::Box2D(box[0], box[1], // x
                    box[2], box[3]); // y
        }

    };

    template<>
    struct ToPlb<3>
    {
        template<typename Vector>
        static plb::Dot3D convertVector(const Vector& pos)
        {
            return plb::Dot3D(pos[0], pos[1], pos[2]);
        }

        template<typename Box>
        static plb::Box3D convertBox(const Box& box)
        {
            return plb::Box3D(box[0], box[1], // x
                    box[2], box[3],  // y
                    box[4], box[5]); // z
        }
    };

    inline plb::Dot2D makeDot(const av::offset<2>& idx)
    {
        return plb::Dot2D(idx[0], idx[1]);
    }

    inline plb::Dot3D makeDot(const av::offset<3>& idx)
    {
        return plb::Dot3D(idx[0], idx[1], idx[2]);
    }

    template<typename A, typename B>
    inline const plb::Array<A, 3> toPlb(const Array<B, 3>& arr)
    {
        return plb::Array<double, 3>(arr[0], arr[1], arr[2]);
    }

    inline const plb::Array<double, 3> toPlb(const IPPVector3D& arr)
    {
        return plb::Array<double, 3>(arr[0], arr[1], arr[2]);
    }

    inline const plb::Array<double, 2> toPlb2D(const IPPVector3D& arr)
    {
        return plb::Array<double, 2>(arr[0], arr[1]);
    }

    inline const plb::Box2D toPlb(const IPPBox2DInt& box)
    {
        return plb::Box2D(box.lower[0], box.upper[0],
                box.lower[1], box.upper[1]);
    }

    inline const plb::Box3D toPlb(const IPPBox3DInt& box)
    {
        return plb::Box3D(box.lower[0], box.upper[0],
                box.lower[1], box.upper[1],
                box.lower[2], box.upper[2]);
    }

    inline const plb::Box2D toPlb2D(const IPPBox3DInt& box)
    {
        return plb::Box2D(box.lower[0], box.upper[0],
                box.lower[1], box.upper[1]);
    }


    template<typename Vector>
    inline auto toPlb(const Vector& vec, plb::Dot2D& dot)
    {
        dot = plb::Dot2D(vec[0], vec[1]);
    }

    template<typename Vector>
    inline auto toPlb(const Vector& vec, plb::Dot3D& dot)
    {
        dot = plb::Dot3D(vec[0], vec[1], vec[2]);
    }

    inline void convertBox(const IPPBox3DInt& box, plb::Box2D& result)
    {
        result = PalabosConversionTools::toPlb2D(box);
    }

    inline void convertBox(const IPPBox3DInt& box, plb::Box3D& result)
    {
        result = PalabosConversionTools::toPlb(box);
    }

    inline void convertBox(const plb::Box2D& box, IPPBox3DInt& result)
    {
        result[0] = box.x0;
        result[1] = box.x1;
        result[2] = box.y0;
        result[3] = box.y1;

        result[4] = 1;
        result[5] = 1;
    }

    inline void convertBox(const plb::Box3D& box, IPPBox3DInt& result)
    {
        result[0] = box.x0;
        result[1] = box.x1;
        result[2] = box.y0;
        result[3] = box.y1;
        result[4] = box.z0;
        result[5] = box.z1;

    }


    inline plb::Dot3D convertOffset(const plb::Dot2D& offset)
    {
        return plb::Dot3D(offset.x, offset.y, 0);
    }

    inline plb::Dot3D convertOffset(const plb::Dot3D& offset)
    {
        return offset;
    }

    template<typename T>
    inline void convertVector(const IPPVector3D& vec, plb::Array<T,2>& result)
    {
        result = toPlb2D(vec);
    }

    template<typename T>
    inline void convertVector(const IPPVector3D& vec, plb::Array<T,3>& result)
    {
        result = toPlb(vec);
    }

    template<typename S, size_t N, typename T>
    inline void convertVector(const Array<S, N>& vec, plb::Array<T,N>& result)
    {
        result = toPlb<T>(vec);
    }


    template<typename T>
    inline void convertVector(const plb::Array<T,3>& vec, plb::Array<T,2>& result)
    {
        result = plb::Array<T,2>(vec[0], vec[1]);
    }

    template<typename T>
    inline void convertVector(const plb::Array<T,3>& vec, plb::Array<T,3>& result)
    {
        result = vec;
    }


    inline const plb::Box2D toBox(const plb::Dot2D& dot)
    {
        return plb::Box2D(dot.x, dot.x, dot.y, dot.y);
    }

    inline const plb::Box3D toBox(const plb::Dot3D& dot)
    {
        return plb::Box3D(dot.x, dot.x, dot.y, dot.y, dot.z, dot.z);
    }

    template<typename Vector>
    inline void convertToIPP(const plb::Dot2D& dot, Vector& vec)
    {
        vec = {{ (typename Vector::value_type)dot.x,
                 (typename Vector::value_type)dot.y,
                 0 }};
    }

    template<typename Vector>
    inline void convertToIPP(const plb::Dot3D& dot, Vector& vec)
    {
        vec = {{ (typename Vector::value_type)dot.x,
                 (typename Vector::value_type)dot.y,
                 (typename Vector::value_type)dot.z }};
    }

    template<typename DotList, typename Vector>
    inline void convertTo3D(const DotList& dots, std::vector<Vector>& arrVec)
    {
        for(int i = 0; i < dots.getN(); ++i)
        {
            const auto& dot = dots.getDot(i);
            Vector vec;
            convertToIPP(dot, vec);
            arrVec.push_back(std::move(vec));
        }
    }


    inline plb::DotList2D makeDotList(const plb::Dot2D& pos)
    {
        plb::DotList2D list;
        list.addDot(pos);
        return list;
    }

    inline plb::DotList3D makeDotList(const plb::Dot3D& pos)
    {
        plb::DotList3D list;
        list.addDot(pos);
        return list;
    }


} // end of namespcae


}


#endif // PALABOSCONVERSIONTOOLS_H
