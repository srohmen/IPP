#include "inputgeometryscale.h"
#include "ippexception.h"
#include <cmath>

namespace IPP
{

namespace InputGeometryScale
{


template <typename T>
bool isInRange(const T& value, const T& low, const T& high)
{
    return !(value < low) && (value < high);
}

const size_t srcBoundarySize = 1;
const size_t dstBoundarySize = srcBoundarySize * 2;
const size_t dstBulkSize = 3;

inline size_t calcBulkDstBegin(const size_t src)
{
    return (src - 1) * dstBulkSize + dstBoundarySize;
}

void srcToDst(const ColumnMajorArrayDimensionConvert& dstIndexConv,
              const size_t dstXBegin, const size_t dstYBegin, const size_t dstZBegin,
              const size_t dstXEnd, const size_t dstYEnd, const size_t dstZEnd,
              const size_t srcIndex,
              const std::vector<size_t>& input,
              std::vector<size_t>& output)
{
    for(size_t dstX = dstXBegin; dstX < dstXEnd; ++dstX)
    {
        for(size_t dstY = dstYBegin; dstY < dstYEnd; ++dstY)
        {
            for(size_t dstZ = dstZBegin; dstZ < dstZEnd; ++dstZ)
            {
                const size_t dstIndex = dstIndexConv.calcIndex(dstX, dstY, dstZ);
                output[dstIndex] = input[srcIndex];
            }
        }
    }
}

void scaleBulk(const ColumnMajorArrayDimensionConvert& srcIndexConv,
               const ColumnMajorArrayDimensionConvert& dstIndexConv,
               const std::vector<size_t>& input, std::vector<size_t>& output)
{
    for(size_t srcX = 1; srcX+1 < srcIndexConv.getNx(); ++srcX)
    {
        const size_t dstXBegin = calcBulkDstBegin(srcX);
        const size_t dstXEnd = dstXBegin + dstBulkSize;

        for(size_t srcY = 1; srcY+1 < srcIndexConv.getNy(); ++srcY)
        {
            const size_t dstYBegin = calcBulkDstBegin(srcY);
            const size_t dstYEnd = dstYBegin + dstBulkSize;

            if(srcIndexConv.getNz() > 1)
            {
                for(size_t srcZ = 1; srcZ+1 < srcIndexConv.getNz(); ++srcZ)
                {
                    const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, srcZ);

                    const size_t dstZBegin = calcBulkDstBegin(srcZ);
                    const size_t dstZEnd = dstZBegin + dstBulkSize;

                    srcToDst(dstIndexConv,
                             dstXBegin, dstYBegin, dstZBegin,
                             dstXEnd, dstYEnd, dstZEnd,
                             srcIndex, input, output);
                }
            }
            else
            {
                const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);

                const size_t dstZBegin = 0;
                const size_t dstZEnd = 1;

                srcToDst(dstIndexConv,
                         dstXBegin, dstYBegin, dstZBegin,
                         dstXEnd, dstYEnd, dstZEnd,
                         srcIndex, input, output);
            }

        }

    }
}

void scaleEdges(const ColumnMajorArrayDimensionConvert& srcIndexConv,
                const ColumnMajorArrayDimensionConvert& dstIndexConv,
                const std::vector<size_t>& input, std::vector<size_t>& output)
{
    {
        for(size_t srcX = 1; srcX+1 < srcIndexConv.getNx(); ++srcX)
        {
            const size_t dstXBegin = calcBulkDstBegin(srcX);
            const size_t dstXEnd = dstXBegin + dstBulkSize;

            // bottom
            {
                const size_t srcY = 0;
                const size_t dstYBegin = 0;


                const size_t dstYEnd = dstYBegin + dstBoundarySize;

                IPPCheck::assertCheck(srcIndexConv.getNz() == 1, "TODO: implement 3d scaling");
                const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);

                const size_t dstZBegin = 0;
                const size_t dstZEnd = 1;

                srcToDst(dstIndexConv,
                         dstXBegin, dstYBegin, dstZBegin,
                         dstXEnd, dstYEnd, dstZEnd,
                         srcIndex, input, output);
            }

            // top
            {
                const size_t srcY = srcIndexConv.getNy() - 1;
                const size_t dstYBegin = calcBulkDstBegin(srcY);


                const size_t dstYEnd = dstYBegin + dstBoundarySize;

                IPPCheck::assertCheck(srcIndexConv.getNz() == 1, "TODO: implement 3d scaling");
                const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);

                const size_t dstZBegin = 0;
                const size_t dstZEnd = 1;

                srcToDst(dstIndexConv,
                         dstXBegin, dstYBegin, dstZBegin,
                         dstXEnd, dstYEnd, dstZEnd,
                         srcIndex, input, output);
            }
        }

        for(size_t srcY = 1; srcY+1 < srcIndexConv.getNy(); ++srcY)
        {
            const size_t dstYBegin = calcBulkDstBegin(srcY);
            const size_t dstYEnd = dstYBegin + dstBulkSize;

            // bottom
            {
                const size_t srcX = 0;
                const size_t dstXBegin = 0;


                const size_t dstXEnd = dstXBegin + dstBoundarySize;

                IPPCheck::assertCheck(srcIndexConv.getNz() == 1, "TODO: implement 3d scaling");
                const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);

                const size_t dstZBegin = 0;
                const size_t dstZEnd = 1;

                srcToDst(dstIndexConv,
                         dstXBegin, dstYBegin, dstZBegin,
                         dstXEnd, dstYEnd, dstZEnd,
                         srcIndex, input, output);
            }

            // top
            {
                const size_t srcX = srcIndexConv.getNx() - 1;
                const size_t dstXBegin = calcBulkDstBegin(srcX);


                const size_t dstXEnd = dstXBegin + dstBoundarySize;

                IPPCheck::assertCheck(srcIndexConv.getNz() == 1, "TODO: implement 3d scaling");
                const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);

                const size_t dstZBegin = 0;
                const size_t dstZEnd = 1;

                srcToDst(dstIndexConv,
                         dstXBegin, dstYBegin, dstZBegin,
                         dstXEnd, dstYEnd, dstZEnd,
                         srcIndex, input, output);
            }
        }
    }
}

void scaleCorner(const ColumnMajorArrayDimensionConvert& srcIndexConv,
                 const ColumnMajorArrayDimensionConvert& dstIndexConv,
                 const std::vector<size_t>& input,
                 std::vector<size_t>& output)
{
    {
        // lower left
        {
            const size_t srcX = 0;
            const size_t srcY = 0;

            const size_t dstXBegin = 0;
            const size_t dstYBegin = 0;

            const size_t dstXEnd = dstXBegin + dstBoundarySize;
            const size_t dstYEnd = dstYBegin + dstBoundarySize;

            IPPCheck::assertCheck(srcIndexConv.getNz() == 1, "TODO: implement 3d scaling");
            const size_t dstZBegin = 0;
            const size_t dstZEnd = 1;

            const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);
            srcToDst(dstIndexConv,
                     dstXBegin, dstYBegin, dstZBegin,
                     dstXEnd, dstYEnd, dstZEnd,
                     srcIndex, input, output);
        }

        // lower right
        {
            const size_t srcX = srcIndexConv.getNx() - 1;
            const size_t srcY = 0;

            const size_t dstXBegin = calcBulkDstBegin(srcX);
            const size_t dstYBegin = 0;

            const size_t dstXEnd = dstXBegin + dstBoundarySize;
            const size_t dstYEnd = dstYBegin + dstBoundarySize;

            IPPCheck::assertCheck(srcIndexConv.getNz() == 1, "TODO: implement 3d scaling");
            const size_t dstZBegin = 0;
            const size_t dstZEnd = 1;

            const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);
            srcToDst(dstIndexConv,
                     dstXBegin, dstYBegin, dstZBegin,
                     dstXEnd, dstYEnd, dstZEnd,
                     srcIndex, input, output);
        }

        // upper left
        {
            const size_t srcX = 0;
            const size_t srcY = srcIndexConv.getNy() - 1;

            const size_t dstXBegin = 0;
            const size_t dstYBegin = calcBulkDstBegin(srcY);

            const size_t dstXEnd = dstXBegin + dstBoundarySize;
            const size_t dstYEnd = dstYBegin + dstBoundarySize;

            IPPCheck::assertCheck(srcIndexConv.getNz() == 1, "TODO: implement 3d scaling");
            const size_t dstZBegin = 0;
            const size_t dstZEnd = 1;

            const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);
            srcToDst(dstIndexConv,
                     dstXBegin, dstYBegin, dstZBegin,
                     dstXEnd, dstYEnd, dstZEnd,
                     srcIndex, input, output);
        }

        // upper right
        {
            const size_t srcX = srcIndexConv.getNx() - 1;
            const size_t srcY = srcIndexConv.getNy() - 1;

            const size_t dstXBegin = calcBulkDstBegin(srcX);
            const size_t dstYBegin = calcBulkDstBegin(srcY);

            const size_t dstXEnd = dstXBegin + dstBoundarySize;
            const size_t dstYEnd = dstYBegin + dstBoundarySize;

            IPPCheck::assertCheck(srcIndexConv.getNz() == 1, "TODO: implement 3d scaling");
            const size_t dstZBegin = 0;
            const size_t dstZEnd = 1;

            const size_t srcIndex = srcIndexConv.calcIndex(srcX, srcY, 0);
            srcToDst(dstIndexConv,
                     dstXBegin, dstYBegin, dstZBegin,
                     dstXEnd, dstYEnd, dstZEnd,
                     srcIndex, input, output);
        }
    }
}

void scale(const ColumnMajorArrayDimensionConvert& srcIndexConv,
           const ColumnMajorArrayDimensionConvert& dstIndexConv,
           const std::vector<size_t>& input,
           std::vector<size_t>& output)
{
    const size_t nxyz = srcIndexConv.getNxyz();
    IPPCheck::assertCheck(nxyz == input.size(), "input geometry does not fit to scaled size");
    output.resize(dstIndexConv.getNxyz());

    // corner
    scaleCorner(srcIndexConv, dstIndexConv, input, output);

    // edges
    scaleEdges(srcIndexConv, dstIndexConv, input, output);

    // bulk
    scaleBulk(srcIndexConv, dstIndexConv, input, output);

}

size_t calcUpscaledSize(const size_t n, const size_t scalePasses)
{
    if(n == 1)
    {
        return 1;
    }
    else if(scalePasses == 0)
    {
        return n;
    }
    else
    {
        const size_t boundary = 2;
        const size_t inner = n - boundary;
        const size_t result = boundary * 2 + inner * 3;
        return calcUpscaledSize(result, scalePasses - 1);
    }
}

size_t calcDownscaledSize(const size_t n, const size_t scalePasses)
{
    if(n == 1)
    {
        return 1;
    }
    else if(scalePasses == 0)
    {
        return n;
    }
    else
    {
        const size_t boundary = 2;
        const size_t result = (n - (boundary * 2)) / 3 + boundary;
        return calcDownscaledSize(result, scalePasses - 1);
    }
}

}

}
