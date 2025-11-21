#ifndef PLOTTER2D_H
#define PLOTTER2D_H

#include <vector>
#include <memory>

class vtkTable;
class vtkFloatArray;
class vtkContextView;
class vtkChartXY;

namespace IPP
{

struct TableData2D
{
    TableData2D(const size_t nCurves, const size_t nx, const double dx)
        : nCurves(nCurves)
        , nx(nx)
        , dx(dx)
    {

    }

    const size_t nCurves;
    const size_t nx;
    const double dx;
};


class Plotter2DImpl;

class Plotter2D
{
public:
    Plotter2D();
    ~Plotter2D();
    void init();
    void initPointData(const TableData2D& pointData);
    void initLineData(const TableData2D& lineData);

    void updatePointData(const size_t iCol, const std::vector<double>& values);
    void updateLineData(const size_t iCol, const std::vector<double>& values);


    void show(const bool stay);


private:
    std::unique_ptr<Plotter2DImpl> m_impl;


};

}

#endif // PLOTTER2D_H
