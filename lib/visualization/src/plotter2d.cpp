#include "plotter2d.h"


#include "ippexception.h"

#ifndef NO_RENDERING
#include <vtkAutoInit.h>

VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
//VTK_MODULE_INIT(vtkRenderingContextOpenGL);
VTK_MODULE_INIT(vtkRenderingContextOpenGL2)
//VTK_MODULE_INIT(vtkContextDevice2D)

#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkPlotPoints.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>
#include <vtkAxis.h>
#include <vtkGlyph2D.h>

#include <vtkRenderingContextOpenGL2Module.h>

#include <vtkPen.h>



#endif

namespace IPP
{


class Plotter2DImpl
{
public:
    ~Plotter2DImpl() = default;
    void init();
    void initPointData(const TableData2D& pointData);
    void initLineData(const TableData2D& lineData);

    void updatePointData(const size_t iCol, const std::vector<double>& values);
    void updateLineData(const size_t iCol, const std::vector<double>& values);


    void show(const bool stay);


private:

    enum LineStyle
    {
        LS_Point,
        LS_Line
    };

#ifndef NO_RENDERING
    typedef vtkSmartPointer<vtkFloatArray> FloatArrPtr;

    void initTable(const TableData2D& data, const LineStyle &style,
                   std::vector<FloatArrPtr>& values, vtkSmartPointer<vtkTable> &table);
    auto toVtkStyle(const Plotter2DImpl::LineStyle &style);
    void addCurvesAll();
    void addCurves(const size_t nCurves, const vtkSmartPointer<vtkTable> &table, const LineStyle &style);

    vtkSmartPointer<vtkTable> m_pointDataTable;
    vtkSmartPointer<vtkTable> m_lineDataTable;

    FloatArrPtr arrXPoints;
    FloatArrPtr arrXLine;

    std::vector<FloatArrPtr> m_pointValues;
    std::vector<FloatArrPtr> m_lineValues;

    vtkSmartPointer<vtkContextView> m_view;
    vtkSmartPointer<vtkChartXY> m_chart;
#endif

};

#ifndef NO_RENDERING
auto Plotter2DImpl::toVtkStyle(const LineStyle& style)
{

    switch(style)
    {

    case LS_Point:
    {
        return vtkChart::POINTS;
        break;
    }

    case LS_Line:
    {
        return vtkChart::LINE;
        break;
    }

    default:
    {
        throw std::runtime_error("unknown line style");
        break;
    }


    }

}


void Plotter2DImpl::initTable(const TableData2D& data, const LineStyle& style,
                          std::vector<FloatArrPtr>& values, vtkSmartPointer<vtkTable>& table)
{

    for(size_t iCol = 0; iCol < data.nCurves; ++iCol)
    {
        FloatArrPtr col = vtkSmartPointer<vtkFloatArray>::New();
        values.push_back(col);
        col->SetName(std::to_string(iCol).c_str());
        table->AddColumn(col);
    }

    table->SetNumberOfRows(data.nx);
    for (size_t i = 0; i < data.nx; ++i)
    {
        table->SetValue(i, 0, i * data.dx);

        for(size_t iCol = 0; iCol < data.nCurves; ++iCol)
        {
            table->SetValue(i, iCol + 1, 0.0);
        }

    }

    // Add multiple scatter plots, setting the colors etc
    //            vtkPlot *points = chart->AddPlot(vtkChart::POINTS);
    //            points->SetInputData(table, 0, 1);
    //            points->SetColor(0, 0, 0, 255);
    //            points->SetWidth(1.0);
    //            vtkPlotPoints::SafeDownCast(points)->SetMarkerStyle(vtkPlotPoints::CROSS);
    addCurves(values.size(), table, style);


}


void Plotter2DImpl::addCurvesAll()
{
    m_view->GetScene()->ClearItems();
    m_view->GetScene()->AddItem(m_chart);

    addCurves(m_pointValues.size(), m_pointDataTable, LS_Point);
    addCurves(m_lineValues.size(), m_lineDataTable, LS_Line);
}

void Plotter2DImpl::addCurves(const size_t nCurves, const vtkSmartPointer<vtkTable>& table, const LineStyle& style)
{

    for(size_t iCurve = 0; iCurve < nCurves; ++iCurve)
    {
        const auto vtkStyle = toVtkStyle(style);
        vtkPlot *plot = m_chart->AddPlot(vtkStyle);
        // line->GetPen()->SetLineType(vtkPen::DASH_LINE);
        // line = chart->AddPlot(vtkChart::LINE);

        plot->SetInputData(table, 0, iCurve + 1);

        double rgb[3];
        vtkMath::HSVToRGB((double)iCurve/nCurves, 1.0, 0.7, &rgb[0], &rgb[1], &rgb[2]);
        plot->SetColor(rgb[0], rgb[1], rgb[2]);
        //                line->SetColor(255, 255, 0, 255);
        plot->SetWidth(3.0);
    }

}

#endif

void Plotter2DImpl::init()
{
#ifndef NO_RENDERING
    // Set up a 2D scene, add an XY chart to it
    m_view = vtkSmartPointer<vtkContextView>::New();
    m_view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
    m_view->GetRenderWindow()->SetSize(800, 800);

    m_chart = vtkSmartPointer<vtkChartXY>::New();
    m_view->GetScene()->AddItem(m_chart);
    m_chart->SetShowLegend(true);

    // FIXME: this is hardcoded range
    vtkAxis* axis = m_chart->GetAxis(vtkAxis::LEFT);
    axis->SetBehavior(vtkAxis::FIXED);
    axis->SetRange(-0.1, 2.0);



    //Finally render the scene
    m_view->GetRenderWindow()->SetMultiSamples(0);
    m_view->GetInteractor()->Initialize();
    // view->GetInteractor()->Start();


#endif
}


void Plotter2DImpl::initPointData(const TableData2D& pointData)
{
#ifndef NO_RENDERING
    m_pointDataTable = vtkSmartPointer<vtkTable>::New();
    arrXPoints = vtkSmartPointer<vtkFloatArray>::New();
    arrXPoints->SetName("x");
    m_pointDataTable->AddColumn(arrXPoints);
    initTable(pointData, LS_Point, m_pointValues, m_pointDataTable);
#endif
}

void Plotter2DImpl::initLineData(const TableData2D &lineData)
{
#ifndef NO_RENDERING
    m_lineDataTable = vtkSmartPointer<vtkTable>::New();
    arrXLine = vtkSmartPointer<vtkFloatArray>::New();
    arrXLine->SetName("x");
    m_lineDataTable->AddColumn(arrXLine);
    initTable(lineData, LS_Line, m_lineValues, m_lineDataTable);
#endif
}

void Plotter2DImpl::updatePointData(const size_t iCol, const std::vector<double>& values)
{
#ifndef NO_RENDERING
    for (size_t i = 0; i < values.size(); ++i)
    {
        const double& val = values[i];
        m_pointDataTable->SetValue(i, iCol+1, val);
    }
#endif
}

void Plotter2DImpl::updateLineData(const size_t iCol, const std::vector<double> &values)
{
#ifndef NO_RENDERING
    for (size_t i = 0; i < values.size(); ++i)
    {
        const double& val = values[i];
        m_lineDataTable->SetValue(i, iCol+1, val);
    }
#endif
}


void Plotter2DImpl::show(const bool stay)
{
#ifndef NO_RENDERING
    m_chart->ClearPlots();
    addCurvesAll();


    //view->GetRenderWindow()->SetMultiSamples(0);

    // Start interactor
    if(stay)
    {
        m_view->GetInteractor()->Start();
    }
    else
    {
        m_view->GetInteractor()->Render();
    }
#endif
}

Plotter2D::Plotter2D()
    : m_impl(new Plotter2DImpl)
{

}

Plotter2D::~Plotter2D() = default;

void Plotter2D::init()
{
    m_impl->init();
}

void Plotter2D::initPointData(const TableData2D &pointData)
{
    m_impl->initPointData(pointData);
}

void Plotter2D::initLineData(const TableData2D &lineData)
{
    m_impl->initLineData(lineData);
}

void Plotter2D::updatePointData(const size_t iCol, const std::vector<double> &values)
{
    m_impl->updatePointData(iCol, values);
}

void Plotter2D::updateLineData(const size_t iCol, const std::vector<double> &values)
{
    m_impl->updateLineData(iCol, values);
}

void Plotter2D::show(const bool stay)
{
    m_impl->show(stay);
}



}
