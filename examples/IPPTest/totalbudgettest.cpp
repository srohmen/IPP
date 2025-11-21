#include "totalbudgettest.h"

#include <numeric>
#include <cassert>

TotalBudgetTest::TotalBudgetTest()
{

}

void TotalBudgetTest::setPorosity(const std::vector<double> &porosity)
{
    m_porosity = porosity;
}

void TotalBudgetTest::begin(const size_t iteration, const double &time,
                            const std::vector<std::string>& dataNames)
{
    std::ios_base::openmode flag = std::ios_base::trunc;
    if(iteration > 0)
    {
        flag = std::ios_base::app;
    }

    assert(m_file.is_open() == false);
    m_file.open("budget", flag);

    if(iteration == 0)
    {
        m_file << "iteration,time,";
        for(size_t iData = 0; iData < dataNames.size(); ++iData)
        {
            const std::string& name = dataNames[iData];
            m_file << name << ",";
        }

        m_file << std::endl;
    }


    m_file << iteration << "," << time << ",";

}

void TotalBudgetTest::process(const std::string &name, const std::vector<double> &data)
{
    assert(m_file.is_open());

    std::vector<double> total(data.size());

    if(name == "HighVm")
    {
        total = data;
    }
    else
    {
        assert(m_porosity.size() == data.size());

        for(size_t i = 0; i < total.size(); ++i)
        {
            total[i] = m_porosity[i] * data[i];
        }
    }

    const double sum = std::accumulate(total.begin(), total.end(), 0.0);

    m_file << sum << ",";
}


void TotalBudgetTest::end()
{
    assert(m_file.is_open());

    m_file << std::endl;

    m_file.close();
}
