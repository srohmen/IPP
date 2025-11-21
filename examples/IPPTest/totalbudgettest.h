#ifndef TOTALBUDGETTEST_H
#define TOTALBUDGETTEST_H

#include "auxresultprocessing.h"

#include <fstream>

class TotalBudgetTest : public IPP::AuxResultProcessing
{
public:
    TotalBudgetTest();

    virtual void setPorosity(const std::vector<double>& porosity) override;

    virtual void begin(const size_t iteration, const double& time,
                       const std::vector<std::string> &dataNames) override;

    virtual void process(const std::string& name,
                         const std::vector<double>& data) override;

    virtual void end() override;

private:
    std::ofstream m_file;
    std::vector<double> m_porosity;
};

#endif // TOTALBUDGETTEST_H
