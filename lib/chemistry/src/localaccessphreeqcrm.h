#ifndef LOCALACCESSPHREEQCRM_H
#define LOCALACCESSPHREEQCRM_H

#include <PhreeqcRM.h>

namespace IPP
{

class FieldDecomposition;

class LocalAccessPhreeqcRM : public PhreeqcRM
{
public:
    LocalAccessPhreeqcRM(int nxyz_arg, MP_TYPE parallelProcessData, PHRQ_io *io=NULL);
    ~LocalAccessPhreeqcRM();

    void init(const FieldDecomposition* decomp);

    void setGlobalEnabled(const bool state);

    void getEnabledCells(std::vector<char>& enabledCells) const;

    void enableCell(const size_t iCellLocal);
    void disableCell(const size_t iCellLocal);

    IRM_RESULT runCells();
    IRM_RESULT runCellsIfEnabled();
    IRM_RESULT runCellsLocal();

    const std::vector<double>& getPorosityLocal() const;
    void setPorosityLocal(const std::vector<double>& porosities);

    void getSolutionVolumeLocal(std::vector<double>& volume) const;

    IRM_RESULT getConcentrationsLocal(std::vector<double> &c);
    IRM_RESULT setConcentrationsLocal(const std::vector<double> &c);

    IRM_RESULT getSelectedOutputLocal(std::vector<double> &so);

    void setTimeLocal(const double& t);    
    void setFilePrefixLocal(const std::string& prefix);
    void setPrintChemistryMaskLocal(std::vector<int> & m);
    void setPrintChemistryOnLocal(bool worker, bool ip, bool utility);
    IRM_RESULT runStringLocal(bool runWorkers, bool runInitialPhreeqc,
                              bool runUtility, const std::string& input);

    IRM_RESULT runStringLocal(bool runWorkers, bool runInitialPhreeqc,
                              bool runUtility, std::istream& input);


    size_t getNumberCells() const;



private:
    void applyDecomp(const FieldDecomposition& decomp);
    void enableAllCellsTemp();
    void resetToEnabledStatus();

    void prepareErrorStreams(std::vector<int>& r_vector);
    void runCellsLocal(std::vector<int> &r_vector);
    IRM_RESULT runStringStreamThread(int n, std::istream& input);
    IRM_RESULT retrieveErrorsLocal(std::vector<int>& r_vector);
    IRM_RESULT retrieveErrors(std::vector<int>& r_vector);
    IRM_RESULT collectErrors(const std::vector<IRM_RESULT>& errors);

    bool m_enabled;

    const FieldDecomposition* m_decomp;

    bool m_isAllEnabled;
    std::vector<double> m_saturationCache;
    std::vector<double> m_saturationEnabled;

    std::vector<double> c_chem;
    std::vector<double> m_solns;


};

}

#endif // LOCALACCESSPHREEQCRM_H
