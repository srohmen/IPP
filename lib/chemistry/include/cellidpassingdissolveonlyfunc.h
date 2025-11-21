#ifndef CELLIDPASSINGDISSOLVEONLYFUNC_H
#define CELLIDPASSINGDISSOLVEONLYFUNC_H

//#include "cellidpassingfuncwrapper.h"

#include "abstractdissolveonlycalcfwd.h"
#include "abstractdissolveonlycalc.h"
#include "cellid_dissolveonly_func.h"
#include "cellneighborinfo.h"
#include "phasenametoinfos.h"

namespace IPP
{

//typedef CellIDPassingFuncWrapper<CellID_DissolveOnly_Func, AbstractDissolveOnlyCalc> CellIDPassingDissolveOnlyFunc;

class CellIDPassingDissolveOnlyFunc : public CellID_DissolveOnly_Func
{
public:
    CellIDPassingDissolveOnlyFunc(AbstractDissolveOnlyCalcPtr dissOnlyFunc,
                                  const std::vector<double>& porosity,
                                  const std::vector<CellNeighborInfo>& neighInfos,
                                  const std::vector<double>& currentSI,
                                  const std::vector<size_t>& phaseToNuclPhaseID,
                                  const std::vector<size_t>& nuclPhaseToMonomerID,
                                  const std::vector<double>& monomerConc,
                                  const PhaseNameToInfos& phaseInfos);

    void init(const double& cellVol,
              const double& cellArea,
              const double& dt);

    virtual double getLowerPorosityThresh() const override;

    virtual bool needsNeighInfos() const override;
    virtual bool needsSaturationIndices() const override;

    virtual bool evaluate(const size_t iCell, const iterator& begin, const iterator& end) const override;

private:
    AbstractDissolveOnlyCalcPtr m_dissOnlyFunc;
    const std::vector<double>& m_porosity;
    const std::vector<CellNeighborInfo>& m_neighInfos;
    const std::vector<double>& m_currentSI;

    size_t m_nPhases;

    const std::vector<size_t>& m_phaseToNuclPhaseID;
    const std::vector<size_t>& m_nuclPhaseToMonomerID;
    const std::vector<double>& m_monomerConc;

    size_t m_nMonomers;

    const PhaseNameToInfos& m_phaseInfos;


};


}

#endif // CELLIDPASSINGDISSOLVEONLYFUNC_H
