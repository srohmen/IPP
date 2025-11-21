#include "retrieveauxvalues.h"

#include "localaccessphreeqcrm.h"
#include "mpitools.h"
#include "ippexception.h"
#include "phreeqcconstants.h"

namespace IPP
{


RetrieveAuxValues::RetrieveAuxValues(LocalAccessPhreeqcRM& phreeqc,
                                     const std::vector<char> &enabledCells)
    : m_phreeqc(phreeqc)
    , m_enabledCells(enabledCells)
{

}

void RetrieveAuxValues::retrieve(AuxDataVec& auxDataVec) const
{

    if(auxDataVec.empty() == false)
    {
        IRM_RESULT status;

        status = m_phreeqc.SetCurrentSelectedOutputUserNumber(PhreeqcConstants::SO_Aux);
        IPPCheck::assertCheck(status == IRM_OK);


        std::vector<double> so;
        status = m_phreeqc.getSelectedOutputLocal(so);
        IPPCheck::assertCheck(status == IRM_OK);


        const size_t nCellsPhreeqc = m_phreeqc.getNumberCells();

        for(size_t iDataSet = 0; iDataSet < auxDataVec.size(); ++iDataSet)
        {
            AuxDataName& auxData = auxDataVec[iDataSet];
            std::vector<double>& data = auxData.data;

            std::string header;
            status = m_phreeqc.GetSelectedOutputHeading(iDataSet, header);
            IPPCheck::assertCheck(status == IRM_OK);

            data.resize(nCellsPhreeqc, 0.0);


            assert(nCellsPhreeqc == m_enabledCells.size());
            const size_t soIndex = nCellsPhreeqc * iDataSet;
            std::vector<double>::const_iterator it = so.cbegin() + soIndex;

            for(size_t iCell = 0; iCell < nCellsPhreeqc; ++iCell)
            {
                assert(it != so.end());

                if(m_enabledCells[iCell])
                {
                    data[iCell] = *it;
                }

                ++it;
            }


        }
    }
}


}
