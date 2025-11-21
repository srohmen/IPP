#include "retrievephaseamounts.h"

#include <chrono>
#include <cmath>

#include "localaccessphreeqcrm.h"

#include "ippexception.h"
#include "ippconstants.h"
#include "phreeqcconstants.h"
#include "applyfuncontupleelement.h"
#include "stringutils.h"
#include "isenabledincontainer.h"
#include "handledisabledreactioncells.h"

#include "ippstream.h"

namespace IPP
{

RetrievePhaseAmounts::RetrievePhaseAmounts(LocalAccessPhreeqcRM& phreeqc, const std::vector<char>& enabledCells)
    : m_phreeqc(phreeqc),
      m_enabledCells(enabledCells)
{

}

void RetrievePhaseAmounts::retrieve(std::vector<double> &precipMoles) const
{
    //    int n_user = m_phreeqcRM.GetNthSelectedOutputUserNumber(SO_Phases);
    //    (void) n_user;
    IRM_RESULT status = m_phreeqc.SetCurrentSelectedOutputUserNumber(PhreeqcConstants::SO_PhasesAmount);
    IPPCheck::assertCheck(status == IRM_OK);

    const size_t nCells = m_enabledCells.size();
    const size_t nCols = m_phreeqc.GetSelectedOutputColumnCount();

    // equi phases have differences too, thus multiply them by 2.
    // this does not calculate correctly when solid solutions are involved. disabling check for now
    // assert(nCols == 2*precipMoles.size() / nCells);

    // phreeqc does give global row count instead of local
    // const size_t nRows = m_phreeqc.GetSelectedOutputRowCount();
    // assert(nCells == nRows);
    assert(nCells == m_phreeqc.getNumberCells());


    std::vector<double> so;
    {
        status = m_phreeqc.getSelectedOutputLocal(so);
        IPPCheck::assertCheck(status == IRM_OK);

        // this does not calculate correctly when solid solutions are involved. disabling check for now
        // assert(so.size() == 2*precipMoles.size());

        //        if(so.size() != 2*precipMoles.size())
        //        {
        //            int* foo = nullptr;
        //            int bar = *foo;
        //            std::cout << bar;
        //        }
    }

    // each second col is the difference to prev timestep
    typedef std::vector<std::pair<std::string, size_t>> SrcMap;
    SrcMap srcMap;
    for(size_t iCol = 0; iCol < nCols; ++iCol)
    {
        std::string heading;
        status = m_phreeqc.GetSelectedOutputHeading(iCol, heading);
        IPPCheck::assertCheck(status == IRM_OK);
        srcMap.push_back(std::make_pair(heading, iCol));
    }

    // remove difference columns
    typedef ApplyFuncOnTupleElement<StringUtils::StringBeginsWith, 0> StringBeginsWithIndexed;
    srcMap.erase(std::remove_if(srcMap.begin(), srcMap.end(), StringBeginsWithIndexed(StringUtils::StringBeginsWith("d_"))),
                 srcMap.end());

    srcMap.erase(std::remove_if(srcMap.begin(), srcMap.end(), StringBeginsWithIndexed(StringUtils::StringBeginsWith("dk_"))),
                 srcMap.end());

    // remove k_ from kinetic reactants
    for(std::pair<std::string, size_t>& nameToIndex : srcMap)
    {
        std::string& phaseName = nameToIndex.first;
        if(phaseName[0] == 'k' && phaseName[1] == '_')
        {
            phaseName = phaseName.substr(2);
        }
    }

    // remove s_ from solid solutions
    for(std::pair<std::string, size_t>& nameToIndex : srcMap)
    {
        std::string& phaseName = nameToIndex.first;
        if(phaseName[0] == 's' && phaseName[1] == '_')
        {
            phaseName = phaseName.substr(2);
        }
    }

    std::sort(srcMap.begin(), srcMap.end());

    const size_t nPhasesTimesCells = nCells * srcMap.size();
    std::vector<double> newPrecipMoles(nPhasesTimesCells);

    for(size_t i = 0; i < srcMap.size(); ++i)
    {
        const size_t srcCol = srcMap[i].second;
        const size_t srcOffset = nCells * srcCol;
        const size_t dstOffset = nCells * i;

        std::copy(so.begin() + srcOffset,
                  so.begin() + srcOffset + nCells,
                  newPrecipMoles.begin() + dstOffset);
    }




    for(size_t i = 0; i < newPrecipMoles.size(); ++i)
    {
        if(newPrecipMoles[i] < 0.0)
        {
            const double tolerance = 1.0E-7;

            std::ios::fmtflags f( std::cout.flags() );

            std::cout << "rank " << MPIManager::getInstance().getRank()
                      << ": strange negative precip mole occured. index: " << i
                      << "\tiCellLocal: " << i % nCells
                      << "\tiPhase: " << i / nCells
                      << "\t" << std::scientific << newPrecipMoles[i]
                      << "\t" << std::hex << reinterpret_cast<int64_t&>(newPrecipMoles[i])
                      << "\t isEnabled: " << (int)m_enabledCells[i % nCells]
                      << "\t snapping to zero if within tolerance ("
                      << std::scientific << tolerance << ")" << std::endl;

            std::cout.flags( f );

            if(std::abs(newPrecipMoles[i]) < tolerance )
            {
                newPrecipMoles[i] = 0.0;
            }
            else
            {
                throw std::runtime_error("Phreeqc is broken! Error larger than "
                                         + std::to_string(tolerance));
            }

        }
    }

    precipMoles.resize(nPhasesTimesCells);
    typedef IsEnabledInContainer<decltype(m_enabledCells)> IsEnabledFunc;
    typedef WrapIndexFunc<IsEnabledFunc> WrapFunc;
    HandleDisabledReactionCells::updateEnabledValues(newPrecipMoles,
                                                     WrapFunc(IsEnabledFunc(m_enabledCells),
                                                              nCells), precipMoles);
}

}
