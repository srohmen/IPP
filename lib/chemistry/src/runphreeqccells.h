#ifndef RUNPHREEQCCELLS_H
#define RUNPHREEQCCELLS_H

#include <cstddef>
#include <vector>

#include <IrmResult.h>


namespace IPP
{

class LocalAccessPhreeqcRM;
class PhreeqcSetFilePrefix;
class SetupPhreeqcVerbosity;

class RunPhreeqcCells
{
public:
    RunPhreeqcCells(LocalAccessPhreeqcRM& phreeqc,
                    const PhreeqcSetFilePrefix &setPrefix,
                    const SetupPhreeqcVerbosity& setupVerbosity);

    void run() const;

private:
    void checkConvergence(IRM_RESULT status) const;

    LocalAccessPhreeqcRM& m_phreeqc;
    const PhreeqcSetFilePrefix& m_setPrefix;
    const SetupPhreeqcVerbosity& m_setupVerbosity;
};

}

#endif // RUNPHREEQCCELLS_H
