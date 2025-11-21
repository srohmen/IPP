#ifndef SETUPPHREEQCVERBOSITY_H
#define SETUPPHREEQCVERBOSITY_H

namespace IPP
{

class LocalAccessPhreeqcRM;

class SetupPhreeqcVerbosity
{
public:
    SetupPhreeqcVerbosity(LocalAccessPhreeqcRM& phreeqc,
                          const bool isPrintChemistryOn);

    void runLocal() const;

    void run() const;

private:
    LocalAccessPhreeqcRM& m_phreeqc;
    const bool m_isPrintChemistryOn;
};

}

#endif // SETUPPHREEQCVERBOSITY_H
