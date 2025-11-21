#ifndef PHREEQCCONSTANTS_H
#define PHREEQCCONSTANTS_H

#include <string>

namespace IPP
{

namespace PhreeqcConstants
{

enum SelectedOutputBlock
{
    SO_PhasesAmount = 0,
    SO_ComponentTotal = 1,
    SO_PhasesSaturationIndices = 2,
    SO_NucleationMonomers = 3,
    SO_Aux = 4
};

static const bool s_calcH2OSeperately = true;
static const size_t s_primaryCompsBegin = s_calcH2OSeperately ? 4 : 3;

static const double s_minPoros = 0.0;

static const int s_unitSolution = 2;
static const int s_unitPPassemblage = 0;
static const int s_unitKinetics = 0;
static const int s_unitSSassemblage = 0;

static const std::string s_inertPhaseName = "inert";




}

}

#endif // PHREEQCCONSTANTS_H
