#ifndef PHREEQCDUMP_H
#define PHREEQCDUMP_H

#include <string>

class PhreeqcRM;

namespace IPP
{

namespace PhreeqcDump
{

void dump(PhreeqcRM& phreeqc, const std::string& outDir);

}
}

#endif // PHREEQCDUMP_H
