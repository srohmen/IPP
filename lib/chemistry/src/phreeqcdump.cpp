#include "phreeqcdump.h"

#include <PhreeqcRM.h>

#include "phreeqcconstants.h"
#include "ippexception.h"

namespace IPP
{


void PhreeqcDump::dump(PhreeqcRM& phreeqc, const std::string& outDir)
{
    IRM_RESULT status;

    status = phreeqc.SetDumpFileName(outDir);
    IPPCheck::assertCheck(status == IRM_OK);

    status = phreeqc.DumpModule(true, false);
    IPPCheck::assertCheck(status == IRM_OK);
}


}
