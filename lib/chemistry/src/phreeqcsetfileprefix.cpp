#include "phreeqcsetfileprefix.h"

#include "createfilename.h"
#include "mpidelegate.h"
#include "mpitools.h"

#include "localaccessphreeqcrm.h"

namespace IPP
{


PhreeqcSetFilePrefix::PhreeqcSetFilePrefix(LocalAccessPhreeqcRM &phreeqc,
                                           const boost::filesystem::path &chemDir,
                                           const size_t iteration)
    : m_phreeqc(phreeqc)
{
    const std::string fileName = CreateFileName::create("mainrun_", iteration);
    const boost::filesystem::path outPath = chemDir / fileName;
    m_filePrefix = outPath.string();
}

void PhreeqcSetFilePrefix::runLocal() const
{
    MPI_CHECK_SYNC;

    if(m_phreeqc.GetMpiMyself() == 0)
    {
        const std::string currPrefix = m_phreeqc.GetFilePrefix();

        if(currPrefix != m_filePrefix)
        {
            m_phreeqc.CloseFiles();
            m_phreeqc.setFilePrefixLocal(m_filePrefix);
            m_phreeqc.OpenFiles();
        }
    }

}


void PhreeqcSetFilePrefix::run() const
{
    assert(m_phreeqc.GetMpiMyself() == 0);

    const std::string currPrefix = m_phreeqc.GetFilePrefix();

    if(currPrefix != m_filePrefix)
    {
        m_phreeqc.CloseFiles();
        m_phreeqc.SetFilePrefix(m_filePrefix);
        m_phreeqc.OpenFiles();
    }
}


}
