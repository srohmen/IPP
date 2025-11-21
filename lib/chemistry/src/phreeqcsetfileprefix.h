#ifndef PHREEQCSETFILEPREFIX_H
#define PHREEQCSETFILEPREFIX_H

#include <boost/filesystem/path.hpp>

namespace IPP
{

class LocalAccessPhreeqcRM;

class PhreeqcSetFilePrefix
{
public:
    PhreeqcSetFilePrefix(LocalAccessPhreeqcRM& m_phreeqc,
                         const boost::filesystem::path& chemDir,
                         const size_t iteration);

    void run() const;
    void runLocal() const;

private:
    LocalAccessPhreeqcRM& m_phreeqc;
    std::string m_filePrefix;
};



}


#endif // PHREEQCSETFILEPREFIX_H
