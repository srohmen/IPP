#ifndef FILESYSTEMUTILS_H
#define FILESYSTEMUTILS_H

#include <boost/filesystem/path.hpp>

namespace IPP
{

namespace FileSystemUtils
{

void createIfNotExisting(const boost::filesystem::path& outPath);

}

}

#endif // FILESYSTEMUTILS_H
