#include "filesystemutils.h"

#include <boost/filesystem.hpp>

namespace IPP
{

void FileSystemUtils::createIfNotExisting(const boost::filesystem::path& outPath)
{
    using namespace boost::filesystem;
    if(exists(outPath) == false || is_directory(outPath) == false)
    {
        create_directory(outPath);
    }
}

}
