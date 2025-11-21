#include "ippstream.h"

namespace IPP
{

Parallel_referring_ostream pcout(std::cout);
Parallel_referring_ostream pcerr(std::cerr);
Parallel_referring_ostream pclog(std::clog);

}
