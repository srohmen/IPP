#ifndef IPPSTREAM_H
#define IPPSTREAM_H

#include <iostream>
#include <mpi.h>
#include "mpimanager.h"

namespace IPP
{

/// A buffer which reads and writes nothing.
/** This buffer is designed for use in parallel programming.
 *  It is assigned to processors which have potentially no access
 *  to the file system.
 */
class DevNullBuffer : public std::streambuf {
protected:
    virtual int_type overflow(int_type c) {
        return EOF;
    }
    virtual int_type underflow() {
        return EOF;
    }
};

// Parallel Output Streams.

struct Parallel_ostream {
    virtual ~Parallel_ostream() { }
    virtual std::ostream& getOriginalStream() =0;
};

class Parallel_referring_ostream : public Parallel_ostream
{
public:
    Parallel_referring_ostream(std::ostream& original_ostream_)
        : devNullStream(&devNullBuffer),
          original_ostream(original_ostream_)
    {
    }

    virtual std::ostream& getOriginalStream()
    {
        if(MPIManager::getInstance().isMainProc())
        {
            return original_ostream;
        }
        else
        {
            return devNullStream;
        }
    }
private:
    DevNullBuffer devNullBuffer;
    std::ostream  devNullStream;
    std::ostream& original_ostream;
};

template<typename Value>
Parallel_ostream& operator<< (Parallel_ostream& lhs, Value const& rhs) {
    lhs.getOriginalStream() << rhs;
    return lhs;
}

inline Parallel_ostream& operator<< (Parallel_ostream& lhs, std::ostream& (*op)(std::ostream&)) {
    lhs.getOriginalStream() << op;
    return lhs;
}

class ofstream : public Parallel_ostream {
public:
    ofstream();
    explicit ofstream(const char* filename,
                          std::ostream::openmode mode = std::ostream::out | std::ostream::trunc );
    ~ofstream();
    virtual std::ostream& getOriginalStream();

    bool is_open();
    void open(const char* filename, std::ostream::openmode mode = std::ostream::out | std::ostream::trunc);
    void close();
private:
    ofstream(ofstream const& rhs);
    ofstream& operator=(ofstream const& rhs);
private:
    DevNullBuffer devNullBuffer;
    std::ostream  devNullStream;
    std::ofstream *original;
};

extern Parallel_referring_ostream pcout;
extern Parallel_referring_ostream pcerr;
extern Parallel_referring_ostream pclog;

}

#endif // IPPSTREAM_H
