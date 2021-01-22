#include <exception>
#include "util.h"

const Log::flush_tag Log::flush = {};
const Log::endl_tag Log::endl = {};

//---------------------------------------------------------------------------
Log::Log(std::ostream &os, int rank) : os(os), rank(rank) {}

//---------------------------------------------------------------------------
Log& Log::operator<<(endl_tag) {
    if (rank == 0) os << std::endl;
    return *this;
}

//---------------------------------------------------------------------------
Log& Log::operator<<(flush_tag) {
    if (rank == 0) os << std::flush;
    return *this;
}
