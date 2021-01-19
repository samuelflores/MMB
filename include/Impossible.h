#ifndef Impossible_H_
#define Impossible_H_

#include <cassert>
#include <cstdlib>

#ifdef NDEBUG
    #define __IMPOSSIBLE__ std::abort()
#else
    #define __IMPOSSIBLE__ assert(false)
#endif // NDEBUG

#endif // Impossible_H_
