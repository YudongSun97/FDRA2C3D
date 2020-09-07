#ifndef INCLUDE_UTIL_DEFS
#define INCLUDE_UTIL_DEFS

#define __HAVE_UINT

#ifdef __HAVE_UINT
#include <sys/types.h>
#endif

namespace util {

#ifndef __HAVE_UINT
  typedef unsigned int uint;
#endif

};

#endif
