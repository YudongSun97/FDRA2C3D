#ifndef INCLUDE_FDRA_DEFS
#define INCLUDE_FDRA_DEFS

#include "util/include/Defs.hpp"

namespace fdra {
  typedef double real;
  typedef float tau_real;
};

#if __STDC_VERSION__ >= 199901L
# include <stdint.h>
namespace fdra {
  typedef int64_t int64;
  typedef int32_t int32;
};
#else
namespace fdra {
  typedef long long int int64;
  typedef int int32;
};
#endif

#endif
