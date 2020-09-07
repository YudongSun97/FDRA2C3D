#ifndef INCLUDE_UTIL_DEFS_PRI
#define INCLUDE_UTIL_DEFS_PRI

#include <vector>
#include "util/include/Defs.hpp"
#include "util/include/Matrix.hpp"

namespace util {

  typedef Matrix<double> Matd;
  typedef Matrix<int> Mati;
  typedef std::vector<double> Vecd;
  typedef std::vector<uint> Vecui;

}

#if __STDC_VERSION__ >= 199901L
# include <stdint.h>
namespace util {
  typedef int64_t int64;
  typedef int32_t int32;
};
#else
namespace util {
  typedef long long int int64;
  typedef int int32;
};
#endif

#endif
