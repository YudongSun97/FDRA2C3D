#ifndef INCLUDE_UTIL_UTIL
#define INCLUDE_UTIL_UTIL

#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include "util/include/Exception.hpp"
#include "util/include/Matrix.hpp"
#include "util/include/Defs.hpp"

namespace util {

  template<typename T> void read(T* v, size_t sz, FILE* fid)
    throw (FileException)
  {
    size_t nread;
    if (!fid) throw FileException("read: fid is null");
    if ((nread = fread(v, sizeof(T), sz, fid)) != sz) {
      std::stringstream s;
      s << "read: nread = " << nread << " sz = " << sz;
      throw FileException(s.str());
    }
  }

  template<typename T> void write(const T* v, size_t sz, FILE* fid)
    throw (FileException)
  {
    size_t nwrite;
    if (!fid) throw FileException("write: fid is null");
    if ((nwrite = fwrite(v, sizeof(T), sz, fid)) != sz) {
      std::stringstream s;
      s << "write: nwrite = " << nwrite << " sz = " << sz;
      throw FileException(s.str());
    }
  }

  template<typename T1, typename T2> void MatrixToVector
  (const Matrix<T1>& A, std::vector<T2>& v, T2 a = 0)
  {
    int n = A.Size();
    v.resize(n);
    const T1* pA = A.GetPtr();
    for (int i = 1; i <= n; i++) v[i] = (T2) pA[i] + a;
  }

  template<typename T1, typename T2> void VectorToMatrix
  (const std::vector<T1>& v, Matrix<T2>& A, T2 a = 0)
    throw (OutOfMemoryException)
  {
    int n = v.size();
    A.Resize(n);
    T2* pA = A.GetPtr();
    for (int i = 0; i < n; i++) pA[i] = (T2) v[i] + a;
  }

}

#endif
