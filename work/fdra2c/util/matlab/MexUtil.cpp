#include <algorithm>
#include "mex.h"
#include "MexUtil.hpp"
using namespace std;

namespace util {

  template<> mxArray* MatrixToMex<double>(const Matrix<double>& m)
    throw (OutOfMemoryException)
  {
    int nr = m.Size(1), nc = m.Size(2);
    mxArray* ma = mxCreateDoubleMatrix(nr, nc, mxREAL);
    if (!ma) throw OutOfMemoryException("MatrixToMex");
    memcpy(mxGetPr(ma), m.GetPtr(), nr*nc*sizeof(double));
    return ma;
  }

  template<> void MexToMatrix<double>(const mxArray* ma, Matrix<double>& m)
    throw (OutOfMemoryException)
  {
    int nr = mxGetM(ma), nc = mxGetN(ma);
    m.Resize(nr,nc);
    memcpy(m.GetPtr(), mxGetPr(ma), nr*nc*sizeof(double));
  }

};
