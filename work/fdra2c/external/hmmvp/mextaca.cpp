// mex -g CXXLIBS="\$CXXLIBS -lmwblas -lmwlapack" -I./src -I./include mextaca.cpp src/Compress.cpp src/Mpi.cpp src/Hmat.cpp src/Hd.cpp

#include <mex.h>
#include <math.h>
#include <vector>
#include <list>
#include "include/Exception.hpp"
#include "include/Matrix.hpp"
#include "src/LinAlg.hpp"
#include "src/CodeAnalysis.hpp"
#include "src/Compress_pri.hpp"
#include "matlab/HmmvpMexUtil.hpp"
using namespace std;
using namespace hmmvp;

class MexMa : public hmmvp::MatrixAccessor {
public:
  MexMa(const mxArray* usrfn);
  virtual ~MexMa();
  virtual bool Call(const MatBlock& blk, const vector<uint>* rs,
                    const vector<uint>* cs, double* B);
  
private:
  mxArray* _usrfn;
};

MexMa::MexMa(const mxArray* usrfn)
{
  _usrfn = mxDuplicateArray(usrfn);
}

MexMa::~MexMa()
{
  mxDestroyArray(_usrfn);
}

static mxArray* GetIndices(const MatBlock& blk, const vector<uint>* is, int dim)
{
  if (is) {
    return VectorToMex(*is, (uint) 1);
  } else {
    mxArray* mis;
    if (dim == 1) {
      mis = mxCreateDoubleMatrix(blk.m, 1, mxREAL);
      double* dmis = mxGetPr(mis);
      for (int i = 0; i < blk.m; i++) dmis[i] = (double) (blk.r0 + i + 1);
    } else {
      mis = mxCreateDoubleMatrix(blk.n, 1, mxREAL);
      double* dmis = mxGetPr(mis);
      for (int i = 0; i < blk.n; i++) dmis[i] = (double) (blk.c0 + i + 1);
    }
    return mis;
  }
}

bool MexMa::Call(const MatBlock& blk, const vector<uint>* rs,
                 const vector<uint>* cs, double* B)
{
  int nlhs = 1, nrhs = 3;
  mxArray *prhs[3], *plhs[1];
  // mexCallMATLAB is supposed to let prhs be const, but it isn't.
  prhs[0] = _usrfn;
  prhs[1] = GetIndices(blk, rs, 1);
  prhs[2] = GetIndices(blk, cs, 2);
  plhs[0] = NULL;
  mxArray* merr = mexCallMATLABWithTrap(nlhs, plhs, nrhs, prhs, "feval");
  for (int i = 1; i < 3; i++) mxDestroyArray(prhs[i]);
  if (plhs[0]) {
    memcpy(B, mxGetPr(plhs[0]), mxGetNumberOfElements(plhs[0])*sizeof(double));
    mxDestroyArray(plhs[0]);
  }
  return merr == NULL;
}

void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  int strlen = mxGetNumberOfElements(prhs[0]) + 1;
  vector<char> vfn(strlen);
  mxGetString(prhs[0], &vfn[0], strlen);
  for (int i = 0; i < 1; i++) vfn[i] = tolower(vfn[i]);
  string fn(&vfn[0]);
  prhs++; // Increment the pointer so we can use base-0 offsets in what follows.
  
  if (fn[0] == 'a') {
    if (nlhs != 2 || !(nrhs == 6 || nrhs == 8))
      mexErrMsgTxt("[U V] = aca([r c (base 1) m n], @(r,c), scale, "
                   "use_rel_tol, err, [U, V])");
    double* db = mxGetPr(prhs[0]);
    MatBlock blk((int) db[0] - 1, (int) db[2], (int) db[1] - 1, (int) db[3]);
    MexMa gf(prhs[1]);
    double scale = mxGetScalar(prhs[2]);
    bool use_rel_err = (bool) mxGetScalar(prhs[3]);
    double err = mxGetScalar(prhs[4]);
    Matd U, V;
    if (nrhs == 8) {
      MexToMatrix(prhs[5], U);
      MexToMatrix(prhs[6], V);
    }
    Aca(blk, gf, scale, use_rel_err, err, U, V);
    plhs[0] = MatrixToMex(U);
    plhs[1] = MatrixToMex(V);

  } else if (fn[0] == 'c') {
    if (nlhs != 2 || nrhs != 5)
      mexErrMsgTxt("[U V] = compressqr(U, V, use_rel_err, err)");
    Matd U, V;
    MexToMatrix(prhs[0], U);
    MexToMatrix(prhs[1], V);
    if (U.Size(2) != V.Size(2))
      mexErrMsgTxt("size(U,2) ~= size(V,2)");
    bool use_rel_err = (bool) mxGetScalar(prhs[2]);
    double err = mxGetScalar(prhs[3]);
    CompressQr(U, V, use_rel_err, err);
    plhs[0] = MatrixToMex(U);
    plhs[1] = MatrixToMex(V);

  } else if (fn[0] == 'l') {
    static const unsigned int ndef = 5;
    if (nlhs != 3 || !(nrhs == ndef || nrhs == ndef + 3))
      mexErrMsgTxt("[B U V] = lra([r c (base 1) m n], @(r,c), abs_tol, realp, "
                   "[U, V, old_abs_tol])");
    double* db = mxGetPr(prhs[0]);
    MatBlock blk((int) db[0] - 1, (int) db[2], (int) db[1] - 1, (int) db[3]);
    MexMa gf(prhs[1]);
    double abs_tol = mxGetScalar(prhs[2]);
    unsigned int realp = (unsigned int) mxGetScalar(prhs[3]);
    LraBlock* lb = NewLraBlock(blk, realp, abs_tol);
    LraBlock* ob = NULL;
    if (nrhs > ndef) {
      double old_tol = mxGetScalar(prhs[6]);
      ob = NewLraBlock(blk, 2, old_tol);
      TypedLraBlock<double>* dob = dynamic_cast<TypedLraBlock<double>*>(ob);
      MexToMatrix(prhs[4], dob->U());
      MexToMatrix(prhs[5], dob->V());
    }
    LraOptions opts;
    ApproxByLowRank(opts, lb, gf, ob);
    if (lb->GetPrec() == 1) {
      TypedLraBlock<float>* sob = dynamic_cast<TypedLraBlock<float>*>(lb);
      plhs[0] = MatrixToMex(sob->B());
      plhs[1] = MatrixToMex(sob->U());
      plhs[2] = MatrixToMex(sob->V());
    } else {
      TypedLraBlock<double>* dob = dynamic_cast<TypedLraBlock<double>*>(lb);
      plhs[0] = MatrixToMex(dob->B());
      plhs[1] = MatrixToMex(dob->U());
      plhs[2] = MatrixToMex(dob->V());
    }
    if (ob) delete ob;
  }
}
