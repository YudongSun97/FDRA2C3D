// mpic++ -DFDRA_MPI -I include -I src test_Ode.cpp src/Ode.cpp src/KeyValueFile.cpp src/StreamDataFile.cpp src/Mpi.cpp ~/dl/lapack-3.2.1/libblas.a -lgfortran -o test_Ode

#include <stdio.h>
#include <iostream>
#include "include/KeyValueFile.hpp"
#include "src/Ode.hpp"
#include "include/StreamDataFile.hpp"
#include "src/NLA.hpp"
#include "src/Mpi.hpp"
#include "src/WorkArray.hpp"
using namespace fdra;

class AOdeFn : public OdeFunction {
public:
  AOdeFn(const ArraySegmenter* as, const Matrix<real>& A,
         const Matrix<real>& y0, StreamDataFile* sdf);

  virtual int GetN() { return _A.Size(2); }

  virtual void GetIc(real* y0)
  { memcpy(y0, _rwrk.GetPtr(), GetN()*sizeof(real)); }

  virtual void Call
    (double t, double dt, const real* y, real* yd, bool& is_error);

  virtual bool ViewOutput
    (real t, real dt, const real* y, Message m);

private:
  const ArraySegmenter* _as;
  bool _am_root;
  Matrix<real> _A;
  WorkArray<real> _rwrk;
  StreamDataFile* _sdf;
};

AOdeFn::AOdeFn(const ArraySegmenter* as, const Matrix<real>& A,
               const Matrix<real>& y0, StreamDataFile* sdf)
  : _as(as), _A(A), _sdf(sdf)
{
  int m = _A.Size(1);
  _rwrk.Reset(m);
  memcpy(_rwrk.GetPtr(), y0.GetPtr(), m*sizeof(real));
  _am_root = mpi::AmRoot();
}

void AOdeFn::Call(double t, double dt, const real* y, real* yd, bool& is_error)
{
  int m = _A.Size(1), n = _A.Size(2);
  _rwrk.Reset(n);
  _as->Allgather(y, _rwrk.GetPtr());
  gemv('n', m, n, 1.0, _A.GetPtr(), m, _rwrk.GetPtr(), 1, 0.0, yd, 1);
}

bool AOdeFn::ViewOutput(double t, double dt, const real* y, Message m)
{
  int n = _A.Size(2);
  if (_am_root)
    _rwrk.Reset(n);
  _as->Gather(y, _rwrk.GetPtr(), mpi::Root());
  if (_am_root) {
    _sdf->Write(&t, 1);
    _sdf->Write(_rwrk.GetPtr(), n);
  }
  return true;
}

int main(int argc, char** argv)
{
  mpi::Init(argc, argv);

  double ti, tf, reltol, abstol, initial_step;
  Matrix<real> y0, A;

  const string read_fn = "test_ode.ser", stream_fn = "test_ode.dat";
  KeyValueFile* kvf = NewKeyValueFile();
  if (!kvf->Read(read_fn)) {
    cout << "Failed to read " << read_fn << endl;
    return -1;
  }

  const Matrix<double> *py0, *pA;
  if (   !kvf->GetDouble("ti",     ti)
      || !kvf->GetDouble("tf",     tf)
      || !kvf->GetMatd  ("A",      pA)
      || !kvf->GetMatd  ("y0",     py0)
      || !kvf->GetDouble("reltol", reltol)
      || !kvf->GetDouble("abstol", abstol)
      || !kvf->GetDouble("initial_step", initial_step)
      ||  pA->Size(1) != pA->Size(2)
      ||  pA->Size(1) != py0->Size()) {
    cerr << "Input file is not right." << endl;
    DeleteKeyValueFile(kvf);
    return -1;
  }

  ArraySegmenter as;
  as.Segment(py0->Size());
  int bds[2];
  as.GetIndexBounds(bds);
  int nelem = as.GetN();
  y0.Resize(nelem);
  int nA = pA->Size(2);
  int os = as.GetOffset();
  A.Resize(nelem, nA);
  for (int i = 0; i < nelem; i++)
    y0(i+1) = (*py0)(os + i + 1);
  for (int j = 0; j < nA; j++)
    for (int i = 0; i < nelem; i++)
      A(i+1, j+1) = (*pA)(os + i + 1, j+1);
  DeleteKeyValueFile(kvf);

  StreamDataFile::Precision precision = sizeof(real) == sizeof(float) ?
    StreamDataFile::p_single : StreamDataFile::p_double;
  StreamDataFile* sdf = NULL;
  if (mpi::AmRoot())
    sdf = NewStreamDataFileForWriting(stream_fn, 1 + A.Size(2), precision);
  AOdeFn ofn(&as, A, y0, sdf);
  Ode23 ode;
  ode.SetTspan(ti, tf);
  ode.SetRelTol(reltol);
  ode.SetAbsTol(abstol);
  ode.SetInitialStep(initial_step);
  try {
    ode.Run(ofn);
  } catch (const OdeException& oe) {
    cerr << oe.GetMsg() << endl;
  }
  if (mpi::AmRoot())
    DeleteStreamDataFile(sdf);

  mpi::Finalize();
}
