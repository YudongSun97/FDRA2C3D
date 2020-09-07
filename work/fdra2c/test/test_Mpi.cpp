// mpic++ -DFDRA_MPI -O -I src test_Mpi.cpp src/Mpi.cpp -o test_Mpi
// mpirun -np 8 ./test_Mpi

#include <stdio.h>
#include <vector>
#include <cmath>
#include "test_Parallel.hpp"
#include "src/Mpi.hpp"
using namespace std;
using namespace fdra;

class Sip {
public:
  Sip(int narray, int nelem, int os);

  real Call();

private:
  vector<real> _data;
  vector<real*> _arrs;
  int _narray, _nelem, _os;
};

Sip::Sip(int narray, int nelem, int os)
  : _narray(narray), _nelem(nelem), _os(os)
{
  _data.resize(_narray * _nelem, 0.0);
  _arrs.resize(_narray);
  for (int i = 0; i < _narray; i++)
    _arrs[i] = &_data[0] + i*_nelem;
}

real Sip::Call()
{
  real sum = 0.0;
  for (int ia = 0; ia < _narray; ia++) {
    real* a = _arrs[ia];
    for (int i = 0; i < _nelem; i++) {
      a[i] = DoSomeSillyStuff(a[i], _os + i);
      double sgn = ((_os + i) % 2 == 0) ? 1.0 : -1.0;
      sum += sgn*a[i];
    }
  }
  return sum;
}

int main(int argc, char** argv)
{
  int narray = 20, nelem = 2 << 14, nrepeat = 30;

  mpi::Init(argc, argv);
  bool am_root = mpi::AmRoot();

  ArraySegmenter as;
  as.Segment(nelem);
  int bds[2];
  as.GetIndexBounds(bds);
  Sip sip(narray, bds[1] - bds[0], bds[0]);

  timeval t1, t2;
  mpi::Barrier();
  if (am_root)
    gettimeofday(&t1, 0);
  
  real sum = 0.0;
  for (int i = 0; i < nrepeat; i++) sum += sip.Call();

  real sum_every;
  if (am_root)
    mpi::Reduce<real>(&sum, &sum_every, 1, MPI_SUM);
  else
    mpi::Reduce<real>(&sum, NULL, 1, MPI_SUM);

  mpi::Barrier();
  if (am_root) {
    gettimeofday(&t2, 0);
    double et = difftime(t1, t2);
    printf("MPI:  et = %1.5e  sum = %1.5e\n", et, sum_every);
  }

  mpi::Finalize();
}
