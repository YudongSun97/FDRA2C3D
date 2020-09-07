/* I want to try out some ideas that may be useful in implementing an efficient
   MPI-parallelized hm('Compress').

   Here I develop some infrastructure to do a parallel for loop.
 */

// mpic++ -DANALYZE_CODE -DFDRA_MPI -I include -I src src/CodeAnalysis.cpp src/Mpi.cpp test_mpi_parfor.cpp

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include "include/Mpi.hpp"
#include "test_Parallel.hpp"
#include "src/CodeAnalysis.hpp"
using namespace fdra;

static long long int SillySum(int n)
{
  long long int s = 0;
  // Obviously the stupid way of doing this, but I want to kill time.
  for (int i = 1; i <= n; i++) {
    s += i;
    double o = DoSomeSillyStuff((double) i, i);
  }
  return s;
}

class TestPm : public mpi::ParforManager {
public:
  TestPm(int n)
  {
    _sums.resize(n);
    _whogot.resize(n);
    _reqs.resize(mpi::GetNproc() - 1);
  }

  virtual int Isend(int i, int dest)
  {
    return 0;
  }

  virtual int Irecv(int i, int src)
  {
    _whogot[i] = src;
    return mpi::Irecv(&_sums[i], 1, src, _tag, &_reqs[src - 1]);
  }

  void Print()
  {
    for (int i = 0; i < _sums.size(); i += 1000) {
      int d = (int) _sums[i] - i*(i + 1)/2;
      if (d != 0 || i % 1 == 0)
        printf("%5d %3d %10d %10d\n", i, _whogot[i], (int) _sums[i], d);
    }
  }

private:
  static const int _tag = 100;
  vector<int> _whogot;
  vector<long long int> _sums;
  vector<MPI_Request> _reqs;
};

class TestPw : public mpi::ParforWorker {
public:
  virtual int Work(int job_idx, int root)
  {
    Ca::GetTimer()->Tic(1);
    long long int s = SillySum(job_idx);
    int ret = mpi::Send(&s, 1, root, _tag);
    Ca::GetTimer()->Toc(1);
    return ret;
  }

private:
  static const int _tag = 100;
};

int main(int argc, char** argv)
{
  mpi::Init(argc, argv);
  int njobs = 20000;
  TestPm* pm = NULL;
  TestPw* pw = NULL;
  if (mpi::AmRoot()) pm = new TestPm(njobs);
  if (mpi::GetNproc() == 1 || !pm) pw = new TestPw();
  Ca::GetTimer()->Tic(0);
  int ret = Parfor(pm, pw, njobs);
  Ca::GetTimer()->Toc(0);
  if (ret != 0) fprintf(stderr, "Parfor failed.\n");
  if (mpi::AmRoot()) {
    pm->Print();
    delete pm;
  }
  if (pw) delete pw;
  mpi::Barrier();
  caprint("pid %d spent %1.1f%% of its time doing real work\n", mpi::Pid(),
          100.0 * Ca::GetTimer()->TotEt(1) / Ca::GetTimer()->TotEt(0));
  mpi::Barrier();
  if (mpi::AmRoot())
    caprint("Elapsed time = %1.2e\n", Ca::GetTimer()->TotEt(0));
  mpi::Finalize();
}
