/* This experiment tests the effect of memory layout on OpenMP performance. I
   identify two basic layout patterns that can be considered careful.

   The first pattern allocates each array in full, regardless of number of
   threads. Then each thread views a piece of each array. If all the memory is
   allocated up front in one long data buffer, then the memory layout is a
   sequence of small blocks, each owned by a thread. This pattern is careful in
   that allocation is done once up front, read and write access to blocks is
   carefully considered, and read and write patterns on a piece of memory does
   not change over time. This pattern relies on OpenMP and caching behavior to
   keep only the pieces needed by each thread in the corresponding core.

   The second pattern is motivated by the fear that the first pattern asks too
   much of OpenMP. It allocates one large block for each thread
   independently. Then each of these blocks is broken up into pieces of each
   array. This is basically what one does when using MPI (though with MPI it is
   at the process, rather than thread, level). This pattern is obviously
   careful.

   Conclusions: I ran this program on both rainier (Intel i7) and redoubt (AMD
   Opteron) for a range of values. I find that the run time of pattern 2 is a
   few percent (~3%) lower than that for pattern 1.

   Next I want to compare OpenMP and MPI implementations. There are a fair
   number of papers comparing the two on shared-memory machines, but either (a)
   the paper is too old to be useful or (b) the paper has some serious flaws in
   the setup and so I don't believe the results. mpic++ and g++ use the same GNU
   compiler. The file is test_Mpi.cpp.

   Conclusions: I find that the MPI version is slower than OpenMP pattern 2 by a
   few percent (again ~3%, sometimes more like 4 or 5, sometimes less) on
   smaller problems. Often it is slower than pattern 1, too. On large-memory
   problems, it seems to do better than pattern 1 and almost as well, if not as
   well, as pattern 2. All of this makes sense.

   There are a few takeaways:

   0. In pattern 1, multithreading comes down to having loop bounds. In pattern
      2, every object must allocate according to number of threads. Hence
      pattern 2 is quite invasive.

   1. As long as one is careful about memory access and do resonable memory
      layout, one is within a few percent of optimal. A few percent is not worth
      making things ridiculously inconvenient, so go with pattern 1.

   2. Same deal with MPI vs OpenMP on the same machine. If you want the few
      percent, don't go with OpenMP pattern 2; rather, just do MPI, since the
      two have the same feel. I personally don't care about anything <10%, so
      OpenMP pattern 1 is the natural choice if I'm sure I won't be running on a
      distributed system.

   3. From the opposite perspective, there's not a meaningful performance hit
      when using MPI on a shared-memory machine rather than OpenMP. So it's
      quite reasonable never to use OpenMP for standalone programs if you want
      to be sure you can always run on a distributed machine if necessary.
 */

// g++ -O -DFDRA_OMP -I src test_OpenMP.cpp -fopenmp -o test_OpenMP

#include <stdio.h>
#include <vector>
#include <cmath>
#include "test_Parallel.hpp"
using namespace std;


#ifdef FDRA_OMP
# include <omp.h>
#else
# define omp_set_num_threads(nthreads);
# define omp_get_max_threads() 1
# define omp_get_thread_num() 0
#endif

class ArraySegmenter {
public:
  void Segment(int n, int nthreads);

  int GetNumberThreads() const { return _bds.size() - 1; }

  // The caller is responsible for indices bounds[0]:bounds[1]-1. Hence a
  // typical for loop appears as
  //     as.GetIndexBounds(bds, 2);
  //     for (int i = bds[0]; i < bds[1]; i++)
  void GetIndexBounds(int bounds[2], int factor = 1) const
  {
    int tid = omp_get_thread_num();
    bounds[0] = _bds[tid] * factor;
    bounds[1] = _bds[tid + 1] * factor;
  }

  int GetOffset(int tid) const
  {
    return _bds[tid];
  }

private:
  vector<int> _bds;
};

void ArraySegmenter::Segment(int n, int nthreads)
{
  omp_set_num_threads(nthreads);
  nthreads = omp_get_max_threads();

  _bds.resize(nthreads + 1);
  for (int i = 0, b = 0, k = n / nthreads, extra = n - k*nthreads;
       i < nthreads; i++) {
    _bds[i] = b;
    b += k;
    if (extra > 0) {
      b++;
      extra--;
    }
  }
  _bds[nthreads] = n;
}

// SIP = silly intensive process
class Sip {
public:
  Sip(int narray, int nelem, int nthreads)
    : _narray(narray), _nelem(nelem), _nthreads(nthreads) {}
  virtual ~Sip() {}

  virtual real Call() = 0;

protected:
  int _narray, _nelem, _nthreads;
};

// Each of the k arrays is allocated in one malloc. Hence the memory used in one
// thread is checkered throughout the total allocated memory.
class Sip1 : public Sip {
public:
  Sip1(int narray, int nelem, int nthreads);
  virtual real Call();

private:
  vector<real> _data, _for_sum;
  vector<real*> _arrs;
  ArraySegmenter _as;
};

Sip1::Sip1(int narray, int nelem, int nthreads)
  : Sip(narray, nelem, nthreads)
{
  _as.Segment(nelem, nthreads);
  _nthreads = _as.GetNumberThreads();
  printf("nthreads = %d\n", _nthreads);
  _data.resize(narray * nelem, 0.0);
  _for_sum.resize(_nthreads);
  _arrs.resize(narray);
  for (int i = 0; i < narray; i++)
    _arrs[i] = &_data[0] + i*nelem;
}

real Sip1::Call()
{
#pragma omp parallel
  {
    real sum = 0.0;
    int bds[2];
    _as.GetIndexBounds(bds);
    for (int ia = 0; ia < _narray; ia++) {
      real* a = _arrs[ia];
      for (int i = bds[0]; i < bds[1]; i++) {
        a[i] = DoSomeSillyStuff(a[i], i);
        double sgn = (i % 2 == 0) ? 1.0 : -1.0;
        sum += sgn*a[i];
      }
    }
    int tid = omp_get_thread_num();
    _for_sum[tid] = sum;
  }

  real sum = 0.0;
  for (int i = 0; i < _nthreads; i++) sum += _for_sum[i];
  return sum;
}

// Each thread allocates its own part of the arrays. The question I have is how
// much faster this is, if at all.
class Sip2 : public Sip {
public:
  Sip2(int narray, int nelem, int nthreads);
  virtual real Call();

private:
  vector< vector<real> > _data;
  vector< vector<real*> > _arrs;
  vector<real> _for_sum;
  vector<int> _tnelems;
  ArraySegmenter _as;
};

Sip2::Sip2(int narray, int nelem, int nthreads)
  : Sip(narray, nelem, nthreads)
{
  _as.Segment(nelem, nthreads);
  _nthreads = _as.GetNumberThreads();
  printf("nthreads = %d\n", _nthreads);
  _data.resize(_nthreads);
  _arrs.resize(_nthreads);
  _for_sum.resize(_nthreads);
  _tnelems.resize(_nthreads);
  for (int tid = 0; tid < _nthreads; tid++) {
    _arrs[tid].resize(_narray);
    int bds[2];
    _as.GetIndexBounds(bds);
    int tnelem = bds[1] - bds[0];
    _tnelems[tid] = tnelem;
    _data[tid].resize(_narray * tnelem, 0.0);
    for (int ia = 0; ia < _narray; ia++)
      _arrs[tid][ia] = &_data[tid][0] + ia*tnelem;
  }
}

real Sip2::Call()
{
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    real sum = 0.0;
    int os = _as.GetOffset(tid);
    for (int ia = 0; ia < _narray; ia++) {
      real* a = _arrs[tid][ia];
      for (int i = 0, n = _tnelems[tid]; i < n; i++) {
        a[i] = DoSomeSillyStuff(a[i], os + i);
        double sgn = ((os + i) % 2 == 0) ? 1.0 : -1.0;
        sum += sgn*a[i];
      }
    }
    _for_sum[tid] = sum;
  }

  real sum = 0.0;
  for (int i = 0; i < _nthreads; i++)
    sum += _for_sum[i];
  return sum;
}

int main(int argc, char** argv)
{
  int narray = 20, nelem = 2 << 24, nthreads = 16, nrepeat = 30;
  real sum1 = 0.0, sum2 = 0.0;
  timeval t1, t2;
  double et1, et2;

  for (int k = 0; k < 2; k++) {
    {
      Sip1 sip1(narray, nelem, nthreads);
      gettimeofday(&t1, 0);
      for (int i = 0; i < nrepeat; i++) sum1 += sip1.Call();
      gettimeofday(&t2, 0);
      et1 = difftime(t1, t2);
    }
    printf("1.  et = %1.5e  sum = %1.5e\n", et1, sum1);

    {
      Sip2 sip2(narray, nelem, nthreads);
      gettimeofday(&t1, 0);
      for (int i = 0; i < nrepeat; i++) sum2 += sip2.Call();
      gettimeofday(&t2, 0);
      et2 = difftime(t1, t2);
    }
    printf("2.  et = %1.5e  sum = %1.5e\n", et2, sum2);
  }
  printf("sum1 - sum2 = %e\n", sum1 - sum2);
}
