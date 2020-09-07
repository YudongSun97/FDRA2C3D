// Test how things work.

// mpic++ -DFDRA_MPI -g -I src try_Mpi.cpp src/Mpi.cpp -o try_Mpi

#include <vector>
using namespace std;
#include <stdio.h>
#include "src/Mpi.hpp"
using namespace fdra;

int main(int argc, char** argv)
{
  mpi::Init(argc, argv);
  int pid = mpi::Pid(), nproc = mpi::GetNproc();
  bool am_root = mpi::AmRoot();

  {
    int foo, bar;
    if (am_root) foo = 10, bar = 20;
    if ((am_root && !mpi::IsTrue(false)) ||
        !mpi::IsTrue())
      printf("pid %d correctly got false.\n", pid);
    if (!mpi::IsTrue(true)) printf("pid %d got false!\n", pid);
    mpi::Bcast(&foo, 1);
    if (mpi::IsTrue(false)) printf("pid %d got true!\n", pid);
    mpi::Bcast(&bar, 1);
    printf("pid %d gets foo = %d bar = %d\n", pid, foo, bar);
  }

#if 1
  int nelem = 10;
  ArraySegmenter as;
  as.Segment(nelem);
  int nmax = nelem / nproc;
  if (nmax * nproc < nelem) nmax++;
  int n = as.GetN();
  vector<int> snd(n);
  for (int i = 0; i < snd.size(); i++) snd[i] = 10*pid + i + 1;
  if (am_root) printf("nelem = %d\n", nelem);
  vector<int> rcv;
  int ret;
  if (nelem >= nproc) { // Gather seems to have a problem with size-0 snd.
    // This shows that a lot of MPI convenience routines require the blocks to
    // be equally sized.
    if (am_root) rcv.resize(nproc * nmax);
    ret = mpi::Gather(&snd[0], snd.size(), &rcv[0], nmax, mpi::Root());
    if (am_root) {
      printf("1 (%d). ", ret);
      for (int i = 0; i < rcv.size(); i++)
        printf("%02d ", rcv[i]);
      printf("\n");
    }
  }
  // However, ArraySegmenter provides an implementation of what we want.
  if (am_root) rcv.resize(nelem);
  ret = as.Gather(&snd[0], &rcv[0], mpi::Root());
  if (am_root) {
    printf("2 (%d). ", ret);
    for (int i = 0; i < rcv.size(); i++)
      printf("%02d ", rcv[i]);
    printf("\n");
  }
  rcv.resize(nelem);
  ret = as.Allgather(&snd[0], &rcv[0]);
  if (pid == nproc - 1) {
    printf("3 (%d). ", ret);
    for (int i = 0; i < rcv.size(); i++)
      printf("%02d ", rcv[i]);
    printf("\n");
  }

  rcv.resize(as.GetN());
  snd.resize(0);
  if (am_root) {
    printf("-> Testing Scatter\n");
    snd.resize(nelem);
    for (int i = 0; i < nelem; i++) snd[i] = i + 1;
  }
  ret = as.Scatter(&snd[0], &rcv[0], mpi::Root());
  printf("pid %d rcv = ", pid);
  for (int i = 0; i < rcv.size(); i++) printf("%02d ", rcv[i]);
  printf("\n");
#endif

  mpi::Finalize();
}
