/* Basic usage in C++.
   
     This is a basic build command. Key parts are linking to the libraries
   LAPACK, BLAS, libgfortran, and of course lib/libhmmvp_s.a and using the flag
   -fopenmp.

     g++ exmvp.cpp -I . -L lib -lhmmvp_s ~/dl/lapack-3.2.1/liblapack.a ~/dl/lapack-3.2.1/libblas.a -lgfortran -fopenmp -o exmvp
 */

#include <stdlib.h>
#include <iostream>
#include "hmmvp/include/Hmat.hpp"
using namespace std;

int main (int argc, char** argv) {
  if (argc != 2) {
    cerr << argv[0] << " [H-matrix filename]" << endl;
    return -1;
  }

  // MVP for one x vector at a time using 4 OpenMP threads.
  hmmvp::Hmat* hm = hmmvp::NewHmat(argv[1], 1, 4);
  size_t m = hm->GetM(), n = hm->GetN();
  cout << "H-matrix is " << m << " by " << m << " with " << hm->GetNnz()
       << " nonzeros for a compression factor of "
       << (double) m * n / hm->GetNnz() << endl;

  // Optionally reorganize memory to encourage continguous memory.
  hm->ReorganizeMemory();

  // Make a random x.
  double* x = new double[n];
  for (size_t i = 0; i < n; i++) x[i] = drand48();

  // Do an MVP.
  double* y = new double[m];  
  hm->Mvp(x, y, 1);

  delete[] x;
  delete[] y;
  hmmvp::DeleteHmat(hm);
}
