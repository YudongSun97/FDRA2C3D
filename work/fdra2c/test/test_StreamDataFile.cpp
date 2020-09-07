// g++ -I include test_StreamDataFile.cpp src/StreamDataFile.cpp

#include <iostream>
#include "StreamDataFile.hpp"
#include "Defs.hpp"
using namespace std;

typedef double real;

int main(int argc, char** argv)
{
  const string test_fn = "test.sdf";
  int nr = 100;
  fdra::StreamDataFile* sdf = fdra::NewStreamDataFileForWriting
    (test_fn, nr, fdra::StreamDataFile::p_single);
  if (!sdf) {
    cerr << "Failed to open " << test_fn << " for writing." << endl;
    return -1;
  }

  fdra::Matrix<real> m(nr);
  real* rm = m.GetPtr();
  for (int i = 0; i < 1000; i++) {
    for (int j = 0; j < nr; j++) rm[j] = i + (real) j / nr;
    if (!sdf->Write(m)) {
      cout << "Failed to write." << endl;
      break;
    }
  }

  fdra::DeleteStreamDataFile(sdf);
}
