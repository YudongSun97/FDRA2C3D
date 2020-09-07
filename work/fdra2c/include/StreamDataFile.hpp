#ifndef INCLUDE_FDRA_STREAMDATAFILE
#define INCLUDE_FDRA_STREAMDATAFILE

#include <string>
#include <vector>
#include "util/include/Exception.hpp"
#include "util/include/Matrix.hpp"
#include "Defs.hpp"

namespace fdra {
  using namespace std;
  using namespace util;

  class StreamDataFile {
  public:
    enum Precision { p_single, p_double };

    // Does not check if nr divides m.Size().
    virtual bool Write(const Matrix<float>& m) = 0;
    virtual bool Write(const Matrix<double>& m) = 0;
    virtual bool Write(const float* m, int n) = 0;
    virtual bool Write(const double* m, int n) = 0;
    virtual void Flush() = 0;

    friend void DeleteStreamDataFile(StreamDataFile* sdf);

  protected:
    StreamDataFile() {}
    StreamDataFile(const StreamDataFile&);
    StreamDataFile& operator=(const StreamDataFile&);
  };

  StreamDataFile* NewStreamDataFileForWriting
    (const string& filename, int nr, StreamDataFile::Precision p);
  void DeleteStreamDataFile(StreamDataFile* sdf);

};

#endif
