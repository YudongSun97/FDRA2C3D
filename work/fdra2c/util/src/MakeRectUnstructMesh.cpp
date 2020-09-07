#include <stdio.h>
#include <iostream>
#include "util/include/Util.hpp"
#include "util/include/KeyValueFile.hpp"
#include "util/include/RectMeshUnstruct.hpp"
using namespace std;
using namespace util;

namespace {

  void PrintHelp()
  {
  }

  struct Inputs {
    // Minimum and maximum side length of a square.
    double min_len, max_len;
    // Tensor mesh points relative to an arbitrary origin. The domain extends
    // from x(1) to x(end), y(1) to y(end), and so has the dimensions
    //   xlen = x(end) - x(1)
    //   ylen = y(end) - y(1).
    // We use the word "square" to mean the rectangle that has side lengths
    //    xlen / ceil(xlen/max_len)  and  ylen / ceil(ylen/max_len).
    // Hence to get geometric squares, xlen and ylen should be divisible by
    // max_len.
    Matrix<double> x, y;
    // ny x nx array of f(x, y). f should be sampled enough that linear
    // interpolation is sufficient to represent the continuous function. Then
    // the mesh is refined so that max f in a square is equal to the square's
    // side length.
    Matrix<double> f;
    // Output file.
    string save_filename;
  };

  bool ProcessKvf(KeyValueFile* kvf, Inputs& in, string& missing)
  {
    const Matrix<double>* m;
    const string* s;
    if (!kvf->GetMatd("x", m)) { missing = "x"; return false; }
    in.x = *m;
    if (!kvf->GetMatd("y", m)) { missing = "y"; return false; }
    in.y = *m;
    if (!kvf->GetMatd("f", m)) { missing = "f"; return false; }
    in.f = *m;
    if (!kvf->GetDouble("min_len", in.min_len))
      { missing = "min_len"; return false; }
    if (!kvf->GetDouble("max_len", in.max_len))
      { missing = "max_len"; return false; }
    if (!kvf->GetString("save_filename", s))
      { missing = "save_filename"; return false; }
    in.save_filename = *s;
    return true;
  }

  bool TestFileWrite(const string& fn)
  {
    FILE* fid = fopen(fn.c_str(), "w");
    bool ret = fid;
    if (ret) fclose(fid);
    return ret;
  }

  void WriteRects(const rmesh::RectMeshUnstruct& rm, const string& filename)
    throw (FileException)
  {
    FILE* fid = fopen(filename.c_str(), "w");
    if (!fid) throw FileException("Can't read " + filename);
    rm.GetDomain().Serialize(fid);
    rm.GetRectOpts().Serialize(fid);
    const vector<rmesh::Rect>& rs = rm.GetRects();
    long long int n = rs.size();
    write(&n, 1, fid);
    for (size_t i = 0; i < rs.size(); i++) rs[i].Serialize(fid);
    fclose(fid);
  }

}

int main(int argc, char** argv)
{
  if (argc != 2) { PrintHelp(); return -1; }

  Inputs in;
  { KeyValueFile* kvf = NewKeyValueFile();
    if (!kvf->Read(argv[1])) {
      cerr << "Can't read " << argv[1] << endl;
      DeleteKeyValueFile(kvf);
      return -1;
    }
    string missing;
    if (!ProcessKvf(kvf, in, missing)) {
      cerr << "Missing " << missing << endl;
      return -1;
    }
    DeleteKeyValueFile(kvf); }

  if (!TestFileWrite(in.save_filename + ".rect")) {
    cerr << "Can't write " << in.save_filename << endl;
    return -1;
  }

  double xlen = in.x(in.x.Size()) - in.x(1);
  double ylen = in.y(in.y.Size()) - in.y(1);
  rmesh::Rect domain(0.0, 0.0, xlen, ylen);
  rmesh::TensorMeshLinInterpRF rf(in.x, in.y, in.f);
  double check_res = max(xlen, ylen);
  for (int i = 1; i < in.x.Size(); i++)
    check_res = min(check_res, in.x(i+1) - in.x(i));
  for (int i = 1; i < in.y.Size(); i++)
    check_res = min(check_res, in.y(i+1) - in.y(i));
  rmesh::RectOpts ro(in.min_len, in.max_len, check_res);
  rmesh::RectMeshUnstruct rmu(domain, ro, &rf);
  try {
    WriteRects(rmu, in.save_filename + ".rect");
    rmu.Serialize(in.save_filename + ".ser");
  } catch (const FileException& e) {
    cerr << e.GetMsg() << endl;
  }

  return 0;
}
