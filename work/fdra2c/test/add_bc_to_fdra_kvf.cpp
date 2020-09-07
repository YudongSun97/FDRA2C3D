/* Read a fdra kvf and a dc3dm bc file. If compatible, stick sum(bc, 2) into the
   kvf and write it. */

#include <stdio.h>
#include "util/include/KeyValueFile.hpp"
#include "util/include/IO.hpp"
using namespace util;

static void add_bc_to_fdra_kvf (const string& fn) throw (Exception) {
  const string* s;
  const Matd* m;
  double d;

  KeyValueFile* kvf = NewKeyValueFile();
  if (!kvf->Read(fn)) {
    DeleteKeyValueFile(kvf);
    throw FileException("Can't read kvf.");
  }

  if (!kvf->GetMatd("a", m)) {
    DeleteKeyValueFile(kvf);
    throw Exception("No field 'a'.");
  }
  const size_t n = m->Size();

  kvf->GetString("hm_filename", s);

  Matd bc(n, 4);
  FILE* fid = fopen((*s + ".bc").c_str(), "rb");
  if (!fid) {
    string msg = string("Can't read bc file: " + *s + string(".bc"));
    DeleteKeyValueFile(kvf);
    throw FileException(msg);
  }
  try {
    read(bc.GetPtr(), 4*n, fid);
  } catch (const FileException& e) {
    fclose(fid);
    DeleteKeyValueFile(kvf);
    throw FileException(string("bc file is not correct: ") + e.GetMsg());
  }
  fclose(fid);

  bool use_edc = false;
  if (kvf->GetDouble("use_edc", d)) use_edc = (bool) d;

  Matd bc_sum(n);
  if (use_edc) {
    for (size_t i = 1; i <= n; i++) bc_sum(i) = bc(i,2) + bc(i,4);
    Matd bc_edc(n);
    for (size_t i = 1; i <= n; i++) bc_edc(i) = bc(i,1) + bc(i,3);
    kvf->AddMatd("hm_bc_edc", bc_edc);
  } else {
    for (size_t i = 1; i <= n; i++) {
      bc_sum(i) = 0;
      for (size_t j = 1; j <= 4; j++) bc_sum(i) += bc(i,j);
    }
  }
  kvf->AddMatd("hm_bc", bc_sum);

  kvf->Write(fn);
  DeleteKeyValueFile(kvf);
}

int main (int argc, char** argv) {
  if (argc != 2) {
    fprintf(stderr, "%s fdra-kvf\n", argv[0]);
    return -1;
  }

  try {
    add_bc_to_fdra_kvf(argv[1]);
  } catch (const Exception& e) {
    fprintf(stderr, "%s: %s\n", argv[0], e.GetMsg().c_str());
    return -1;
  }

  return 0;
}
