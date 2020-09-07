// g++ -I include test_KeyValueFile.cpp src/KeyValueFile.cpp

#include <iostream>
#include "KeyValueFile.hpp"

int main(int argc, char** argv)
{
  using namespace fdra;

  const string read_fn = "test.ser", write_fn = "testc.ser";
  KeyValueFile* kvf = NewKeyValueFile();
  if (!kvf->Read(read_fn))
    cout << "Failed to read " << read_fn << endl;
  if (!kvf->Write(write_fn))
    cout << "Failed to write " << write_fn << endl;

  const string key = "use_mixed";
  switch (kvf->GetType(key)) {
  case KeyValueFile::vt_none:
    cout << key << ": type couldn't be deduced." << endl;
    break;
  case KeyValueFile::vt_string:
    const string *s;
    kvf->GetString(key, s);
    cout << key << ": " << *s << endl;
    break;
  case KeyValueFile::vt_Matd:
    const Matd *m;
    kvf->GetMatd(key, m);
    cout << key << ": " << *m << endl;
    break;
  }

  DeleteKeyValueFile(kvf);
}
