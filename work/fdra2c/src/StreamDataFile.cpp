//todo Does fstream have machine-independent I/O?

#include <fstream>
#include "util/include/Exception.hpp"
#include "StreamDataFile.hpp"

namespace fdra {
  using namespace util;

  template<typename T>
  class TypedStreamDataFile : public StreamDataFile {
  public:
    TypedStreamDataFile(const string& fn, int nr)
      throw (FileException);
    virtual ~TypedStreamDataFile();

    virtual bool Write(const Matrix<float>& m);
    virtual bool Write(const Matrix<double>& m);
    virtual bool Write(const float* m, int n);
    virtual bool Write(const double* m, int n);
    virtual void Flush();

  private:
    ofstream _fs;
    int _nr;
    vector<T> _work;

  protected:
    TypedStreamDataFile(const TypedStreamDataFile&);
    TypedStreamDataFile& operator=(const TypedStreamDataFile&);
  };

  template<typename T>
  TypedStreamDataFile<T>::
  TypedStreamDataFile(const string& fn, int nr) throw (FileException)
    : _nr(nr)
  {
    _fs.open(fn.c_str(), ios::out | ios::binary);
    if (!_fs.is_open()) throw FileException();
    // Write header.
    int32 hnr = nr, realp;
    if (sizeof(T) == sizeof(float)) realp = 1;
    else realp = 2;
    if (_fs.write((char*) &hnr, sizeof(int32)).bad() ||
        _fs.write((char*) &realp, sizeof(int32)).bad())
      throw FileException();
  }

  template<typename T>
  TypedStreamDataFile<T>::~TypedStreamDataFile()
  {
    _fs.close();
  }

  template<typename T>
  void TypedStreamDataFile<T>::Flush()
  {
    _fs.flush();
  }

  template<typename T>
  static inline bool WriteSameType(ofstream& fs, const T* m, int n)
  {
    return !fs.write((char*) m, n*sizeof(T)).bad();
  }

  template<typename T1, typename T2>
  static inline bool WriteDiffType
  (ofstream& fs, const T2* m, size_t n, vector<T1>& work)
  {
    if (work.size() < n) work.resize(n);
    for (size_t i = 0; i < n; i++) work[i] = (T1) m[i];
    return !fs.write((char*) &work[0], n*sizeof(T1)).bad();
  }

  template<> bool TypedStreamDataFile<float>::Write(const Matrix<float>& m)
  {    
    return WriteSameType(_fs, m.GetPtr(), m.Size());
  }

  template<> bool TypedStreamDataFile<float>::Write(const float* m, int n)
  {    
    return WriteSameType(_fs, m, n);
  }

  template<> bool TypedStreamDataFile<double>::Write(const Matrix<double>& m)
  {
    return WriteSameType(_fs, m.GetPtr(), m.Size());
  }

  template<> bool TypedStreamDataFile<double>::Write(const double* m, int n)
  {    
    return WriteSameType(_fs, m, n);
  }

  template<> bool TypedStreamDataFile<float>::Write(const Matrix<double>& m)
  {    
    return WriteDiffType<float, double>(_fs, m.GetPtr(), m.Size(), _work);
  }

  template<> bool TypedStreamDataFile<float>::Write(const double* m, int n)
  {    
    return WriteDiffType<float, double>(_fs, m, n, _work);
  }

  template<> bool TypedStreamDataFile<double>::Write(const Matrix<float>& m)
  {    
    return WriteDiffType<double, float>(_fs, m.GetPtr(), m.Size(), _work);
  }

  template<> bool TypedStreamDataFile<double>::Write(const float* m, int n)
  {    
    return WriteDiffType<double, float>(_fs, m, n, _work);
  }

  StreamDataFile* NewStreamDataFileForWriting
  (const string& filename, int nr, StreamDataFile::Precision p)
  {
    if (p == StreamDataFile::p_single)
      return new TypedStreamDataFile<float>(filename, nr);
    else
      return new TypedStreamDataFile<double>(filename, nr);
  }

  void DeleteStreamDataFile(StreamDataFile* sdf)
  {
    if (sdf) delete sdf;
  }

};
