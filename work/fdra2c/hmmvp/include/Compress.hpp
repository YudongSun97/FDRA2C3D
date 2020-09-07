#ifndef INCLUDE_HMMVP_COMPRESS
#define INCLUDE_HMMVP_COMPRESS

#include <vector>
#include <string>
#include "util/include/Exception.hpp"
#include "hmmvp/include/Hd.hpp"

namespace hmmvp {
  using namespace util;

  class GreensFn {
  public:
    virtual ~GreensFn() {}
    /* Compute B(rs,cs). Indexing starts at 1. B is preallocated.
         Return true if all is well; false if there is an error in computing the
       Green's function and you want compression to stop*/
    virtual bool Call(const vector<uint>& rs, const vector<uint>& cs,
                      double* B) = 0;
  };

  class Compressor {
  private:
    virtual ~Compressor();

  public:
    enum TolMethod { tm_Bfro, tm_abs };

    void SetTolMethod(TolMethod tm);
    void SetTol(double tol) throw (Exception);
    void SetBfroEstimate(double Bfro) throw (Exception);
    // Use another H-matrix file to speed up compressing this one. Returns false
    // if file is incompatible.
    void UseHmatFile(const string& hmat_filename)
      throw (Exception, FileException);
    bool HaveOldHmat() const;
    double GetOldHmatBfro() const throw (Exception);
    void SetOutputLevel(uint lev);

    double EstimateBfro()
      throw (OutOfMemoryException, UserReqException);

    void Compress()
      throw (OutOfMemoryException, UserReqException, FileException);

  private:
    Compressor();
    Compressor(const Compressor&);
    Compressor& operator=(const Compressor&);
  };

  Compressor* NewCompressor
   (const Hd* hd, GreensFn* gf, const string& hmat_filename)
    throw (Exception, FileException);

  void DeleteCompressor(Compressor* c);

}

#endif
