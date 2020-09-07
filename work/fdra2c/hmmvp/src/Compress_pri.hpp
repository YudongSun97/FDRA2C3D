#ifndef INCLUDE_HMMVP_COMPRESS_PRI
#define INCLUDE_HMMVP_COMPRESS_PRI

#include <vector>
#include <string>
#include "Exception.hpp"
#include "Mpi.hpp"
#include "Hd_pri.hpp"

namespace hmmvp {
  using namespace util;

  // Row and column global indices for a matrix block. Indexing is 0-based.
  struct MatBlock {
    uint r0, m, c0, n;

    MatBlock();
    MatBlock(uint r0, uint m, uint c0, uint n);
    MatBlock& operator=(const Block& b);
  };

  // The linear algebra routines use this class to access the raw matrix.
  class MatrixAccessor {
  public:
    virtual ~MatrixAccessor() {}
    /* Compute B(rs,cs). Indexing starts at 0. B is preallocated.
         If rs is NULL, then blk.r0 to blk.r0 + blk.m - 1 is requested, and
       similarly for cs.
         Order B with the fastest index going down the columns.
         Return true if all is well; false if there is an error in computing the
       Green's function and you want compression to stop. */
    virtual bool Call(const MatBlock& blk, const vector<uint>* rs,
                      const vector<uint>* cs, double* B) = 0;
  };

  // This is a wrapper to a MatrixAccessor. It accounts for permutations. The
  // user sees this interface. Indexing is 1-based because it's a user-side
  // interface, and I keep all user-side indexing 1-based.
  class GreensFn {
  public:
    virtual ~GreensFn() {}
    virtual bool Call(const vector<uint>& rs, const vector<uint>& cs,
                      double* B) = 0;
  };

  // Essentially a duplicate of ACAr in Bebendorf's AHMED/ACA.h.
  template<typename T>
  void Aca(const MatBlock& blk, MatrixAccessor& ma, double scale,
           const bool use_rel_err, const double err,
           Matrix<T>& U, Matrix<T>& V)
    throw (OutOfMemoryException, UserReqException);

  // Recompression using the QR factorization.
  template<typename T>
  void CompressQr(Matrix<T>& U, Matrix<T>& V,
                  const bool use_rel_err, const double err);

  struct LraOptions {
    enum LraMethod { lra_svd, lra_aca };

    LraMethod method;
    double aca_tol_factor;
    static const double qr_alpha;
    static const uint min_aca_size;

    LraOptions() : method(lra_aca), aca_tol_factor(0.1) {}
  };

  template<typename T>
  struct NumberBlock {
    Matrix<T> B, U, V;
  };

  // An LraBlock packages up everything we need for an H-matrix block. We use it
  // to ferry block data without having to know its type.
  class LraBlock {
  public:
    virtual ~LraBlock() {}

    virtual uint GetPrec() const = 0;
    double GetTol() const { return _tol; }
    const MatBlock& GetBlock() const { return _b; }

    virtual void WriteToFile(FILE* fid) const = 0;
    virtual void Pack(mpi::ByteBufferWriter& bw) const = 0;
    virtual void Unpack(mpi::ByteBufferReader& br) = 0;

  protected:
    LraBlock(const MatBlock& b, double tol) : _b(b), _tol(tol) {}

  private:
    MatBlock _b;
    // This is the abs tol on the block. We are implementing MREM only, and so
    // there is no need to deal with block-wise relative error.
    double _tol;
  };

  template<typename T>
  class TypedLraBlock : public LraBlock {
  public:
    TypedLraBlock(const MatBlock& b, double tol);

    virtual uint GetPrec() const;

    bool HaveB()  const { return _nb.B.Size()  > 0; }
    bool HaveUV() const { return _nb.U.Size()  > 0; }
    uint NcolUV() const { return _nb.U.Size(2); }

    Matrix<T>& B() { return _nb.B; }
    Matrix<T>& U() { return _nb.U; }
    Matrix<T>& V() { return _nb.V; }
    const Matrix<T>& B() const { return _nb.B; }
    const Matrix<T>& U() const { return _nb.U; }
    const Matrix<T>& V() const { return _nb.V; }

    virtual void WriteToFile(FILE* fid) const;
    virtual void Pack(mpi::ByteBufferWriter& bw) const;
    virtual void Unpack(mpi::ByteBufferReader& br);

  private:
    NumberBlock<T> _nb;
  };

  // We use this factory method so that we don't have to worry about type until
  // we get to the details
  LraBlock* NewLraBlock(const MatBlock& b, uint realp, double tol);

  // Makes the high-level low-rank approximation decisions: whether and how to
  // use an old H-matrix block, and whether to use SVD or ACA.
  int ApproxByLowRank(const LraOptions& opts, LraBlock* b, MatrixAccessor& ma,
                      const LraBlock* ob = NULL)
    throw (OutOfMemoryException, UserReqException);

  // Wrapper to the client's Green function to account for H-matrix
  // permutations.
  class GreensFnMa : public MatrixAccessor {
  public:
    GreensFnMa(GreensFn* gf) : _gf(gf) {}

    vector<uint>& pp() { return _pp; }
    vector<uint>& pq() { return _pq; }

    virtual bool Call(const MatBlock& blk, const vector<uint>* rs,
                      const vector<uint>* cs, double* B);

  private:
    GreensFn* _gf;
    vector<uint> _pp, _pq; // 0-based
  };

  class HmatReader;

  // Handles everything. This is exposed to the client as an opaque data type.
  class Compressor {
  public:
    enum TolMethod { tm_Bfro, tm_abs };

    // hd can be NULL for non-root threads.
    Compressor(const Hd* hd, GreensFn* gf, const string& hmat_filename)
      throw (Exception, FileException);
    virtual ~Compressor();

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

    void CompressBlockToFile(FILE* fptr, const MatBlock& b,
                             const LraBlock* ob = NULL);
    const string& GetFilename() const { return _filename; }
    const vector<MatBlock>& GetBlocks() const { return _blocks; }
    uint GetOutputLevel() const { return _output_lvl; }
    HmatReader* GetHmatReader() { return _ohm; }

  private:
    GreensFnMa _ma;

    vector<MatBlock> _blocks;

    string _filename; // Save to this file.

    // Use an old H-matrix to build this one.
    HmatReader* _ohm;
    bool _have_old_hmat;

    // Error tolerance data.
    uint _prec;
    TolMethod _tm;
    double _Bfro;
    double _tol;

    uint _output_lvl;

    vector<int> _ablr_stats; // Record ApproxByLowRank return states.

    static const uint _send_sz;

  private:
    void GetHdData(const Hd* hd);
    void WriteHmatHeader(FILE* fid);
    void CompressSerial();
    void FinalizeTol();
    void ConcatFiles();
    void RecvAndWriteFile(FILE* fid, mpi::ByteBufferReaderMpi& br, uint pid);
    void SendFile();
  };
}

#include "Compress_inl.hpp"

#endif
