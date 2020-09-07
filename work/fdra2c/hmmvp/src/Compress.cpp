/* todo
   - Allow for UserReqException. Right now the code will just abort.
 */

#include <stdio.h>
#include <sstream>
#include "Mpi.hpp"
#include "Hmat.hpp"
#include "Util.hpp"
#include "HmatIo.hpp"
#include "Compress_pri.hpp"
#include "CodeAnalysis.hpp"

namespace hmmvp {
  using namespace util;

  // ---------------------------------------------------------------------------

  const double LraOptions::qr_alpha = 0.5;
  const uint LraOptions::min_aca_size = 251;

  MatBlock::MatBlock() : r0(0), m(0), c0(0), n(0) {}

  MatBlock::MatBlock(uint r0_, uint m_, uint c0_, uint n_)
      : r0(r0_), m(m_), c0(c0_), n(n_)
  {
    assert(r0 >= 0 && c0 >= 0 && m >= 1 && n >= 1);
  }

  MatBlock& MatBlock::operator=(const Block& b)
  {
    r0 = b.r0 - 1; m = b.m; c0 = b.c0 - 1; n = b.n;
    return *this;
  }
  
  template<typename T>
  void TypedLraBlock<T>::WriteToFile(FILE* fid) const
  {
    const MatBlock& mb = GetBlock();
    Blint n;
    n = mb.r0; write(&n, 1, fid);
    n = mb.m;  write(&n, 1, fid);
    n = mb.c0; write(&n, 1, fid);
    n = mb.n;  write(&n, 1, fid);
    if (HaveB()) {
      char s = 'B'; write(&s, 1, fid);
      write(B().GetPtr(), B().Size(), fid);
    } else {
      char s = 'U'; write(&s, 1, fid);
      n = U().Size(2); write(&n, 1, fid);
      write(U().GetPtr(), U().Size(), fid);
      Matrix<T> Vt;
      Transpose(V(), Vt);
      write(Vt.GetPtr(), Vt.Size(), fid);
    }
  }

  template<typename T>
  void TypedLraBlock<T>::Pack(mpi::ByteBufferWriter& bw) const
  {
    if (HaveB()) {
      bw.WriteScalar((char) 'B');
      PackMatrix(bw, B());
    } else {
      bw.WriteScalar((char) 'U');
      PackMatrix(bw, U());
      PackMatrix(bw, V());
    }
  }

  template<typename T>
  void TypedLraBlock<T>::Unpack(mpi::ByteBufferReader& br)
  {
    char code = br.ViewScalarAndAdvance<char>();
    if (code == 'B') {
      UnpackMatrix(br, B());
    } else {
      UnpackMatrix(br, U());
      UnpackMatrix(br, V());
    }
  }

  LraBlock* NewLraBlock(const MatBlock& b, uint realp, double tol)
  {
    if (realp == 1) return new TypedLraBlock<float>(b, tol);
    else return new TypedLraBlock<double>(b, tol);
  }

  // ---------------------------------------------------------------------------
  // Low-rank approximation. No need to carry over all of the options in
  // hm('Compress'). In particular, I'm going to implement only MREM here.

  template<typename T>
  class FullBMa : public MatrixAccessor {
  public:
    FullBMa(const Matrix<T>* B) : _B(B) {}
    virtual bool Call(const MatBlock& blk, const vector<uint>* rs,
                      const vector<uint>* cs, double* B);
  private:
    const Matrix<T>* _B;
  };

  template<typename T>
  bool FullBMa<T>::Call(const MatBlock& blk, const vector<uint>* rs,
                        const vector<uint>* cs, double* B)
  {
    if (!rs && !cs) {
      Matrix<double> wB;
      wB.SetPtr(blk.m, blk.n, B);
      wB = *_B;
      return true;
    }
    if (!rs) {
      for (uint ic = 0, nc = cs->size(), k = 0; ic < nc; ic++) {
        const T* p_B = _B->GetPtr() + (*cs)[ic] * blk.m;
        for (uint ir = 0; ir < blk.m; ir++, k++) B[k] = (double) p_B[ir];
      }
      return true;
    }
    if (!cs) {
      const T* p_B = _B->GetPtr();
      uint nr = rs->size();
      for (uint ic = 0, k = 0; ic < blk.n; ic++) {
        for (uint ir = 0; ir < nr; ir++, k++) B[k] = (double) p_B[(*rs)[ir]];
        p_B += blk.m;
      }
      return true;
    }
    // Should never have both rs and cs.
    assert(false);
    return false;
  }

  template<typename T>
  inline uint SvdChooseRank(const Matrix<T>& s, double tol)
  {
    uint n = s.Size();
    double tol2 = tol*tol;
    T sum = 0.0;
    uint rank;
    for (rank = n; rank >= 1; rank--) {
      sum += s(rank)*s(rank);
      if (sum >= tol2) break;
    }
    if (rank == 0) rank = 1;
    return rank;
  }

  template<typename T, typename OT>
  inline int ApproxByLowRank(const LraOptions& opts, TypedLraBlock<T>& b,
                             MatrixAccessor& ma, const TypedLraBlock<OT>* ob)
    throw (OutOfMemoryException, UserReqException)
  {
    const MatBlock& blk = b.GetBlock();

    int ret = 0;
    if (ob) {
      double old_tol = ob->GetTol();
      double new_tol = b.GetTol();
      // ||B_single - B_double||_F ~= eps(single) sqrt(m n)
      bool trust_ob =
        new_tol >= numeric_limits<OT>::epsilon() * sqrt(blk.m * blk.n);
      if (trust_ob) {
        if (new_tol < old_tol && ob->HaveB()) {
          // We need a full B, and we can use the old one.
          b.B() = ob->B();
          return 1;
        } else if (new_tol > old_tol && !ob->HaveB()) {
          // We can use the old and U and V.
          b.U() = ob->U();
          b.V() = ob->V();
          if (b.NcolUV() > 1) {
            // We have ||C - A|| <= old_tol and want ||C - B|| <= new_tol. Hence
            //     ||C - B|| = ||(C - A) + (A - B)|| <= ||C - A|| + ||A - B||
            //       <= old_tol + tol <= new_tol
            // and so tol = new_tol - old_tol.
            double tol = new_tol - old_tol;
            b.U() = ob->U();
            b.V() = ob->V();
            CompressQr(b.U(), b.V(), false, tol);
          }
          return 2;
        }
      } else {
        ob = NULL;
        ret = 100;
      }
    }

    if (blk.m == 1 || blk.n == 1 || (blk.m <= 2 && blk.n <= 2)) {
      // To small, so just grab B.
      b.U().Resize(0);
      b.V().Resize(0);
      if (ob && ob->HaveB()) {
        ret = 10;
        b.B() = ob->B();
        return 1;
      } else {
        Matrix<double> B(blk.m, blk.n);
        ma.Call(blk, NULL, NULL, B.GetPtr());
        b.B() = B;
        return ret + 0;
      }
    }

    LraOptions::LraMethod method = opts.method;
    if (std::max(blk.m, blk.n) < opts.min_aca_size)
      method = LraOptions::lra_svd;

    switch (method) {
    case LraOptions::lra_svd: {
      Matrix<T> s, Vt;
      Matrix<double> B;
      if (ob && ob->HaveB()) {
        ret = 10;
        Svd(ob->B(), b.U(), s, Vt);
      } else {
        B.Resize(blk.m, blk.n);
        ma.Call(blk, NULL, NULL, B.GetPtr());
        Svd(B, b.U(), s, Vt);
      }
      uint rank = SvdChooseRank(s, b.GetTol());
      if ((blk.m + blk.n) * rank < blk.m*blk.n) {
        T* pU = b.U().Reshape(blk.m, rank).GetPtr();
        for (uint j = 0, k = 0; j < rank; j++) {
          T sj = s(j+1);
          for (uint i = 0; i < blk.m; i++, k++) pU[k] = pU[k] * sj;
        }
        T* pV = b.V().Resize(blk.n, rank).GetPtr();
        for (uint j = 0, k = 0; j < rank; j++)
          for (uint i = 0; i < blk.n; i++, k++) pV[k] = Vt(j+1, i+1);
      } else {
        // Better just to keep the whole block.
        b.U().Resize(0);
        b.V().Resize(0);
        if (ob && ob->HaveB()) b.B() = ob->B();
        else b.B() = B;
      }
      return ret + 3;
    }

    case LraOptions::lra_aca: {
      MatrixAccessor* pma = &ma;
      bool ma_alloced = false;
      // Make use of old data, if available.
      if (ob) {
        if (ob->HaveB()) {
          ret = 20;
          ma_alloced = true;
          pma = new FullBMa<OT>(&ob->B());
        } else {
          ret = 30;
          b.U() = ob->U();
          b.V() = ob->V();
        }
      }
      // alpha splits up the error between Aca and CompressQr.
      double alpha = 0.5;
      // In addition, we set this first error to ~1/10 of the requested error
      // because ACA sometimes is not accurate enough by ~<= 10x. CompressQr
      // generally makes up for this. (aca_factor is nominally 1/10, but can
      // have another value.)
      double tol = alpha * opts.aca_tol_factor * b.GetTol();
      // -1.0 just forces the first irow to be accepted (unless ob).
      double scale = -1.0;
      Aca(blk, *pma, scale, false, tol, b.U(), b.V());
      if ((blk.m + blk.n) * b.NcolUV() >= blk.m*blk.n) {
        // Better just to keep the whole block.
        b.U().Resize(0);
        b.V().Resize(0);
        if (ob && ob->HaveB())
          b.B() = ob->B();
        else {
          Matrix<double> B(blk.m, blk.n);
          pma->Call(blk, NULL, NULL, B.GetPtr());
          b.B() = B;
        }
      } else {
        double tol = (1.0 - alpha) * b.GetTol();
        CompressQr(b.U(), b.V(), false, tol);
      }
      if (ma_alloced) delete pma;
      return ret + 4;
    }}

    // Should never get here.
    assert(false);
    return -1;
  }

  //todo I think I should use a visitor pattern here to implement double
  // dispatch on b and ob. For now, I'm going to commit the crime of using RTTI
  // in the following two functions. I could implement single dispatch on b just
  // by making a method wrapper to ApproxByLowRank for b, but as long as I'm
  // commiting one crime, I'll just go ahead and commit two.
  template<typename T>
  inline int ApproxByLowRank(const LraOptions& opts, LraBlock* b,
                             MatrixAccessor& ma, const TypedLraBlock<T>* ob)
  {
    if (b->GetPrec() == 1)
      return ApproxByLowRank(opts, *(static_cast<TypedLraBlock<float>*>(b)),
                             ma, ob);
    else
      return ApproxByLowRank(opts, *(static_cast<TypedLraBlock<double>*>(b)),
                             ma, ob);
  }

  int ApproxByLowRank(const LraOptions& opts, LraBlock* b, MatrixAccessor& ma,
                      const LraBlock* ob)
    throw (OutOfMemoryException, UserReqException)
  {
    if (ob) {
      if (ob->GetPrec() == 1)
        return ApproxByLowRank(opts, b, ma,
                               static_cast<const TypedLraBlock<float>*>(ob));
      else
        return ApproxByLowRank(opts, b, ma,
                               static_cast<const TypedLraBlock<double>*>(ob));
    } else {
      return ApproxByLowRank(opts, b, ma, (TypedLraBlock<float>*) NULL);
    }
  }

  // ---------------------------------------------------------------------------
  // Compressor.

  // Print to stdout "msg: 1 2 3 ... 99\n" to indicate progress.
  class ProgressBar {
  public:
    ProgressBar(const string& msg, uint total, uint output_lvl = 1);
    void Incr(uint amount = 1);

  private:
    uint _i, _total, _progress, _output_lvl;
  };

  ProgressBar::ProgressBar(const string& msg, uint total, uint output_lvl)
    : _i(0), _total(total), _progress(0), _output_lvl(output_lvl)
  {
    if (_output_lvl > 0) {
      fprintf(stdout, "%s: ", msg.c_str());
      fflush(stdout);
    }
  }

  void ProgressBar::Incr(uint amount)
  {
    if (_output_lvl > 0) {
      _i += amount;
      double p = 100.0 * (double) _i / _total;
      if ((uint) p > _progress) {
        _progress = (uint) p;
        fprintf(stdout, "%d ", _progress);
        fflush(stdout);
        if (_progress == 100) fprintf(stdout, "\n");
      }
    }
  }

  // Read an H-matrix sequentially one block at a time, return an LraBlock for
  // each block as requested.
  class HmatReader {
  private:
    string _filename;
    FILE* _fid;
    double _tol;
    Blint _tol_denom;
    Blint _realp, _nb;
    bool _is_compatible;
    vector<FileBlint> *_p, *_q;

  public:
    HmatReader(const string& filename) throw (FileException);
    ~HmatReader();

    const string& GetFilename() const { return _filename; }
    double GetTol() const { return _tol; }

    // This checks everything.
    bool IsCompatible(const vector<uint>& p, const vector<uint>& q,
                      const vector<const MatBlock*>& sorted_blocks);
    LraBlock* ReadNextBlock() const;
  };

  HmatReader::HmatReader(const string& filename) throw (FileException)
    : _filename(filename), _is_compatible(false)
  {
    _fid = fopen(filename.c_str(), "r");
    if (!_fid) throw FileException("HmatReader");
    Blint m, n;
    _p = new vector<FileBlint>();
    _q = new vector<FileBlint>();
    ReadHmatHeader(_fid, m, n, _realp, _nb, _tol, *_p, *_q);
    _tol_denom = m*n;
  }

  HmatReader::~HmatReader()
  {
    fclose(_fid);
    if (_p) delete _p;
    if (_q) delete _q;
  }

  // Comparison of MatBlock* based on r0, then c0.
  bool MatBlockLessThanRC(const MatBlock* b1, const MatBlock* b2)
  {
    if (b1->r0 < b2->r0) return true;
    if (b1->r0 > b2->r0) return false;
    return b1->c0 < b2->c0;
  }

  bool HmatReader::IsCompatible(const vector<uint>& p, const vector<uint>& q,
                                const vector<const MatBlock*>& sbs)
  {
    if (_p) {
      // First check the header data.
      for (;;) { // breakable
        if (!(_is_compatible = p.size() == _p->size() &&
              q.size() == _q->size() && sbs.size() == (size_t) _nb)) break;
        for (size_t i = 0; i < p.size(); i++)
          if (!(_is_compatible = p[i] == (*_p)[i] + 1)) break;
        for (size_t i = 0; i < q.size(); i++)
          if (!(_is_compatible = q[i] == (*_q)[i] + 1)) break;
        break;
      }
      // We don't need these permutations anymore.
      delete _p; _p = NULL;
      delete _q; _q = NULL;
      if (!_is_compatible) return _is_compatible;

      // Checks so far, so now look at the blocks.
      for (uint i = 0; i < sbs.size(); i++) {
        Blint r0, c0, m, n, rank;
        ReadHmatBlockInfo(_fid, _realp, r0, c0, m, n, rank);
        MatBlock ob(r0, m, c0, n);
        vector<const MatBlock*>::const_iterator it =
          lower_bound(sbs.begin(), sbs.end(), &ob, MatBlockLessThanRC);
        const MatBlock* pb = *it;
        if (!(_is_compatible = pb->r0 == r0 && pb->c0 == c0 &&
              pb->m == m && pb->n == n)) break;
      }

      // Rewind to just after the header.
      fclose(_fid);
      _fid = fopen(_filename.c_str(), "r");
      Blint m, n;
      ReadHmatHeader(_fid, m, n, _realp, _nb, _tol);
    }
    return _is_compatible;
  }

  template<typename T>
  static LraBlock* ReadNextBlock(FILE* fid, double tol, Blint tol_denom)
  {
    Blint r0, c0, m, n, rank;
    T *B, *U, *Vt;
    ReadHmatBlock<T>(fid, r0, c0, m, n, rank, B, U, Vt);
    MatBlock b(r0, m, c0, n);
    tol *= sqrt((double) (m * n) / tol_denom);
    TypedLraBlock<T>* lb = new TypedLraBlock<T>(b, tol);
    if (B) {
      lb->B().Resize(m, n);
      memcpy(lb->B().GetPtr(), B, m*n*sizeof(T));
      delete[] B;
    } else {
      lb->U().Resize(m, rank);
      memcpy(lb->U().GetPtr(), U, m*rank*sizeof(T));
      delete[] U;
      lb->V().Resize(n, rank);
      Matrix<T> mVt(rank, n, Vt);
      Transpose(mVt, lb->V());
      delete[] Vt;
    }
    return lb;
  }

  LraBlock* HmatReader::ReadNextBlock() const
  {
    if (_realp == 1)
      return hmmvp::ReadNextBlock<float>(_fid, _tol, _tol_denom);
    else
      return hmmvp::ReadNextBlock<double>(_fid, _tol, _tol_denom);
  }

  // Pack and unpack basic objects into/from ByteBuffers.
  template<typename T>
  inline void PackVec(mpi::ByteBufferWriter& bw, const vector<T>& v)
  {
    bw.WriteScalar((size_t) v.size());
    bw.Write(&v[0], v.size());
  }

  template<typename T>
  inline void UnpackVec(mpi::ByteBufferReader& br, vector<T>& v)
  {
    v.resize(br.ViewScalarAndAdvance<size_t>());
    br.Read(&v[0], v.size());
  }

  template<typename T>
  inline void PackMatrix(mpi::ByteBufferWriter& bw, const Matrix<T>& A)
  {
    bw.WriteScalar((size_t) A.Size(1));
    bw.WriteScalar((size_t) A.Size(2));
    bw.Write(A.GetPtr(), A.Size());
  }

  template<typename T>
  inline void UnpackMatrix(mpi::ByteBufferReader& br, Matrix<T>& A)
  {
    size_t m = br.ViewScalarAndAdvance<size_t>();
    size_t n = br.ViewScalarAndAdvance<size_t>();
    A.Resize(m, n);
    br.Read(A.GetPtr(), m*n);
  }

  // Manager for the Parfor.
  class CompressParforMgr : public mpi::ParforManager {
  private:
    Compressor* _c;
    vector<mpi::ByteBufferWriter> _bw;
    ProgressBar _pb;

  public:
    CompressParforMgr(Compressor* c);

    virtual int Isend(int job_idx, int dest);
    virtual int Irecv(int job_idx, int src) { return 0; }
    virtual void IsDone(int job_idx);
  };

  CompressParforMgr::CompressParforMgr(Compressor* c)
    : _c(c), _bw(mpi::GetNproc() - 1),
      _pb("Compress", c->GetBlocks().size(), c->GetOutputLevel())
  {}

  int CompressParforMgr::Isend(int job_idx, int dest)
  {
    /*todo Order blocks from largest to smallest. Group them into tasks. The
      number of blocks in a task is inversely proportional to their size. Then
      job_idx indexes the task rather than the block. */
    int idx = dest - 1;
    const MatBlock* b;
    LraBlock* lb = NULL;
    if (_c->HaveOldHmat()) {
      lb = _c->GetHmatReader()->ReadNextBlock();
      b = &lb->GetBlock();
    } else {
      b = &_c->GetBlocks()[job_idx];
    }
    _bw[idx].Reset();
    _bw[idx].WriteScalar<uint>(b->r0);
    _bw[idx].WriteScalar<uint>(b->m);
    _bw[idx].WriteScalar<uint>(b->c0);
    _bw[idx].WriteScalar<uint>(b->n);
    _bw[idx].WriteScalar<uint>(lb ? lb->GetPrec() : 0);
    if (lb) {
      _bw[idx].WriteScalar((double) lb->GetTol());
      lb->Pack(_bw[idx]);
      delete lb;
    }
    _bw[idx].Isend(dest);
    return 0;
  }

  void CompressParforMgr::IsDone(int job_idx)
  {
    _pb.Incr();
  }

  // Worker for the Parfor.
  class CompressParforWorker : public mpi::ParforWorker {
  private:
    Compressor* _c;
    FILE* _fid;
    mpi::ByteBufferReaderMpi _br;

  public:
    CompressParforWorker(Compressor* c, string fn) throw (FileException);
    virtual ~CompressParforWorker() { fclose(_fid); }

    virtual int Work(int job_idx, int root);
  };

  CompressParforWorker::CompressParforWorker(Compressor* c, string fn)
    throw (FileException)
    : _c(c)
  {
    _fid = fopen(fn.c_str(), "w");
    if (!_fid)
      throw FileException(string("CompressParforWorker: Could not open ") +
                          fn);
  }

  int CompressParforWorker::Work(int job_idx, int root)
  {
    if (mpi::Pid() == 1) Ca::GetTimer()->Tic(1);
    MatBlock b;
    _br.Recv(root);
    b.r0 = _br.ViewScalarAndAdvance<uint>();
    b.m  = _br.ViewScalarAndAdvance<uint>();
    b.c0 = _br.ViewScalarAndAdvance<uint>();
    b.n  = _br.ViewScalarAndAdvance<uint>();
    uint realp = _br.ViewScalarAndAdvance<uint>();
    LraBlock* ob = NULL;
    if (realp) {
      double tol = _br.ViewScalarAndAdvance<double>();
      ob = NewLraBlock(b, realp, tol);
      ob->Unpack(_br);
    }
    if (mpi::Pid() == 1) Ca::GetTimer()->Tic(2);
    _c->CompressBlockToFile(_fid, b, ob);
    if (mpi::Pid() == 1) Ca::GetTimer()->Toc(2);
    if (ob) delete ob;
    if (mpi::Pid() == 1) Ca::GetTimer()->Toc(1);
    return 0;
  }

  bool GreensFnMa::Call(const MatBlock& blk, const vector<uint>* irs,
                        const vector<uint>* ics, double* B)
  {
    vector<uint> rs, cs;
    if (irs) {
      uint n = irs->size();
      rs.resize(n);
      for (uint i = 0; i < n; i++) rs[i] = _pp[(*irs)[i]];
    } else {
      rs.resize(blk.m);
      for (uint i = 0; i < blk.m; i++) rs[i] = _pp[blk.r0 + i];
    }
    if (ics) {
      uint n = ics->size();
      cs.resize(n);
      for (uint i = 0; i < n; i++) cs[i] = _pq[(*ics)[i]];
    } else {
      cs.resize(blk.n);
      for (uint i = 0; i < blk.n; i++) cs[i] = _pq[blk.c0 + i];
    }
    return _gf->Call(rs, cs, B);
  }

  void Compressor::GetHdData(const Hd* hd)
  {
    if (mpi::AmRoot()) {
      hd->Permutations(_ma.pp(), _ma.pq());
      if (mpi::GetNproc() > 1) {
        // Broadcast permutations.
        mpi::ByteBufferWriter bw;
        PackVec(bw, _ma.pp());
        PackVec(bw, _ma.pq());
        bw.Bcast();
      }
      // Only the root stores the blocks.
      _blocks.resize(hd->NbrBlocks());
      uint i = 0;
      for (Hd::iterator it = hd->Begin(), end = hd->End(); it != end; ++it)
        _blocks[i++] = *it;
    } else {
      // Get broadcasted permutations.
      mpi::ByteBufferReaderMpi br;
      br.Bcast();
      UnpackVec(br, _ma.pp());
      UnpackVec(br, _ma.pq());
    }
  }

  static const vector<Blint>&
  CopyAndDecr(const vector<uint>& c, vector<Blint>& d)
  {
    d.resize(c.size());
    for (uint i = 0; i < c.size(); i++) d[i] = c[i] - 1;
    return d;
  }
    
  void Compressor::WriteHmatHeader(FILE* fid)
  {
    // See hm('WriteHeader') for more.
    Blint n;
    n = _ma.pp().size(); write(&n, 1, fid);
    n = _ma.pq().size(); write(&n, 1, fid);
    n = _prec; write(&n, 1, fid);
    // Switch from base-1 indexing (used by Hd) to base-0 indexing (used by
    // Hmat).
    vector<Blint> p;
    write(&CopyAndDecr(_ma.pp(), p)[0], _ma.pp().size(), fid);
    write(&CopyAndDecr(_ma.pq(), p)[0], _ma.pq().size(), fid);
    n = _blocks.size(); write(&n, 1, fid);
    n = 1; /* use abs tol */ write(&n, 1, fid);
    double d;
    d = _tol; write(&d, 1, fid);
    d = 14.0; write(&d, 1, fid);
  }

  Compressor::Compressor(const Hd* hd, GreensFn* gf,
                         const string& filename)
    throw (Exception, FileException)
    : _ma(gf), _filename(filename), _ohm(NULL), _have_old_hmat(false), _prec(2),
      _tm(Compressor::tm_abs), _Bfro(-1.0), _tol(1.0e-6), _output_lvl(1),
      _ablr_stats(105, 0)
  {
    if (!gf) throw Exception("gf is NULL.");
    
    // Let's check now whether we'll have a problem.
    FILE* fid = fopen(_filename.c_str(), "w");
    if (!fid) throw FileException(_filename + string(" can't be written."));
    fclose(fid);

    GetHdData(hd);
  }

  Compressor::~Compressor()
  {
    if (_ohm) delete _ohm;
    if (mpi::Pid() == 1)
      printf
        ("pid %2d:  "
         "old B %d  old QR %d  scratch B %d  scratch SVD %d  scratch ACA %d\n"
         "         B help SVD %d  B help ACA %d  UV help ACA %d\n"
         "         no help B %d  no help SVD %d  no help ACA %d\n",
         mpi::Pid(),
         _ablr_stats[1], _ablr_stats[2], _ablr_stats[0], _ablr_stats[3],
         _ablr_stats[4],
         _ablr_stats[13], _ablr_stats[24], _ablr_stats[34],
         _ablr_stats[100], _ablr_stats[103], _ablr_stats[104]);
  }

  Compressor* NewCompressor
   (const Hd* hd, GreensFn* gf, const string& hmat_filename)
    throw (Exception, FileException)
  {
    return new Compressor(hd, gf, hmat_filename);
  }

  void DeleteCompressor(Compressor* c) { delete c; }

  void Compressor::SetTolMethod(Compressor::TolMethod tm) { _tm = tm; }

  void Compressor::SetTol(double tol)
    throw (Exception)
  {
    if (tol <= 0.0) throw Exception("tol must be >= 0");
    _tol = tol;
  }

  void Compressor::SetOutputLevel(uint lev) { _output_lvl = lev; }

  void Compressor::SetBfroEstimate(double Bfro)
    throw (Exception)
  {
    if (Bfro <= 0.0) throw Exception("Bfro must be >= 0");
    _Bfro = Bfro;
  }

  void Compressor::UseHmatFile(const string& hmat_filename)
    throw (Exception, FileException)
  {
    bool can_read = true, is_compatible = true;
    if (mpi::AmRoot()) {
      try {
        _ohm = new HmatReader(hmat_filename);
      } catch (const FileException& e) {}
      if (!mpi_IsTrue(_ohm))
        can_read = false;
      else {
        // Sort the blocks for use in binary search.
        vector<const MatBlock*> sbs(_blocks.size());
        for (size_t i = 0; i < _blocks.size(); i++)
          sbs[i] = &_blocks[i];
        sort(sbs.begin(), sbs.end(), MatBlockLessThanRC);
        if (!mpi_IsTrue(_ohm->IsCompatible(_ma.pp(), _ma.pq(), sbs)))
          is_compatible = false;
        else {
          _have_old_hmat = true;
          mpi::Bcast(&_have_old_hmat, 1);
        }
      }
    } else {
      can_read = mpi_IsTrue();
      if (can_read) {
        is_compatible = mpi_IsTrue();
        if (is_compatible)
          mpi::Bcast(&_have_old_hmat, 1);
      }
    }
    if (!can_read)
      throw FileException("UseHmatFile: Could not read old H-matrix file.");
    if (!is_compatible)
      throw Exception("old H-matrix is not compatible with new one");
  }

  bool Compressor::HaveOldHmat() const { return _have_old_hmat; }

  double Compressor::GetOldHmatBfro() const throw (Exception)
  {
    double Bfro;
    if (!HaveOldHmat())
      throw Exception("GetOldHmatBfro: No old H-matrix.");
    if (mpi::AmRoot()) {
      Hmat* hm = NewHmat(_ohm->GetFilename());
      Bfro = sqrt(hm->NormFrobenius2());
      DeleteHmat(hm);
      mpi::Bcast(&Bfro, 1);
    } else {
      mpi::Bcast(&Bfro, 1);
    }
    return Bfro;
  }

  static inline bool BlockIsOnDiag(const MatBlock& b)
  {
    return (b.r0 <= b.c0 + b.n - 1 &&
            b.c0 <= b.r0 + b.m - 1);
  }

  static double FullBlockNormFro2(MatrixAccessor& ma, const MatBlock& b)
  {
    Matrix<double> B(b.m, b.n);
    double* pB = B.GetPtr();
    ma.Call(b, NULL, NULL, pB);
    double nf2 = 0.0;
    for (uint i = 0, mn = b.m*b.n; i < mn; i++) nf2 += pB[i]*pB[i];
    return nf2;
  }

  // Get a lower bound on ||B||_F based just on the diag blocks.
  double Compressor::EstimateBfro()
    throw (OutOfMemoryException, UserReqException)
  {
    bool am_root = mpi::AmRoot();

    double nf2_me = 0.0, nf2 = 0.0;
    if (mpi::GetNproc() == 1) { // serial
      ProgressBar pb("Estimate ||B||_F", _blocks.size(), GetOutputLevel());
      for (uint i = 0, nb = _blocks.size(); i < nb; i++) {
        if (BlockIsOnDiag(_blocks[i]))
          nf2 += FullBlockNormFro2(_ma, _blocks[i]);
        pb.Incr();
      }
    } else { // parallel
      if (am_root) {
        // Assign blocks to procs.
        vector<uint> idxs_me;
        uint np = mpi::GetNproc();
        vector< mpi::ByteBufferWriter > _idxs(np - 1);
        for (uint i = 0, k = 0, nb = _blocks.size(); i < nb; i++)
          if (BlockIsOnDiag(_blocks[i])) {
            if (k == (uint) mpi::Root())
              idxs_me.push_back(i);
            else
              _idxs[k - 1].Write(&_blocks[i], 1);
            k = (k + 1) % np;
          }

        // Send out jobs.
        for (uint i = 1; i < np; i++)
          _idxs[i - 1].Isend(i);

        // Do my part.
        ProgressBar pb("Estimate ||B||_F", idxs_me.size(), GetOutputLevel());
        for (size_t i = 0, n = idxs_me.size(); i < n; i++) {
          nf2_me += FullBlockNormFro2(_ma, _blocks[idxs_me[i]]);
          pb.Incr();
        }
        mpi::Barrier(); // Stay in scope until all msgs are delivered.
      } else {
        // Receive my job.
        mpi::ByteBufferReaderMpi br;
        br.Recv(mpi::Root());
        uint nb = br.Size() / sizeof(MatBlock);
        for (uint i = 0; i < nb; i++) {
          const MatBlock* b = br.ViewAndAdvance<MatBlock>(1);
          nf2_me += FullBlockNormFro2(_ma, *b);
        }
        mpi::Barrier(); // Match root's Barrier.
      }
      // Get the total.
      mpi::Allreduce(&nf2_me, &nf2, 1, MPI_SUM);
    }

    if (am_root) printf("est %e\n", sqrt(nf2));
    return sqrt(nf2);
  }

  static void GetTmpFilename(const string& base, string& fn)
  {
    stringstream ss;
    ss << base << "_" << mpi::Pid();
    fn = ss.str();
  }

  // Determine final abs tol and prec from _tol and _Bfro.
  void Compressor::FinalizeTol()
  {
    // If the tolerance involves ||B||_F, include it. In either case, the final
    // _tol is an absolute tolerance on the error ||E||_F.
    if (_tm == tm_Bfro) _tol *= _Bfro;
    // For each block, we need
    //   (tol Bfro)^2 blk.m blk.n / (M N) >= eps(_prec)^2 blk.m blk.n
    //    => (tol Bfro)^2 / (M N) >= eps(_prec)^2
    //    => _tol >= eps(_prec).
    // We use a factor of 10 to be conservative.
    _prec = (_tol < 10 * numeric_limits<float>::epsilon()) ?
      GetPrecCode<double>() : GetPrecCode<float>();
  }

  // Send this many bytes at a time when concat'ing the files. I'm not sure
  // whether one number is much better than any other.
  const uint Compressor::_send_sz = 1L << 24;

  void Compressor::SendFile()
  {
    string fn;
    GetTmpFilename(_filename, fn);
    FILE* fid = fopen(fn.c_str(), "r");

    mpi::ByteBufferWriter bw(sizeof(size_t) + _send_sz);
    for (;;) {
      bw.Reset();
      size_t nread =
        fread(bw.CopyTo<char>(_send_sz), sizeof(char), _send_sz, fid);
      bw.Resize<char>(nread);
      bw.Send(mpi::Root());
      if (nread < _send_sz * sizeof(char)) {
        // That was the last segment.
        if (nread > 0) {
          // Send a final empty message to signal I'm done.
          bw.Reset();
          bw.Send(mpi::Root());
        }
        break;
      }
    }

    fclose(fid);
  }

  void Compressor::
  RecvAndWriteFile(FILE* fid, mpi::ByteBufferReaderMpi& br, uint src)
  {
    for (;;) {
      br.Recv(src);
      if (br.Size() == 0) break;
      write(br.View<char>(0), br.Size(), fid);
    }
  }

  void Compressor::ConcatFiles()
  {
    if (mpi::GetNproc() == 1) return;
    if (mpi::AmRoot()) {
      FILE* fid = fopen(_filename.c_str(), "a");
      if (!fid) throw FileException(_filename + string(" can't be updated."));

      mpi::ByteBufferReaderMpi br;
      for (int pid = 1; pid < mpi::GetNproc(); pid++)
        RecvAndWriteFile(fid, br, pid);

      fclose(fid);
    } else {
      SendFile();
    }
  }

  void Compressor::
  CompressBlockToFile(FILE* fid, const MatBlock& b, const LraBlock* ob)
  {
    // MREM formula for tol.
    double tol = _tol * sqrt((double) (b.m * b.n) /
                             (_ma.pp().size() * _ma.pq().size()));
    LraBlock* lb = NewLraBlock(b, _prec, tol);
    
    LraOptions opts;
    if (mpi::Pid() == 1) Ca::GetTimer()->Tic(3);
    int ret = ApproxByLowRank(opts, lb, _ma, ob);
    if (mpi::Pid() == 1) Ca::GetTimer()->Toc(3);
    _ablr_stats[ret]++;
    lb->WriteToFile(fid);
    delete lb;
  }

  void Compressor::CompressSerial()
  {
    FILE* fid = fopen(_filename.c_str(), "a");
    if (!fid) throw FileException(_filename + string(" can't be updated."));
    ProgressBar pb("Compress", _blocks.size(), GetOutputLevel());
    if (_have_old_hmat) {    
      for (uint i = 0, nb = _blocks.size(); i < nb; i++) {
        LraBlock* ob = _ohm->ReadNextBlock();
        CompressBlockToFile(fid, ob->GetBlock(), ob);
        delete ob;
        pb.Incr();
      }
    } else {
      for (uint i = 0, nb = _blocks.size(); i < nb; i++) {
        CompressBlockToFile(fid, _blocks[i]);
        pb.Incr();
      }
    }
    fclose(fid);
  }

  void Compressor::Compress()
    throw (OutOfMemoryException, UserReqException, FileException)
  {
    bool am_root = mpi::AmRoot();
    FinalizeTol();

    // Write the header.
    if (am_root) {
      FILE* fid = fopen(_filename.c_str(), "w");
      WriteHmatHeader(fid);
      fclose(fid);
    }

    if (mpi::GetNproc() == 1) { // serial
      CompressSerial();
    } else { // parallel
      CompressParforMgr* cpm = NULL;
      CompressParforWorker* cpw = NULL;
      bool all_ok;
      string msg = "";
      if (mpi::AmRoot()) {
        cpm = new CompressParforMgr(this);
        all_ok = mpi::AllOk(cpm);
      }
      else {
        string fn;
        GetTmpFilename(GetFilename(), fn);
        try {
          cpw = new CompressParforWorker(this, fn);
          if (mpi::Pid() == 1) {
            Ca::GetTimer()->Reset(0);
            Ca::GetTimer()->Reset(1);
            Ca::GetTimer()->Reset(2);
            Ca::GetTimer()->Tic(0);
          }
        } catch (const Exception& e) {
          stringstream ss;
          ss <<  "Pid " << mpi::Pid() << " got this when creating a worker: " <<
            e.GetMsg();
          msg = ss.str();
        }
        all_ok = mpi::AllOk(cpw);
      }
      if (all_ok) Parfor(cpm, cpw, _blocks.size());
      if (mpi::Pid() == 1) {
        Timer* t = Ca::GetTimer();
        t->Toc(0);
        caprint("\n-> %e %e %e %e %1.2f\n", t->TotEt(0),
                t->TotEt(1), t->TotEt(2), t->TotEt(3),
                100.0*t->TotEt(3)/t->TotEt(0));
      }
      if (cpm) delete cpm;
      if (cpw) delete cpw;
      if (!all_ok) throw FileException(msg);
      ConcatFiles();
    }
  }
}
