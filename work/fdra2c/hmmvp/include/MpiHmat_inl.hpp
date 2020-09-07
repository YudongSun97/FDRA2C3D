// Routines to parallelize Hmat::Mvp().

#include <stdio.h>
#include <algorithm>
#include "util/include/CodeAnalysis.hpp"
#include "hmmvp/include/Hmat.hpp"
#include "hmmvp/include/HmatIo.hpp"

namespace hmmvp {

  template<typename T> Blint MpiHmat<T>::GetM() const { return hm->GetM(); }
  template<typename T> Blint MpiHmat<T>::GetN() const { return hm->GetN(); }
  template<typename T> Blint MpiHmat<T>::GetNnz() const { return nnz; }
  template<typename T> const Blint* MpiHmat<T>::GetQ() const
  { return hm->GetQ(); }
  template<typename T> const Blint* MpiHmat<T>::GetP() const
  { return hm->GetP(); }

  namespace {

    struct BlockData {
      Blint i, nnz, r0, c0;
    };

    template<bool sort_by_col> bool CmpBlockData
    (const BlockData& bd1, const BlockData& bd2)
    {
      if (sort_by_col) {
        // Sort by col first, but if cols are equal, then sort by row.
        if (bd1.c0 == bd2.c0)
          return bd1.r0 < bd2.r0;
        else
          return bd1.c0 < bd2.c0;
      } else {
        // Opposite of above.
        if (bd1.r0 == bd2.r0)
          return bd1.c0 < bd2.c0;
        else
          return bd1.r0 < bd2.r0;
      }
    }

    Blint SetupMpi
    (uint nproc, const string& filename, vector<vector<Blint> >& block_idxs)
      throw (FileException)
    {
      // Get nnz, row, and column data.
      Blint m, n, realp, nb;
      double tol;
      Blint nnz = 0;
      FILE* fid = fopen(filename.c_str(), "rb");
      hmmvp::ReadHmatHeader(fid, m, n, realp, nb, tol);
      vector<BlockData> bds(nb);
      for (uint i = 0; i < nb; i++) {
        Blint bm, bn, brank;
        BlockData& bd = bds[i];
        bd.i = i;
        hmmvp::ReadHmatBlockInfo(fid, realp, bd.r0, bd.c0, bm, bn, brank);
        if (brank == std::min(bm, bn)) bd.nnz = bm*bn;
        else bd.nnz = (bm + bn)*brank;
        nnz += bd.nnz;
      }
      fclose(fid);

      // Sort the blocks by increasing column. This should encourage locality in
      // x space (in y = B*x in permuted space) in each task.
      std::sort(bds.begin(), bds.end(), CmpBlockData<true>);

      // Allocate blocks among tasks. 'seps' separates bds into groups by task.
      vector<Blint> seps(nproc+1);
      seps[0] = 0;
      for (uint i = 0, cnnz = 0, pid = 0, nnzpt = nnz / nproc + 1;
           i < nb; i++) {
        cnnz += bds[i].nnz;
        if (cnnz > nnzpt) {
          seps[pid+1] = i;
          pid++;
          cnnz = 0;
        }
      }
      seps[nproc] = nb;
      block_idxs.clear();
      block_idxs.resize(nproc);
      for (uint i = 0; i < nproc; i++) {
        // Sort each task's blocks by increasing row. This should encourage
        // cache coherence in y space in each task.
        std::sort(&bds[seps[i]], &bds[seps[i+1]], CmpBlockData<false>);
        // Create block_idxs array.
        block_idxs[i].clear();
        block_idxs[i].resize(seps[i+1] - seps[i]);
        for (uint j = seps[i], k = 0; j < seps[i+1]; j++, k++)
          block_idxs[i][k] = bds[j].i;
      }

      return nnz;
    }

    static const int tag0 = 100;

  }

  template<typename T>
  MpiHmat<T>::MpiHmat(const string& filename, uint ncol)
    throw (FileException)
    : workr(0), am_root(false), do_perms(true)
  {
    int pid = mpi::Pid();

    am_root = mpi::AmRoot();
    vector<Blint> my_block_idxs;
    if (am_root) {
      int nproc = mpi::GetNproc();
      // Determine how to distribute blocks among processes.
      vector<vector<Blint> > block_idxs;
      nnz = SetupMpi(nproc, filename, block_idxs);
      // Assign everyone else their blocks.
      for (int i = 1; i < nproc; i++) {
        Blint nb = block_idxs[i].size();
        mpi::Send(&nb, 1, i, tag0);
        mpi::Send(&block_idxs[i][0], nb, i, tag0);
      }
      // Assign myself my blocks.
      my_block_idxs = block_idxs[0];
    } else {
      // Wait until I receive my block indices.
      Blint nb;
      mpi::Recv(&nb, 1, 0, tag0);
      my_block_idxs.resize(nb);
      mpi::Recv(&my_block_idxs[0], nb, 0, tag0);
    }
  
    // Now read in just my blocks.
    try {
      hm = hmmvp::NewHmat(filename, ncol, 1, &my_block_idxs);
    } catch (const hmmvp::FileException& fe) {
      throw FileException(fe.GetMsg());
    }
    hm->TurnOffPermute();

    uint nw = ncol*std::max(hm->GetM(), hm->GetN());
    workr.resize(nw);
  }

  template<typename T> MpiHmat<T>::~MpiHmat() { delete hm; }

  template<typename T>
  void MpiHmat<T>::Mvp(const T* x, T* y, uint ncol) const
    throw (Exception)
  {
    if (do_perms) {
      // For all but root, x can be NULL. (Not advertised in interface.)
      if (am_root) hm->ApplyQ(x, &workr[0], ncol);
      // For now, send Q*x to all procs.
      mpi::Bcast(&workr[0], ncol*hm->GetN());
      // Do my part.
      hm->Mvp(&workr[0], y, ncol);
      // Sum all the y's.
      mpi::Reduce(y, &workr[0], ncol*hm->GetM(), MPI_SUM);
      if (am_root) {
        // Permute back to user's basis.
        hm->ApplyPt(&workr[0], y, ncol);
      }
    } else {
      // Do my part. (NB: x must be valid in every task, unlike the other case.)
      hm->Mvp(x, &workr[0], ncol);
      // Sum all the y's.
      mpi::Reduce(&workr[0], y, (int) ncol*hm->GetM(), MPI_SUM);
    }
  }

}
