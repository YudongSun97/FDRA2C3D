// ncomp = 1 is assumed throughout

#include <stdio.h>
#include <string.h>
#include "util/include/WorkArray.hpp"
#include "util/include/ValueSetter.hpp"
#include "util/include/Mpi.hpp"
#include "MyMpiHmat.hpp"
#include "Fdra.hpp"
#include "MyCodeAnalysis.hpp"
using namespace util;

namespace fdra {

typedef hmmvp::MpiHmat<tau_real> MpiHmat;
typedef hmmvp::Blint Blint;

class HmatrixStressFn : public StressFn {
public:
  enum DomainType { dt_Vanilla,  // y = B x; nothing fancy
                    dt_SymDomain // domain is geometrically symmetric; B is
                    // applied twice
  };
    
  HmatrixStressFn (const mpi::ArraySegmenter* as, const vector<MpiHmat*>& hm,
                   int ncomp)
  //Yudong May 6 2025
  //: _as(as), _ncomp(ncomp), _hm(hm) {}
    : _as(as), _ncomp(ncomp), _hm(hm), _IncludeNormal(false) {}
  virtual ~HmatrixStressFn();

  void SetBc (const vector<real>& bc, real v_creep) {
    _bc = bc; _v_creep = v_creep; }
  void SetScale (real scale) { _scale = scale; }

  //try-edc
  void SetBcEdc (const vector<real>& bc_edc) { _bc_edc = bc_edc; }

  //Yudong May 2 2025 Set boundary condition for normal stressing (loading)
  // virtual bool IncludeNormalComponent () { return false; }
  void SetNormalBc (const vector<real>& bcn, bool IncludeNormal) {
    _bcn = bcn; _IncludeNormal = IncludeNormal; }
  virtual bool IncludeNormalComponent () { return _IncludeNormal; }

  virtual void Call(int deriv, double t, const real* x,
                    tau_real* tau, tau_real* taun);
    
protected:
  const mpi::ArraySegmenter* _as;
  int _ncomp;
  vector<MpiHmat*> _hm;
  vector<real> _bc;
  real _v_creep, _scale;
  WorkArray<tau_real> _call_rwrk, _rwrk;

  //try-edc Funky elastic decoupling experiment.
  vector<real> _bc_edc;

  //Yudong May 2 2025
  vector<real> _bcn;
  bool _IncludeNormal;

  // Yudong May 13 2025
  // To use Green's function for normal stress (Gn), matrix index is 1
  // Green's function for shear stress (Gs) is matrix index 0

  // virtual void Mvp(const tau_real* x, tau_real* y) = 0;
  virtual void Mvp(const tau_real* x, tau_real* y, int matrix_index = 0) = 0;

};

HmatrixStressFn::~HmatrixStressFn () {
  for (size_t i = 0; i < _hm.size(); i++) delete _hm[i];
}

// yudong: This function is used in OdeFn::Call() 
// tau is shear stress and taun is normal stress
void HmatrixStressFn::Call (int deriv, double t, const real* x,
                            tau_real* tau, tau_real* taun) {
  const int my_n = _as->GetN();

  // Yudong May 2 2025, !tau means we are not interested in shear stress, 
  // if (!tau) return;
  
  const tau_real* xc;
  if (sizeof(tau_real) == sizeof(real))
    xc = (tau_real*) x;
  else {
    _call_rwrk.Reset(my_n);
    tau_real* xs = _call_rwrk.GetPtr();
    for (int i = 0; i < my_n; i++) xs[i] = (tau_real) x[i];
    xc = _call_rwrk.GetPtr();
  }

  //Yudong May 2 2025 implement normal stress taun
  if (taun) {
    // Yudong May 6 2025 testing: do not need this
    // Mvp(xc, taun);

    // Add loading term
    if (!_bcn.empty()) {

      // Yudong May 13 2025
      if (_hm.size() == 1) {
        // deriv == 1 means we're doing a time derivative and using it for v.
        if (deriv == 1) t = 1;
        real s = _scale*t*_v_creep;
        // Yudong comments: add loading from boundary condition. _bcn is a Green's function
        //Yudong May 6 2025 testing: taun starts from 0
        // for (int i = 0; i < my_n; i++) taun[i] = _scale*taun[i] + s*_bcn[i];
        for (int i = 0; i < my_n; i++) taun[i] = s*_bcn[i];
      } else if (_hm.size() > 1) {
        Mvp(xc, taun,1);  // Calculate normal stress from H-matrix
        if (deriv == 1) t = 1;
        real s = _scale*t*_v_creep;
        for (int i = 0; i < my_n; i++) taun[i] = _scale*taun[i] + s*_bcn[i];
      }
    } else {
      for (int i = 0; i < my_n; i++) taun[i] *= _scale;
    }
  }
  // Yudong May 2 2025,
  if (tau) {
    Mvp(xc, tau);

    //Add loading term
    if (!_bc.empty()) {
      // deriv == 1 means we're doing a time derivative and using it for v.
      if (deriv == 1) t = 1;
      real s = _scale*t*_v_creep;
      // Yudong comments: add loading from boundary condition. _bc is a Green's function
      for (int i = 0; i < my_n; i++) tau[i] = _scale*tau[i] + s*_bc[i];
      //try-edc
      if (!_bc_edc.empty())
        for (int i = 0; i < my_n; i++) tau[i] += _bc_edc[i]*x[i];
    } else {
      for (int i = 0; i < my_n; i++) tau[i] *= _scale;
    }
  }
}

class VanillaHsf : public HmatrixStressFn {
public:
  VanillaHsf(const mpi::ArraySegmenter* as, const vector<MpiHmat*>& hm,
             int ncomp);

private:
  // Yudong May 13 2025
  // Added matrix_index parameter
  // virtual void Mvp(const tau_real* x, tau_real* y);
  virtual void Mvp(const tau_real* x, tau_real* y, int matrix_index = 0);
};

VanillaHsf::VanillaHsf (const mpi::ArraySegmenter* as,
                        const vector<MpiHmat*>& hm, int ncomp)
  : HmatrixStressFn(as, hm, ncomp)
{
  int n = _hm[0]->GetM();

  // Yudong May 13 2025
  // If we have a second H-matrix, ensure we have enough space for it too
  if (_hm.size() > 1) {
    // Take the maximum size needed by either matrix
    n = std::max(n, (int)_hm[1]->GetM());
  }

  // Yudong May 13 2025
  // if (mpi::AmRoot()) n *= 2;
  if (mpi::AmRoot()) {
    if (_hm.size() > 1) {
      n *= 4;  // For two matrices: x_full, y_full for each matrix
    } else {
      n *= 2;  // Original allocation for one matrix
    }
  }

  _rwrk.Reset(n);
}

// Yudong May 13 2025
// This is the original Mvp function
/*
void VanillaHsf::Mvp (const tau_real* x, tau_real* y) {
  _rwrk.Reset();
  const int ntotc = _hm[0]->GetM();
  bool am_root = mpi::AmRoot();
  // I'm using permutations, so x_full has to be valid only on the root.
  tau_real* x_full = NULL;
  if (am_root) x_full = _rwrk.AllocWork(ntotc);
  _as->Gather(x, x_full, mpi::Root(), _ncomp);
  // y_full has to be allocated on all nodes, but it's valid only on the root.
  tau_real* y_full = _rwrk.AllocWork(ntotc);
  _hm[0]->Mvp(x_full, y_full, 1);
  // y_full is only valid on the root, so scatter the pieces.
  _as->Scatter(y_full, y, mpi::Root(), _ncomp);
}
*/

// Yudong May 13 2025
// This is the new Mvp function with matrix_index parameter
void VanillaHsf::Mvp (const tau_real* x, tau_real* y, int matrix_index) {
  // Check for valid matrix index
  if (matrix_index < 0 || matrix_index >= static_cast<int>(_hm.size())) {
    if (mpi::AmRoot())
      fprintf(stderr, "Invalid matrix index: %d\n", matrix_index);
    return;
  }

  _rwrk.Reset();
  const int ntotc = _hm[matrix_index]->GetM();
  bool am_root = mpi::AmRoot();

  // I'm using permutations, so x_full has to be valid only on the root.
  tau_real* x_full = NULL;
  if (am_root) x_full = _rwrk.AllocWork(ntotc);
  _as->Gather(x, x_full, mpi::Root(), _ncomp);
  
  // y_full has to be allocated on all nodes, but it's valid only on the root.
  tau_real* y_full = _rwrk.AllocWork(ntotc);
  
  // Use the specified matrix index
  _hm[matrix_index]->Mvp(x_full, y_full, 1);
  
  // y_full is only valid on the root, so scatter the pieces.
  _as->Scatter(y_full, y, mpi::Root(), _ncomp);
}
// Yudong May 13 2025

class SymDomainHsf : public HmatrixStressFn {
public:
  SymDomainHsf(const mpi::ArraySegmenter* as, const vector<MpiHmat*>& hm,
               int ncomp, const vector<int>& p, const vector<int>& q1,
               const vector<int>& q2);

private:
  vector<int> _p, _q;

  // Yudong May 13 2025
  // Added matrix_index parameter
  // virtual void Mvp(const tau_real* x, tau_real* y);
  virtual void Mvp(const tau_real* x, tau_real* y, int matrix_index = 0);
};

#define DO_PERMUTE 0

static void ApplyP2P1 (const int* p1b1, const Blint* p2b0, int* pb0, int n) {
#if DO_PERMUTE
  for (int i = 0; i < n; i++) pb0[i] = p1b1[i] - 1;
#else
  for (int i = 0; i < n; i++) pb0[i] = p1b1[(int) p2b0[i]] - 1;
#endif
}

static void ApplyP1P2 (const int* p1b1, const int* p2b0, int* pb0, int n) {
#if DO_PERMUTE
  for (int i = 0; i < n; i++) pb0[i] = p1b1[i] - 1;
#else
  for (int i = 0; i < n; i++) pb0[i] = (int) p2b0[p1b1[i] - 1];
#endif
}

// P2' is applied to each half.
static void ApplyP1P2h (const int* p1b1, const int* p2b0, int* pb0, int n) {
#if DO_PERMUTE
  for (int i = 0; i < n; i++) pb0[i] = p1b1[i] - 1;
#else
  const int n_o_2 = n / 2;
    
  for (int i = 0; i < n_o_2; i++)
    pb0[i] = p1b1[(int) p2b0[i]] - 1;
  for (int i = 0; i < n_o_2; i++)
    pb0[n_o_2 + i] = p1b1[n_o_2 + (int) p2b0[i]] - 1;
#endif
}

SymDomainHsf::SymDomainHsf (
  const mpi::ArraySegmenter* as, const vector<MpiHmat*>& hm, int ncomp,
  const vector<int>& p, const vector<int>& q1, const vector<int>& q2)
  : HmatrixStressFn(as, hm, ncomp)
{

// Yudong May 13 2025
// TODO: implement function for more than one H-matrix for SymDomainHsf
if (_hm.size() > 1) {
  if (mpi::AmRoot())
    fprintf(stderr, "SymDomainHsf: Only one H-matrix is supported for now.\n");
  return;
}


#if !DO_PERMUTE
  _hm[0]->TurnOffPermute();
#endif
  int ntot = _as->GetNtot();
  bool am_root = mpi::AmRoot();

  _rwrk.Reset(3 * _hm[0]->GetN());

  int *q_full, *p_full;
  q_full = p_full = NULL;
  if (am_root) {
    q_full = new int[2 * ntot];
    p_full = new int[ntot];
  }
  _as->Gather(&q1[0], &q_full[0], mpi::Root());
  _as->Gather(&q2[0], &q_full[0] + ntot, mpi::Root());
  _as->Gather(&p[0], &p_full[0], mpi::Root());
  if (am_root) {
    _q.resize(2 * ntot);
    ApplyP2P1(q_full, _hm[0]->GetQ(), &_q[0], ntot);
    ApplyP2P1(q_full + ntot, _hm[0]->GetQ(), &_q[0] + ntot, ntot);

    // Get the inverse permutation of the H-matrix's P.
    int* pi = q_full;
    const Blint* p2 = _hm[0]->GetP();
    int ntot_o_2 = ntot / 2;
    for (int i = 0; i < ntot_o_2; i++) pi[(int) p2[i]] = i;
    // Duplicate it with an offset.
    for (int i = 0; i < ntot_o_2; i++) pi[ntot_o_2 + i] = ntot_o_2 + pi[i];
    _p.resize(ntot);
    // Compose the two permutations.
    ApplyP1P2(p_full, pi, &_p[0], ntot);

    delete[] q_full;
    delete[] p_full;
  }

  Ca::GetTimer()->Reset(20);
}

// Yudong May 13 2025
// TODO: Added matrix_index parameter and check if it work for SymDomainHsf
//void SymDomainHsf::Mvp (const tau_real* x, tau_real* y) {
void SymDomainHsf::Mvp (const tau_real* x, tau_real* y, int matrix_index) {

  Ca::GetTimer()->Tic(catmr_hsf_mvp);
  _rwrk.Reset();

  // Yudong May 13 2025 
  // const int ntot = _hm[0]->GetN();
  const int ntot = _hm[matrix_index]->GetN();
  
  bool am_root = mpi::AmRoot();
  // Gather domain-ordered x onto root. We'll do the whole permutation to
  // H-matrix ordering here. Probably better than parallelizing the
  // permutation, since that entails more comm. I might reexamine this,
  // though.
  tau_real* x_full_p = _rwrk.AllocWork(2 * ntot);
  tau_real* x_full = _rwrk.AllocWork(ntot);
  Ca::GetTimer()->Tic(catmr_hsf_comm);
  _as->Gather(x, x_full, mpi::Root());
  Ca::GetTimer()->Toc(catmr_hsf_comm);
  // Permute two copies of x, each in its own way.
  if (am_root) {
    Ca::GetTimer()->Tic(catmr_hsf_perm);
    for (int ic = 0, os = 0; ic < 2; ic++, os += ntot)
      for (int i = 0; i < ntot; i++)
        x_full_p[os + i] = x_full[_q[os + i]];
    Ca::GetTimer()->Toc(catmr_hsf_perm);
  }
  // Tell everyone the good news.
  Ca::GetTimer()->Tic(catmr_hsf_comm);
  mpi::Bcast(x_full_p, 2 * ntot);
  Ca::GetTimer()->Toc(catmr_hsf_comm);
  // MVP time.
  tau_real* y_full_p = x_full;

  // Yudong May 13 2025
  // _hm[0]->Mvp(x_full_p, y_full_p, 2);
  _hm[matrix_index]->Mvp(x_full_p, y_full_p, 2);

  // Permute back to domain ordering. y has two columns of length ntot/2. Each
  // represents half the domain. Interpret it as one vector of length
  // ntot. _p threads them together.
  tau_real* y_full = NULL;
  if (am_root) {
    Ca::GetTimer()->Tic(catmr_hsf_perm);
    y_full = x_full_p;
    for (int i = 0; i < ntot; i++)
      y_full[i] = y_full_p[_p[i]];
    Ca::GetTimer()->Toc(catmr_hsf_perm);
  }
  // Spread the good word.
  Ca::GetTimer()->Tic(catmr_hsf_comm);
  _as->Scatter(y_full, y, mpi::Root());
  Ca::GetTimer()->Toc(catmr_hsf_comm);
  Ca::GetTimer()->Toc(catmr_hsf_mvp);
}

StressFn* NewHmatrixStressFn (const Model* m, KeyValueFile* kvf) {
  ValueSetter vs(m->GetArraySegmenter(), kvf);

  int symm;
  HmatrixStressFn::DomainType dt;
  if (!vs.SetScalar("hm_symmetric", symm) || !symm)
    dt = HmatrixStressFn::dt_Vanilla;
  else
    dt = HmatrixStressFn::dt_SymDomain;
  int ncol = (dt == HmatrixStressFn::dt_Vanilla) ? 1 : 2;

  string hm_fn;
  if (!vs.SetString("hm_filename", hm_fn)) {
    if (mpi::AmRoot())
      fprintf(stderr, "NewHmatrixStressFn: No hm_filename.\n");
    return NULL;
  }

  vector<MpiHmat*> hm;
  hm.reserve(4);
  try {
    hm.push_back(new MpiHmat(hm_fn + string("_comp11.hmat"), ncol));

    

    if (hm[0]->GetN() != m->GetArraySegmenter()->GetNtot() * m->GetNcomp()) {
      delete hm[0];
      if (mpi::AmRoot())
        fprintf(stderr, "NewHmatrixStressFn: H-matrix size does not agree "
                "with model size.\n");
      return NULL;
    }
    // Yudong May 13 2025
    // try to read in the second H-matrix
    try {
      // read in the second H-matrix (Green's function for normal stress Gn)
      hm.push_back(new MpiHmat(hm_fn + string("_comp1n.hmat"), ncol));

      // Yudong May 13 2025
      // check if the second H-matrix size agrees with model size
      if (hm.size() > 1) {
        if (hm[1]->GetN() != m->GetArraySegmenter()->GetNtot() * m->GetNcomp()) {
          delete hm[1];
          if (mpi::AmRoot())
            fprintf(stderr, "NewHmatrixStressFn: the second H-matrix size does not agree "
                  "with model size.\n");
          return NULL;
        }
      }
    } catch (const Exception& fe2) {
    // Just print a warning and continue with only the first matrix
    if (mpi::AmRoot())
      fprintf(stderr, "Warning: Could not load second H-matrix (Gn for normal stress): %s\nContinuing with only one matrix.\n", fe2.GetMsg().c_str());
    }

  } catch (const Exception& fe) {
    if (mpi::AmRoot())
      fprintf(stderr, "MpiHmat gave this exception: %s\n", fe.GetMsg().c_str());
    return NULL;
  }

  HmatrixStressFn* sf;
  switch (dt) {
  case HmatrixStressFn::dt_Vanilla:
    sf = new VanillaHsf(m->GetArraySegmenter(), hm, m->GetNcomp());
    break;
  case HmatrixStressFn::dt_SymDomain: {
    vector<int> p, q1, q2;
    if (!(vs.SetArray("hm_perm_p", p) && vs.SetArray("hm_perm_q1", q1) &&
          vs.SetArray("hm_perm_q2", q2))) {
      if (mpi::AmRoot())
        fprintf(stderr, "NewHmatrixStressFn: No p and q permutations.\n");
      return NULL;
    }
    sf = new SymDomainHsf(m->GetArraySegmenter(), hm, m->GetNcomp(),
                          p, q1, q2);
    break;
  }}

  { real v_creep, scale;
    vs.SetScalar("v_creep", v_creep, (real) 1.0);
    vs.SetScalar("hm_scale", scale, (real) 1.0);
    sf->SetScale(scale);
    vector<real> bc;
    // Yudong May 2 2025
    vector<real> bcn;
    int ok2;
    bool IncludeNormal;
    vs.SetScalar("IncludeNormal", IncludeNormal, false);

    int uniload, ok;
    vs.SetScalar("uniload", uniload, (int) 0);

    //if uniload==1 (uniform loading), set bc to "tdot_o_vcreep"
    //Later loading will be calculated as: taudot = _bc*_vcreep
    if (uniload) ok=vs.SetArray("tdot_o_vcreep", bc);
    //hm_bc is a (composite) Green's function relating stressing rates to BC slip rates
    //Later loading will be calculated as: taudot = _bc*_vcreep
    else ok=vs.SetArray("hm_bc", bc, m->GetNcomp());
    if (ok)
      sf->SetBc(bc, v_creep);
    else {
      delete sf;
      if (mpi::AmRoot())
        if (uniload) fprintf(stderr, "NewHmatrixStressFn: No tdot_o_vcreep.\n");
        else fprintf(stderr, "NewHmatrixStressFn: No hm_bc.\n");
      return NULL;
    }
    if (vs.SetArray("hm_bc_edc", bc, m->GetNcomp()))
      sf->SetBcEdc(bc);

    // Yudong May 2 2025, Include loading in the normal direction from BC
    if (IncludeNormal) {
      if (uniload) ok2=vs.SetArray("ndot_o_vcreep", bcn);
      else ok2=vs.SetArray("hm_bc_n", bcn, m->GetNcomp());
      if (ok2)
        sf->SetNormalBc(bcn, IncludeNormal);
      else {
        delete sf;
        if (mpi::AmRoot())
          if (uniload) fprintf(stderr, "NewHmatrixStressFn: No ndot_o_vcreep.\n");
          else fprintf(stderr, "NewHmatrixStressFn: No hm_bc_n.\n");
        return NULL;
      }
    }

  }

  return sf;
}

}
