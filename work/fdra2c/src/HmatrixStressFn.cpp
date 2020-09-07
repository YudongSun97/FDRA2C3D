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
    : _as(as), _ncomp(ncomp), _hm(hm) {}
  virtual ~HmatrixStressFn();

  void SetBc (const vector<real>& bc, real v_creep) {
    _bc = bc; _v_creep = v_creep; }
  void SetScale (real scale) { _scale = scale; }

  //try-edc
  void SetBcEdc (const vector<real>& bc_edc) { _bc_edc = bc_edc; }

  virtual bool IncludeNormalComponent () { return false; }
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

  virtual void Mvp(const tau_real* x, tau_real* y) = 0;
};

HmatrixStressFn::~HmatrixStressFn () {
  for (size_t i = 0; i < _hm.size(); i++) delete _hm[i];
}

void HmatrixStressFn::Call (int deriv, double t, const real* x,
                            tau_real* tau, tau_real* taun) {
  const int my_n = _as->GetN();
  if (!tau) return;

  const tau_real* xc;
  if (sizeof(tau_real) == sizeof(real))
    xc = (tau_real*) x;
  else {
    _call_rwrk.Reset(my_n);
    tau_real* xs = _call_rwrk.GetPtr();
    for (int i = 0; i < my_n; i++) xs[i] = (tau_real) x[i];
    xc = _call_rwrk.GetPtr();
  }
  Mvp(xc, tau);

  if (!_bc.empty()) {
    if (deriv == 1) t = 1;
    real s = _scale*t*_v_creep;
    for (int i = 0; i < my_n; i++) tau[i] = _scale*tau[i] + s*_bc[i];
    //try-edc
    if (!_bc_edc.empty())
      for (int i = 0; i < my_n; i++) tau[i] += _bc_edc[i]*x[i];
  } else {
    for (int i = 0; i < my_n; i++) tau[i] *= _scale;
  }
}

class VanillaHsf : public HmatrixStressFn {
public:
  VanillaHsf(const mpi::ArraySegmenter* as, const vector<MpiHmat*>& hm,
             int ncomp);

private:
  virtual void Mvp(const tau_real* x, tau_real* y);
};

VanillaHsf::VanillaHsf (const mpi::ArraySegmenter* as,
                        const vector<MpiHmat*>& hm, int ncomp)
  : HmatrixStressFn(as, hm, ncomp)
{
  int n = _hm[0]->GetM();
  if (mpi::AmRoot()) n *= 2;
  _rwrk.Reset(n);
}

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

class SymDomainHsf : public HmatrixStressFn {
public:
  SymDomainHsf(const mpi::ArraySegmenter* as, const vector<MpiHmat*>& hm,
               int ncomp, const vector<int>& p, const vector<int>& q1,
               const vector<int>& q2);

private:
  vector<int> _p, _q;

  virtual void Mvp(const tau_real* x, tau_real* y);
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

void SymDomainHsf::Mvp (const tau_real* x, tau_real* y) {
  Ca::GetTimer()->Tic(catmr_hsf_mvp);
  _rwrk.Reset();
  const int ntot = _hm[0]->GetN();
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
  _hm[0]->Mvp(x_full_p, y_full_p, 2);
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
    if (vs.SetArray("hm_bc", bc, m->GetNcomp()))
      sf->SetBc(bc, v_creep);
    else {
      delete sf;
      if (mpi::AmRoot())
        fprintf(stderr, "NewHmatrixStressFn: No hm_bc.");
      return NULL;
    }
    if (vs.SetArray("hm_bc_edc", bc, m->GetNcomp()))
      sf->SetBcEdc(bc);
  }

  return sf;
}

}
