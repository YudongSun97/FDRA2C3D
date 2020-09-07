#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include "util/include/KeyValueFile.hpp"
#include "Fdra_pri.hpp"
#include "MyCodeAnalysis.hpp"
#include "OdeFn.hpp"
using namespace util;

namespace fdra {
void CalcFriction (Model::EvolutionLaw el, real mu0, real a, real b,
                   real v_o_v0, real theta_v0_o_dc,
                   real* mu, real* mu_psi, real* mu_gamma) {
  real g = 0.5 * v_o_v0 * exp((mu0 + b*log(theta_v0_o_dc))/a);
  *mu = a * asinh(g);
  if (mu_psi || mu_gamma) {
    real gasinhg_g = g / sqrt(g*g + 1);
    if (mu_psi)   *mu_psi   = a * gasinhg_g;
    if (mu_gamma) *mu_gamma = b * gasinhg_g;
  }
}

void CalcFriction (Model::EvolutionLaw el, real mu0, real a, real b,
                   real v0, real d_c, real v1, real v2,
                   real v_o_v0, real theta_v0_o_dc,
                   real* mu, real* mu_psi, real* mu_gamma) {
  // In v-cutoff friction,
  //     psibar   = -log(v1 / v         + 1);
  //     gammabar =  log(v2 theta / d_c + 1);
  // Do some conversions and then use the regular CalcFriction routine.
  real v_o_v0_standin = 1.0 / (v1 / (v0 * v_o_v0) + 1.0);
  real theta_v0_o_dc_standin = v2 * theta_v0_o_dc / v0 + 1;
  CalcFriction(el, mu0, a, b, v_o_v0_standin, theta_v0_o_dc_standin,
               mu, mu_psi, mu_gamma);
  // At this point, mu_psi and mu_gamma are actually mu_psibar and
  // mu_gammabar. Convert.
  *mu_psi *= 1.0 / (1.0 + v0 * v_o_v0 / v1);
  *mu_gamma *= theta_v0_o_dc / (theta_v0_o_dc + v0 / v2);
}

// theta_dot / theta.
void EvolveState (Model::EvolutionLaw el, real v0, real d_c,
                  real v_o_v0, real theta_v0_o_dc,
                  real* thetad_o_theta) {
  switch (el) {
  case Model::el_Aging: {
    *thetad_o_theta = v0/d_c * (1.0/theta_v0_o_dc - v_o_v0);
    break;
  }
  case Model::el_Slip: {
    real echi = v_o_v0 * theta_v0_o_dc;
    *thetad_o_theta = -echi * log(echi) *(v0 / (d_c * theta_v0_o_dc));
    break;
  }
  default:
    assert(false);
  }
}

class MembraneDiffusionIntegrator {
public:
  MembraneDiffusionIntegrator(const Model* m);
                 
  void DoStage(int i, double dt, real bes, real v_o_v0, real phi_dot_wo_seff,
               real mu, real* pd);
  void FinishStep();

  // Valid after DoStage().
  real p (int i) const { return _p1[i]; }
  real T (int i) const { return _T1[i]; }

  const real* p () const { return _p1.empty() ? NULL : &_p1[0]; }
  const real* T () const { return _T1.empty() ? NULL : &_T1[0]; }

private:
  const Model* _m;
  vector<real> _p0, _p1, _T0, _T1;
};

MembraneDiffusionIntegrator::MembraneDiffusionIntegrator (const Model* m)
  : _m(m)
{
  int nelem = _m->nelem();
  if (m->use_p()) {
    _p0.resize(nelem);
    _p1.resize(nelem);
    for (int i = 0; i < nelem; i++) _p0[i] = _m->p_init(i);
  }
  if (m->use_T()) {
    _T0.resize(nelem);
    _T1.resize(nelem);
    for (int i = 0; i < nelem; i++) _T0[i] = _m->T_init(i);
  }
}

void MembraneDiffusionIntegrator::
DoStage (int i, double dt, real bes, real v_o_v0, real phi_dot_wo_seff,
         real mu, real* pd) {
  if (_m->t_f_p(i) < 0) return;
  real pf, Tf = 0;
  pf = -0.5 / _m->beta(i) * phi_dot_wo_seff / bes;
  const bool do_only_p = _m->Lambda(i) == 0 || _m->t_f_T(i) < 0;
  if (do_only_p) {
    _p1[i] = (_p0[i] + dt * bes * pf) / (1 + dt * (_m->t_f_p(i) + pf));
  } else {
    Tf = 0.5 * mu * (v_o_v0 * _m->v0(i)) / (_m->rho_cp(i) * _m->h(i));
    double A[4], b[2];
    A[0] = 1 + dt * (_m->t_f_p(i) + pf);
    A[1] = dt * Tf;
    A[2] = -_m->Lambda(i);
    A[3] = 1 + dt * _m->t_f_T(i);
    b[0] = _p0[i] - _m->Lambda(i) * _T0[i] + dt * bes * pf;
    b[1] = _T0[i] + dt * bes * Tf;
    double det = A[0]*A[3] - A[1]*A[2];
    _p1[i] = (A[3]*b[0] - A[2]*b[1]) / det;
    _T1[i] = (A[0]*b[1] - A[1]*b[0]) / det;
  }
  if (pd) {
    *pd = _m->t_f_p(i) * -_p1[i] + pf * (bes - _p1[i]);
    if (!do_only_p)
      *pd += _m->Lambda(i) * (_m->t_f_T(i) * -_T1[i] + (bes - _p1[i]) * Tf);
  }
}

void MembraneDiffusionIntegrator::FinishStep () {
  _p0 = _p1;
  _T0 = _T1;
}

// State vector is [slip (v / v0) (theta v0 / d_c) [exp(dlt) v0 / d_c]]

OdeStateVector::OdeStateVector (const Model* m) {
  _nelem = m->nelem();
  _n = m->n();
  _use_only_theta = m->use_only_theta();
  _use_dilatancy = m->use_dilatancy();
}

bool OdeStateVector::Decode (const real* py) {
  slip = py;
  v_o_v0 = py + _n;
  theta_v0_o_dc = py + 2*_n;
  dlte_v0_o_dc = NULL;
  if (_use_dilatancy)
    dlte_v0_o_dc = _use_only_theta ? theta_v0_o_dc : py + 3*_n;
  return true;
}

OdeFn::OdeFn (const Model* m) : _m(m), _osv(m), _mdi(NULL), _ctr(0) {
  if (!m->stop_indicator().empty()) {
    _si.check = m->stop_check_frequency() > 0;
    _si.indicator = m->stop_indicator();
    _si.frequency = m->stop_check_frequency();
  }
  _disp_every = m->disp_every();
  if (m->use_membrane_diffusion()) _mdi = new MembraneDiffusionIntegrator(m);

  // [A1] Design decision: If using dilatancy, then evolve porosity separately
  // from theta.
  assert(!_m->use_dilatancy() || !_m->use_only_theta());
}

OdeFn::~OdeFn () { if (_mdi) delete _mdi; }

int OdeFn::GetN () {
  return 2*_m->n() + _m->nelem() + (_m->use_only_theta() ? 0 : _m->nelem());
}

void OdeFn::GetIc (real* y) {
  int nelem = _m->nelem(), n = _m->n();
  memcpy(y, _m->slip_init(), n*sizeof(real));
  real* v_o_v0 = y + n;
  real* theta_v0_o_dc = v_o_v0 + n;
  for (int i = 0; i < n; i++) {
    v_o_v0[i] = _m->v_init(i) / _m->v0(i);
    theta_v0_o_dc[i] = exp(_m->chi_init(i)) / v_o_v0[i];
  }
  if (!_m->use_only_theta()) {
    real* dlte = theta_v0_o_dc + nelem;
    for (int i = 0; i < n; i++) dlte[i] = exp(_m->dlt_init(i));
  }
}

//todo I could determine whether any listener wants the StateVector before
// forming it.
bool OdeFn::ViewOutput (double t, double dt, const real* y, Message m) {
  const int n = _m->n(), nelem = _m->nelem();
  _osv.Decode(y);

  if (_disp_every > 0 && _ctr % _disp_every == 0) {
    if (mpi::AmRoot()) {
      double et = Ca::GetTimer()->Toc(catmr_viewoutput);
      fprintf(stdout, "%8d %1.6e %1.3e %1.3e\n",
              _ctr, t, dt, et);
      caprint("ode %4.1f tau %4.1f\n"
              "hsfmvp %4.1f = hsfcomm %4.1f + hsfperm %4.1f + "
              "mpimvpred %4.1f + hmatmvp %4.1f\n",
              100*Ca::GetTimer()->TotEt(catmr_odefn)/et,
              100*Ca::GetTimer()->TotEt(catmr_odefn_tau)/et,
              100*Ca::GetTimer()->TotEt(catmr_hsf_mvp)/et,
              100*Ca::GetTimer()->TotEt(catmr_hsf_comm)/et,
              100*Ca::GetTimer()->TotEt(catmr_hsf_perm)/et,
              100*Ca::GetTimer()->TotEt(catmr_mpimvp_reduce)/et,
              100*Ca::GetTimer()->TotEt(catmr_hmat_mvp)/et);
    }

#ifdef ANALYZE_CODE
    vector<double> hmat_times;
    if (mpi::AmRoot()) hmat_times.resize(mpi::GetNproc());
    double hmat_time = Ca::GetTimer()->TotEt(catmr_hmat_mvp);
    mpi::Gather(&hmat_time, 1, &hmat_times[0], 1);
    if (mpi::AmRoot()) {
      double mit = hmat_time, mat = hmat_time;
      caprint("mvp: ");
      for (int i = 0; i < mpi::GetNproc(); i++) {
        if (mit > hmat_times[i]) mit = hmat_times[i];
        if (mat < hmat_times[i]) mat = hmat_times[i];
        //caprint("%1.2e ", hmat_times[i]);
      }
      caprint("(%4.1f)\n", 100*mit/mat);
    }
#endif

    Ca::GetTimer()->Tic(catmr_viewoutput);
    Ca::GetTimer()->Reset(catmr_odefn);
    Ca::GetTimer()->Reset(catmr_odefn_tau);
    Ca::GetTimer()->Reset(catmr_hsf_mvp);
    Ca::GetTimer()->Reset(catmr_hsf_comm);
    Ca::GetTimer()->Reset(catmr_hsf_perm);
    Ca::GetTimer()->Reset(catmr_mpimvp_reduce);
    Ca::GetTimer()->Reset(catmr_hmat_mvp);
  }

  if (_mdi) _mdi->FinishStep();

  StateVector sv(_m->ncomp());
  _rwrk.Reset(n + nelem + (_m->use_only_theta() ? 0 : nelem));
  sv._slip = _osv.slip;
  sv._v = _rwrk.AllocWork(n);
  for (int i = 0; i < n; i++)
    sv._v[i] = _m->v0(i) * _osv.v_o_v0[i];
  sv._theta = _rwrk.AllocWork(nelem);
  for (int i = 0; i < nelem; i++)
    sv._theta[i] = _m->d_c(i) * _osv.theta_v0_o_dc[i] / _m->v0(i);
  if (!_m->use_only_theta()) {
    sv._dlte = _rwrk.AllocWork(nelem);
    for (int i = 0; i < nelem; i++)
      sv._dlte[i] = _m->d_c(i) * _osv.dlte_v0_o_dc[i] / _m->v0(i);
  }
  if (_mdi) {
    sv._p = _mdi->p();
    sv._T = _mdi->T();
  }

  typedef list<OutputListener*> OLs;
  const OLs& ols = _m->GetOutputListeners();
  bool ret = true;
  for (OLs::const_iterator it = ols.begin(), end = ols.end();
       it != end; ++it)
    ret = (*it)->Call(t, dt, &sv) && ret;

  if (_si.check) {
    if (_ctr % _si.frequency == 0 && !_si.indicator.empty()) {
      bool cont = true;
      if (mpi::AmRoot())
        cont = access(_si.indicator.c_str(), F_OK) != 0;
      ret = mpi::IsTrue(cont);
    }
  }

  _ctr++;

  if (m == OdeFunction::m_tolfail && mpi::AmRoot())
    fprintf(stderr, "Ode tol failure.\n");

  if (!ret && _disp_every > 0 && mpi::AmRoot())
    fprintf(stdout, "stop requested\n");
  return ret;
}

void OdeFn::ResetY (real* y) {
}

static void SetTols (const Model* m, Ode23& ode) {
  ode.SetRelTol(m->rel_tol());
  ode.SetAbsTol(1.0e-40);
}

void RunOde (const Model* m) {
  OdeFn ofn(m);
  Ode23 ode;
  ode.SetTspan(m->ti(), m->tf());
  SetTols(m, ode);
  try {
    ode.Run(ofn);
  } catch (const OdeException& oe) {
    fprintf(stderr, "%s\n", oe.GetMsg().c_str());
  }
}

static inline void
CalcFriction (const Model* m, const OdeStateVector& osv, const int i,
              real* mu, real* mu_psi, real* mu_gamma) {
  if (m->use_vcutoff())
    CalcFriction(m->GetEvolutionLaw(), m->mu0(i), m->a(i), m->b(i),
                 m->v0(i), m->d_c(i), m->v1(i), m->v2(i),
                 osv.v_o_v0[i], osv.theta_v0_o_dc[i],
                 mu, mu_psi, mu_gamma);
  else
    CalcFriction(m->GetEvolutionLaw(), m->mu0(i), m->a(i), m->b(i),
                 osv.v_o_v0[i], osv.theta_v0_o_dc[i],
                 mu, mu_psi, mu_gamma);
}

void OdeFn::
Call (double t, double dt, const real* y, real* yd, bool& is_error) {
  Ca::GetTimer()->Tic(catmr_odefn);

  const int n = _m->n(), nelem = _m->nelem();
  _osv.Decode(y);

  _rwrk.Reset(n);
  real* v = _rwrk.AllocWork(n);
  for (int i = 0; i < n; i++) v[i] = _osv.v_o_v0[i] * _m->v0(i);

  Ca::GetTimer()->Tic(catmr_odefn_tau);
  const bool inc = _m->stress_fn()->IncludeNormalComponent();
  tau_real *taun = NULL, *taund = NULL;
  if (inc) {
    _trwrk.Reset(n + 2*nelem);
    taun  = _trwrk.AllocWork(nelem);
    taund = _trwrk.AllocWork(nelem);
    _m->stress_fn()->Call(0, t, _osv.slip, NULL, taun);
  } else {
    _trwrk.Reset(n);
  }
  tau_real* taud = _trwrk.AllocWork(n);
  _m->stress_fn()->Call(1, t, v, taud, taund);
  Ca::GetTimer()->Toc(catmr_odefn_tau);
  
  real* slip_d = yd;
  real* v_o_v0_d = yd + n;
  real* theta_v0_o_dc_d = yd + 2*n;
  real* dlte_v0_o_dc_d = _m->use_only_theta() ? NULL : yd + 3*n;
  for (int i = 0; i < nelem; i++) {
    real bes, es, esd, mu, mu_psi, mu_gamma, thetad_o_theta,
      dlted_o_dlt_wo_seff = 0;
    CalcFriction(_m, _osv, i, &mu, &mu_psi, &mu_gamma);
    es = _m->s_normal(i);
    esd = 0.0;
    if (inc) {
      es -= taun[i];
      esd = -taund[i];
    }
    bes = es;
    if (_m->use_dilatancy()) { // See assert [A1].
      EvolveState(_m->GetEvolutionLaw(), _m->v0(i), _m->d_c(i), _osv.v_o_v0[i],
                  _osv.dlte_v0_o_dc[i], &dlted_o_dlt_wo_seff);
      dlted_o_dlt_wo_seff *= -_m->epsilon(i) * _m->h_c(i) / _m->h(i);
    }
    if (_mdi) {
      double pd;
      _mdi->DoStage(i, dt, es, _osv.v_o_v0[i], dlted_o_dlt_wo_seff, mu, &pd);
      es -= _mdi->p(i);
      esd -= pd;
    }
    slip_d[i] = v[i];
    EvolveState(_m->GetEvolutionLaw(), _m->v0(i), _m->d_c(i), _osv.v_o_v0[i],
                _osv.theta_v0_o_dc[i], &thetad_o_theta);
    v_o_v0_d[i] = _osv.v_o_v0[i] *
      (taud[i] - mu*esd - es*mu_gamma*thetad_o_theta) /
      (mu_psi*es + _m->eta(i)*v[i]);
    theta_v0_o_dc_d[i] = _osv.theta_v0_o_dc[i] * thetad_o_theta;
    if (_m->use_dilatancy())
      dlte_v0_o_dc_d[i] = (_osv.dlte_v0_o_dc[i] * (es / bes) *
                           dlted_o_dlt_wo_seff);
  }

  Ca::GetTimer()->Toc(catmr_odefn);
}

}
