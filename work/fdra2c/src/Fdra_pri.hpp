#ifndef INCLUDE_FDRA_FDRA_PRI
#define INCLUDE_FDRA_FDRA_PRI

#include <vector>
#include <list>
#include <string>
#include "util/include/Mpi.hpp"
#include "util/include/KeyValueFile.hpp"
#include "StreamDataFile.hpp"
#include "Defs.hpp"

namespace fdra {
  using namespace std;
  using namespace util;

  class Model;

  class StressFn {
  public:
    virtual ~StressFn() {}

    virtual bool IncludeNormalComponent() = 0;

    // tau is NULL if only the normal component is wanted.
    virtual void Call(int deriv, double t, const real* x,
                      tau_real* tau, tau_real* taun) = 0;
  };

  class StateVector {
  public:
    StateVector(int ncomp) : _ncomp(ncomp) {
      _slip = _p = _T = NULL;
      _v = _theta = _dlte = NULL;
    }

    int Ncomp () const { return _ncomp; }
    const real* slip () const { return _slip; }
    const real* v () const { return _v; }
    const real* theta () const { return _theta; }
    const real* dlte () const { return _dlte; }
    const real* p () const { return _p; }
    const real* T () const { return _T; }

  private:
    friend class OdeFn;
    const real *_slip, *_p, *_T;
    real *_v, *_theta, *_dlte;
    int _ncomp;
  };

  class OutputListener {
  public:
    OutputListener(const mpi::ArraySegmenter* as) : _as(as) {}

    virtual ~OutputListener() {}

    virtual bool Call(double t, double dt, const StateVector* sv) = 0;

  protected:
    const mpi::ArraySegmenter* _as;
  };

  class Model {
  public:
    enum EvolutionLaw { el_Invalid, el_Aging, el_Slip };

  private:
    list<OutputListener*> _ols;

    list<real*> _ptrs;

    struct Control {
      mpi::ArraySegmenter as;
      string stop_indicator;
      int stop_check_frequency;
    } _ctrl;

    struct Ode {
      double rel_tol;
    } _ode;

    struct Mesh {
      int ncomp, nelem, n;
      double ti, tf;
    } _mesh;

    struct Friction {
      EvolutionLaw law;
      bool use_vcutoff, a_of_T;
      real v0, *v1, *v2;
      real *mu0, *a, *b, *d_c;
    } _fr;

    struct MomentumBalance {
      StressFn* stress_fn;
      real eta;
      real *s_normal;
    } _mb;

    struct DiffusionProcesses {
      bool use_p, use_T, use_dilatancy, use_only_theta, use_membrane_diffusion;
      real epsilon, beta, Lambda, rho_cp, h_c, h;
      //todo Later I may switch these to their atomic parts.
      real *t_f_p, *t_f_T;
    } _dp;

    struct InitialConditions {
      real *v, *chi, *slip, *dlt, *p, *T;
    } _ic;

    struct InputOutput {
      bool allow_overwrite;
      string save_filename;
      int v_every, slip_every, state_every;
      int disp_every;

      string lineint_save_filename;
      int lineint_save_every;
      int lineint_nx;
      vector<real> lineint_wts;

      int sample_save_every;
      int sample_tot_size;
      vector<int> sample_idxs;
    } _io;

  public:
    // Public interface.

    int GetN() const;
    int GetNcomp() const;
    int GetNelem() const;

    // Report whether the current state of the model is consistent and
    // sufficient to use in the simulation. It is not safe to call any other
    // method until IsOk returns true.
    bool IsOk(bool disp_messages = true) const;
    // Initialize the listeners only after IsOk() returns true.
    bool InitListeners(bool disp_messages = true);
    // Put contents in a KeyValueFile. Only the root process needs a non-NULL
    // kvf.
    void ToKvf(KeyValueFile* kvf) const;

    void SetStressFn(StressFn* sf);
    StressFn* GetStressFn() const;

    const mpi::ArraySegmenter* GetArraySegmenter() const;

  public:
    // Internal interface.

    Model();
    ~Model();

    int n () const { return _mesh.n; }
    int ncomp () const { return _mesh.ncomp; }
    int nelem () const { return _mesh.nelem; }

    const real* slip_init () const { return _ic.slip; }
    real v_init (int i) const { return _ic.v[i]; }
    real chi_init (int i) const { return _ic.chi[i]; }
    real dlt_init (int i) const { return _ic.dlt[i]; }
    real p_init (int i) const { return _ic.p ? _ic.p[i] : 0; }
    real T_init (int i) const { return _ic.T ? _ic.T[i] : 0; }
    double ti () const { return _mesh.ti; }
    double tf () const { return _mesh.tf; }

    EvolutionLaw GetEvolutionLaw () const { return _fr.law; }
    bool use_vcutoff () const { return _fr.use_vcutoff; }
    real v0 (int i) const { return _fr.v0; }
    real a (int i) const { return _fr.a[i]; }
    real b (int i) const { return _fr.b[i]; }
    real mu0 (int i) const { return _fr.mu0[i]; }
    real d_c (int i) const { return _fr.d_c[i]; }
    real v1 (int i) const { return _fr.v1[i]; }
    real v2 (int i) const { return _fr.v2[i]; }

    StressFn* stress_fn () const { return _mb.stress_fn; }
    real s_normal (int i) const { return _mb.s_normal[i]; }
    real eta (int i) const { return _mb.eta; }

    bool use_membrane_diffusion () const { return _dp.use_membrane_diffusion; }
    bool use_only_theta () const { return _dp.use_only_theta; }
    bool use_dilatancy () const { return _dp.use_dilatancy; }
    bool use_p () const { return _dp.use_p; }
    bool use_T () const { return _dp.use_T; }
    real epsilon (int i) const { return _dp.epsilon; }
    real h (int i) const { return _dp.h; }
    real h_c (int i) const { return _dp.h_c; }
    real beta (int i) const { return _dp.beta; }
    real Lambda (int i) const { return _dp.Lambda; }
    real rho_cp (int i) const { return _dp.rho_cp; }
    real t_f_p (int i) const { return _dp.t_f_p[i]; }
    real t_f_T (int i) const { return _dp.t_f_T[i]; }

    double rel_tol () const { return _ode.rel_tol; }

    const string& stop_indicator () const { return _ctrl.stop_indicator; }
    int stop_check_frequency () const { return _ctrl.stop_check_frequency; }
    int disp_every () const { return _io.disp_every; }

    bool allow_overwrite () const { return _io.allow_overwrite; }
    const string& save_filename () const { return  _io.save_filename; }
    int save_v_every () const { return _io.v_every; }
    int save_slip_every () const { return _io.slip_every; }
    int save_state_every () const { return _io.state_every; }
    const string& lineint_save_filename () const
      { return _io.lineint_save_filename; }
    int lineint_save_every () const { return _io.lineint_save_every; }
    int lineint_nx () const { return _io.lineint_nx; }
    const vector<real>& lineint_wts () const { return _io.lineint_wts; }
    int sample_save_every () const { return _io.sample_save_every; }
    int sample_tot_size () const { return _io.sample_tot_size; }
    const vector<int>& sample_idxs () const { return _io.sample_idxs; }

    const list<OutputListener*>& GetOutputListeners () const { return _ols; }

  private:
    friend Model* BuildModelFromKeyValueFile(const KeyValueFile*, bool);
  };

}

#endif
