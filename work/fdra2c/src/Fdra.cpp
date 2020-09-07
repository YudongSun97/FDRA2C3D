#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <limits>
#include "util/include/Exception.hpp"
#include "util/include/ValueSetter.hpp"
#include "util/include/KeyValueFile.hpp"
#include "StreamDataFile.hpp"
#include "OdeFn.hpp"
#include "Fdra_pri.hpp"
#include "MyCodeAnalysis.hpp"

namespace fdra {
using namespace util;
  
class SpringSliderStressFn : public StressFn {
public:
  SpringSliderStressFn (const vector<real>& k, real v_creep)
    : _k(k), _v_creep(v_creep) {}

  virtual bool IncludeNormalComponent () { return false; }

  virtual void Call (int deriv, double t, const real* x,
                     tau_real* tau, tau_real* taun);

private:
  int _nrep;
  vector<real> _k;
  real _v_creep;
};

void SpringSliderStressFn::Call (int deriv, double t, const real* x,
                                 tau_real* tau, tau_real* taun) {
  for (int i = 0, n = _k.size(); i < n; i++) {
    if (deriv == 0)
      tau[i] = _k[i]*(_v_creep*t - x[i]);
    else
      tau[i] = _k[i]*(_v_creep - x[i]);
  }
}

static SpringSliderStressFn* NewSpringSliderStressFn (ValueSetter& vs) {
  vector<real> k;
  real v_creep;
  if (!(vs.SetArray("ss_k", k) && vs.SetScalar("v_creep", v_creep)))
    return NULL;
  return new SpringSliderStressFn(k, v_creep);
  exit(-1);
}

namespace ol {
  struct DefaultSdfs { enum Enum { slip = 0, v, theta, dlte, p, T }; };

  void InitSdfs (const Model* m, vector<StreamDataFile*>& sdfs,
                 const string& fnpre, int ntot, int ncomp) {
    bool all_ok = true, already_exists = false;
    if (mpi::AmRoot()) {
      static const char* fns[] = {"slip", "v", "theta", "dlte", "p", "T"};

      int nrs[6];
      nrs[DefaultSdfs::slip] = ntot*ncomp; nrs[DefaultSdfs::v] = ntot*ncomp;
      for (int j = 2; j < 6; j++) nrs[j] = ntot;

      bool write_fns[] = {true, true, true, false, false, false};
      write_fns[DefaultSdfs::dlte] = !m->use_only_theta();
      write_fns[DefaultSdfs::p] = m->use_T() || m->use_p();
      write_fns[DefaultSdfs::T] = m->use_T();

      sdfs.resize(12, NULL);
      for (size_t i = 0; i < sizeof(fns) / sizeof(char*); i++) {
        string fn = fnpre + string("_") + string(fns[i]) + string(".sdf");
        if (!m->allow_overwrite() && access(fn.c_str(), F_OK) == 0) {
          for (size_t j = 0; j < sdfs.size(); j++) DeleteStreamDataFile(sdfs[j]);
          already_exists = true;
          break;
        }
        if (!write_fns[i]) continue;
        StreamDataFile *s1, *s2;
        sdfs[2*i] =
          (s1 = NewStreamDataFileForWriting
           (fnpre + string("_") + string(fns[i]) + string(".tsdf"), 1,
            StreamDataFile::p_double));
        sdfs[2*i+1] =
          (s2 = NewStreamDataFileForWriting
           (fn, nrs[i] + 1, StreamDataFile::p_single));
        if (!(s1 && s2)) {
          for (size_t j = 0; j < sdfs.size(); j++) DeleteStreamDataFile(sdfs[j]);
          all_ok = false;
          break;
        }
      }
    }
    if (mpi::IsTrue(already_exists)) throw FileException("OL: File exists.");
    if (!mpi::IsTrue(all_ok)) throw FileException("OL: Can't init sdf.");
  }

  inline bool CheckCtr (int ctr, int save_every) {
    return save_every > 0 && ctr % save_every == 0;
  }
}

class DefaultGatherOL : public OutputListener {
public:
  DefaultGatherOL(const Model* m) throw (FileException);
  virtual ~DefaultGatherOL();

  virtual bool Call(double t, double dt, const StateVector* sv);

private:
  WorkArray<real> _rwrk;
  vector<StreamDataFile*> _sdfs;
  int _v_every, _slip_every, _state_every, _t_every;
  long long int _ctr;

private:
  void Write(int si, double t, int ntot, int factor);
};

DefaultGatherOL::DefaultGatherOL (const Model* m) throw (FileException)
  : OutputListener(m->GetArraySegmenter()), _ctr(0)
{
  _v_every = m->save_v_every();
  _slip_every = m->save_slip_every();
  _state_every = m->save_state_every();
  ol::InitSdfs(m, _sdfs, m->save_filename(), _as->GetNtot(), m->ncomp());
}

DefaultGatherOL::~DefaultGatherOL () {
  for (size_t i = 0; i < _sdfs.size(); i++) DeleteStreamDataFile(_sdfs[i]);
}

inline void DefaultGatherOL::Write (int si, double t, int ntot, int factor) {
  _sdfs[si]->Write(&t, 1);
  _sdfs[si]->Flush();
  _sdfs[si+1]->Write(&t, 1);
  _sdfs[si+1]->Write(_rwrk.GetPtr(), ntot*factor);
  _sdfs[si+1]->Flush();
}

bool DefaultGatherOL::Call (double t, double dt, const StateVector* sv) {
  int ntot = _as->GetNtot(), ncomp = sv->Ncomp();
  bool am_root = mpi::AmRoot();
  if (am_root) _rwrk.Reset(ntot*ncomp);

#define _Write_(fld) do {                                               \
    _as->Gather(sv->fld(), _rwrk.GetPtr(), mpi::Root(), ncomp);         \
    if (am_root) Write(2 * ol::DefaultSdfs::fld, t, ntot, ncomp);       \
  } while (0)

  if (ol::CheckCtr(_ctr, _slip_every)) _Write_(slip);
  if (ol::CheckCtr(_ctr, _v_every)) _Write_(v);
  if (ol::CheckCtr(_ctr, _state_every)) {
    _Write_(theta);
    if (sv->dlte()) _Write_(dlte);
    if (sv->p()) _Write_(p);
    if (sv->T()) _Write_(T);
  }

#undef _Write_
    
  _ctr++;
  return true;
}

class LineintOL : public OutputListener {
public:
  LineintOL(const Model* m) throw (FileException);
  virtual ~LineintOL();

  virtual bool Call(double t, double dt, const StateVector* sv);

private:
  const vector<real>& _w;
  WorkArray<real> _rwrk;
  StreamDataFile* _sdf;
  int _every, _nzc, _nx, _ncomp, _z_offset;
  long long int _ctr;
};

LineintOL::LineintOL (const Model* m) throw (FileException)
  : OutputListener(m->GetArraySegmenter()), _w(m->lineint_wts()), _sdf(NULL),
    _ctr(0)
{
  _every = m->lineint_save_every();
  _ncomp = m->ncomp();
  int ntot = _as->GetNtot(), nz = ntot / m->lineint_nx();
  _nzc = _ncomp * nz;
  _nx = ntot / nz;
  _z_offset = _as->GetOffset() % nz;

  bool all_ok = true, already_exists = false;
  if (mpi::AmRoot()) {
    string fn = m->save_filename() + string("_lix.sdf");
    if (!m->allow_overwrite() && access(fn.c_str(), F_OK) == 0) {
      already_exists = true;
    } else {
      _sdf = NewStreamDataFileForWriting(fn, _nzc + 1,
                                         StreamDataFile::p_double);
      all_ok = _sdf;
    }
  }
  if (mpi::IsTrue(already_exists))
    throw FileException("LineintOL: File exists.");
  if (!mpi::IsTrue(all_ok))
    throw FileException("LineintOL: Can't init sdf.");
}
  
LineintOL::~LineintOL () { if (_sdf) DeleteStreamDataFile(_sdf); }

bool LineintOL::Call (double t, double dt, const StateVector* sv) {
  if (ol::CheckCtr(_ctr, _every)) {
    bool am_root = mpi::AmRoot();
    int my_n = _as->GetN();
    _rwrk.Reset(am_root ? 2 * _nzc : _nzc);
    real* lix = _rwrk.AllocWork(_nzc);
    memset(lix, 0, _nzc*sizeof(real));
    const real* w = &_w[0];
    const real* v = sv->v();
    const int nzc = _nzc * _ncomp;
    for (int i = 0, ix = 0, iz = _z_offset; i < my_n; i++) {
      for (int ic = 0; ic < _ncomp; ic++)
        lix[iz + ic] += w[ix] * v[i];
      iz++;
      if (iz == nzc) {
        iz = 0;
        ix++;
      }
    }
    real* rcv = NULL;
    if (am_root) rcv = _rwrk.AllocWork(_nzc);
    mpi::Reduce(lix, rcv, _nzc, MPI_SUM);
    if (am_root) {
      _sdf->Write(&t, 1);
      _sdf->Write(rcv, _nzc);
    }
  }
  _ctr++;
  return true;
}

class SampleOL : public OutputListener {
public:
  SampleOL(const Model* m) throw (FileException);
  virtual ~SampleOL();

  virtual bool Call(double t, double dt, const StateVector* sv);

private:
  WorkArray<real> _rwrk;
  vector<StreamDataFile*> _sdfs;
  vector<int> _idxs;
  int _every;
  long long int _ctr;
  mpi::ArraySegmenter as;
};

SampleOL::SampleOL (const Model* m) throw (FileException)
  : OutputListener(m->GetArraySegmenter()), _idxs(m->sample_idxs()), _ctr(0)
{
  _every = m->sample_save_every();
  ol::InitSdfs(m, _sdfs, m->save_filename() + "_smpl", m->sample_tot_size(),
               m->ncomp());
  as.ApportionToMe(m->sample_idxs().size());
}

SampleOL::~SampleOL () {
  for (size_t i = 0; i < _sdfs.size(); i++) DeleteStreamDataFile(_sdfs[i]);
}

bool SampleOL::Call (double t, double dt, const StateVector* sv) {
  if (ol::CheckCtr(_ctr, _every)) {
    //here
  }
  return false;
}

Model::Model () {
  _mb.stress_fn = NULL;
  _mesh.nelem = -1;
  _fr.a_of_T = false;
  _ode.rel_tol = 1.0e-3;
}

Model::~Model () {
  for (list<real*>::iterator it = _ptrs.begin(), end = _ptrs.end();
       it != end; ++it)
    delete[] *it;
  for (list<OutputListener*>::iterator it = _ols.begin(), end = _ols.end();
       it != end; ++it)
    delete *it;
  if (_mb.stress_fn) delete _mb.stress_fn;
}

bool Model::IsOk (bool disp) const {
#define require(a)                                                      \
  (!(a) ? (printf("IsOk: (" #a ") does not hold.\n"), false) : true)

  require(_mesh.nelem > 0);
  require(_mesh.ncomp > 0);

  require(_ode.rel_tol > 0);

  require(_fr.law != el_Invalid);
  require(_fr.mu0);
  require(_fr.a);
  require(_fr.b);
  require(_fr.d_c);
  require(!_fr.use_vcutoff || (_fr.v1 && _fr.v2));

  require(_ic.v);
  require(_ic.chi);
  require(_ic.slip);

  require(_dp.use_dilatancy || _dp.use_only_theta);
  require(_dp.use_only_theta == (_ic.dlt == NULL));

  require(_mb.eta > 0.0);
  require(_mb.stress_fn);
  require(_mb.s_normal);

  require(_mesh.ti >= 0.0);
  require(_mesh.tf > _mesh.ti);

  require(!_io.save_filename.empty());
  require(_io.v_every > 0 || _io.slip_every > 0 || _io.state_every > 0);
  if (!_io.lineint_save_filename.empty()) {
    require(!_io.lineint_wts.empty());
    int lineint_div =_ctrl.as.GetNtot() % _io.lineint_nx;
    require(lineint_div == 0);
  }
  if (_io.sample_tot_size > 0) {
    require(_io.sample_save_every > 0);
    for (size_t i = 0; i < _io.sample_idxs.size(); i++)
      require(_io.sample_idxs[i] >= 0 && _io.sample_idxs[i] < _mesh.n);
  }

  require(!_dp.use_dilatancy || !_dp.use_only_theta);
  if (_dp.t_f_p) {
    // For now, if we're using t_f_p, then we're using dilatancy.
    require(_dp.use_membrane_diffusion && _dp.use_dilatancy);
    require(!_dp.Lambda > 0 || _dp.t_f_T);
    require(!_dp.t_f_p || _ic.dlt);
  }

  // For now.
  require(_mesh.ncomp == 1);
  require(!_fr.a_of_T);

  return true;
#undef require
}

bool Model::InitListeners (bool disp_messages) {
  if (!_io.save_filename.empty()) {
    try {
      OutputListener* ol = new DefaultGatherOL(this);
      _ols.push_back(ol);
    } catch (const FileException& fe) {
      if (disp_messages && mpi::AmRoot())
        fprintf(stderr, "InitListeners DefaultGatherOL: %s\n",
                fe.GetMsg().c_str());
      return false;
    }
  } else {
    if (disp_messages && mpi::AmRoot())
      fprintf(stderr, "InitListeners: No save_filename.");
    return false;
  }
  if (!_io.lineint_save_filename.empty()) {
    try {
      OutputListener* ol = new LineintOL(this);
      _ols.push_back(ol);
    } catch (const FileException& fe) {
      if (disp_messages && mpi::AmRoot())
        fprintf(stderr, "InitListeners LineintOL: %s\n", fe.GetMsg().c_str());
      return false;
    }
  }
  if (_io.sample_save_every > 0) {
    try {
      OutputListener* ol = new SampleOL(this);
      _ols.push_back(ol);
    } catch (const FileException& fe) {
      if (disp_messages && mpi::AmRoot())
        fprintf(stderr, "InitListeners SampleOL: %s\n", fe.GetMsg().c_str());
      return false;
    }
  }
  return true;
}

int Model::GetN () const { return _mesh.n; }
int Model::GetNcomp () const { return _mesh.ncomp; }
int Model::GetNelem () const { return _mesh.nelem; }
void Model::SetStressFn(StressFn* sf) { _mb.stress_fn = sf; }
StressFn* Model::GetStressFn () const { return _mb.stress_fn; }
const mpi::ArraySegmenter* Model::GetArraySegmenter () const
{ return &_ctrl.as; }

class ValuePutter {
public:
  ValuePutter(const mpi::ArraySegmenter* as, KeyValueFile* kvf);

  void PutArray(const string& field, const real* v, int factor = 1);

private:
  const mpi::ArraySegmenter* _as;
  KeyValueFile* _kvf;
  int _ntot, _n, _root;
  vector<real> _rwrk;
  bool _am_root;
};
  
ValuePutter::ValuePutter (const mpi::ArraySegmenter* as, KeyValueFile* kvf)
  : _as(as), _kvf(kvf)
{
  _ntot = _as->GetNtot();
  _n = _as->GetN();
  _root = mpi::Root();
  _am_root = mpi::AmRoot();
  if (_am_root)
    _rwrk.resize(2 * _ntot);
}

void ValuePutter::PutArray (const string& field, const real* v, int factor) {
  if (_am_root && factor > 2) _rwrk.resize(factor * _ntot);
  _as->Gather(v, &_rwrk[0], _root, factor);
  if (_am_root) {
    Matrix<real> m(factor * _ntot, &_rwrk[0]);
    _kvf->AddMatd(field, m);
  }
}

void Model::ToKvf (KeyValueFile* kvf) const {
  ValuePutter vp(&_ctrl.as, kvf);
  vp.PutArray("mu0", _fr.mu0);
  vp.PutArray("a", _fr.a);
  vp.PutArray("b", _fr.b);
  vp.PutArray("d_c", _fr.d_c);
  vp.PutArray("s_normal", _mb.s_normal);
  vp.PutArray("v_init", _ic.v);
  vp.PutArray("chi_init", _ic.chi);
  vp.PutArray("slip_init", _ic.slip);
  if (_ic.dlt) vp.PutArray("dlt_init", _ic.dlt);
  if (_ic.p) vp.PutArray("p_init", _ic.p);
  if (_ic.T) vp.PutArray("T_init", _ic.T);
}

class ModelValueSetter : public ValueSetter {
public:
  ModelValueSetter(const mpi::ArraySegmenter* as, const KeyValueFile* kvf,
                   list<real*>* ptrs);

  bool SetArray(const string& field, real*& v, int factor = 1);
  void ZeroArray(real*& v, int factor = 1);

private:
  list<real*>* _ptrs;
};

ModelValueSetter::
ModelValueSetter  (const mpi::ArraySegmenter* as, const KeyValueFile* kvf,
                   list<real*>* ptrs)
  : ValueSetter(as, kvf), _ptrs(ptrs)
{}

bool ModelValueSetter::SetArray (const string& field, real*& v, int factor) {
  bool ret = ValueSetter::SetArray(field, v, factor);
  if (v) _ptrs->push_back(v);
  return ret;
}
  
void ModelValueSetter::ZeroArray (real*& v, int factor) {
  ValueSetter::ZeroArray(v, factor);
  _ptrs->push_back(v);
}
  
namespace mdl {
  template<typename T>
  void SetVectorAllPid (const KeyValueFile* kvf, const string& field,
                        vector<T>& v) {
    const Matd* mv;
    if (!mpi_IsTrue(kvf && kvf->GetMatd(field, mv))) return;
    int sz;
    if (mpi::AmRoot()) {
      sz = mv->Size();
      v.resize(sz);
      const double* pmv = mv->GetPtr();
      for (int i = 0, n = v.size(); i < n; i++) v[i] = (T) pmv[i];
    }
    mpi::Bcast(&sz, 1);
    if (!mpi::AmRoot()) v.resize(sz);
    mpi::Bcast(&v[0], sz);
  }

  void TrimLixWts (const mpi::ArraySegmenter* as, const vector<real>& vo,
                   vector<real>& v) {
    int nx = vo.size();
    if (nx == 0 || as->GetNtot() % nx != 0) return;
    int nz = as->GetNtot() / nx;
    // Index of first relevant column.
    int xi = as->GetOffset() / nz;
    // Index of last relevant column.
    int ospn = as->GetOffset() + as->GetN();
    int xf = ospn / nz;
    if (xf * nz == ospn) xf--;
    v.resize(xf - xi + 1);
    for (int i = 0, ix = xi; ix <= xf; i++, ix++) v[i] = vo[ix];
  }

  void TrimSampleIdxs (const mpi::ArraySegmenter* as, const vector<int>& vo,
                       vector<int>& v) {
    int bds[2];
    as->GetIndexBounds(bds);
    for (size_t i = 0; i < vo.size(); i++)
      if (vo[i] >= bds[0] && vo[i] < bds[1]) v.push_back(vo[i]);
  }
}
 
Model* BuildModelFromKeyValueFile (const KeyValueFile* kvf, bool disp) {
  static const char* fn_name = "BuildModelFromKeyValueFile";
  bool am_root = mpi::AmRoot();
  Model* m = new Model();

  ModelValueSetter vs(&m->_ctrl.as, kvf, &m->_ptrs);

  // as is not actually valid yet. But we can call SetScalar safely.
  int nelem_tot, ncomp;
  if (!(vs.SetScalar("nelem", nelem_tot) && vs.SetScalar("ncomp", ncomp))) {
    delete m;
    if (disp && am_root)
      printf("%s: nelem or ncomp not specified.\n", fn_name);
    return NULL;
  }

  m->_ctrl.as.ApportionN(nelem_tot);

  m->_mesh.nelem = m->_ctrl.as.GetN();
  m->_mesh.ncomp = ncomp;
  m->_mesh.n = m->_mesh.nelem * m->_mesh.ncomp;
    
  string evo;
  vs.SetString("evolution", evo);
  if (evo == string("aging")) m->_fr.law = Model::el_Aging;
  else if (evo == string("slip")) m->_fr.law = Model::el_Slip;
  else m->_fr.law = Model::el_Invalid;
  vs.SetArray("mu0", m->_fr.mu0);
  vs.SetArray("a", m->_fr.a);
  vs.SetArray("b", m->_fr.b);
  vs.SetArray("d_c", m->_fr.d_c);
  vs.SetScalar("v0", m->_fr.v0, (real) -1.0);
  int use_vcutoff;
  vs.SetScalar("use_vcutoff", use_vcutoff, 0);
  m->_fr.use_vcutoff = (bool) use_vcutoff;
  if (m->_fr.use_vcutoff) {
    vs.SetArray("vc_v1", m->_fr.v1);
    vs.SetArray("vc_v2", m->_fr.v2);
  }
  m->_fr.a_of_T = false;

  string stress_fn;
  vs.SetString("stress_fn", stress_fn);
  if (stress_fn == string("ss"))
    m->_mb.stress_fn = NewSpringSliderStressFn(vs);
  else if (stress_fn == string("h-matrix"))
    m->_mb.stress_fn = NULL;
  else m->_mb.stress_fn = NULL;

  vs.SetScalar("eta", m->_mb.eta, (real) -1.0);
  vs.SetArray("s_normal", m->_mb.s_normal);

  vs.SetScalar("ti", m->_mesh.ti, -1.0);
  vs.SetScalar("tf", m->_mesh.tf, -1.0);

  vs.SetScalar("rel_tol", m->_ode.rel_tol, m->_ode.rel_tol);

  vs.SetArray("v_init", m->_ic.v, ncomp);
  vs.SetArray("chi_init", m->_ic.chi);
  if (!vs.SetArray("slip_init", m->_ic.slip, ncomp))
    vs.ZeroArray(m->_ic.slip, ncomp);

  m->_dp.use_dilatancy = false;
  m->_dp.use_only_theta = true;
  m->_dp.use_membrane_diffusion = false;
  m->_dp.use_p = m->_dp.use_T = false;
  m->_ic.dlt = m->_ic.p = m->_ic.T = NULL;
  vs.SetArray("t_f_p", m->_dp.t_f_p, ncomp);
  if (m->_dp.t_f_p) {
    vs.SetArray("dlt_init", m->_ic.dlt);
    if (!vs.SetArray("p_init", m->_ic.p, ncomp)) vs.ZeroArray(m->_ic.p, ncomp);
    if (!vs.SetArray("T_init", m->_ic.T, ncomp)) vs.ZeroArray(m->_ic.T, ncomp);
    // Individual elements are turned off by setting t_f to < 0.
    vs.SetArray("t_f_T", m->_dp.t_f_T, ncomp);
    vs.SetScalar("Lambda", m->_dp.Lambda, (real) 0);

    //todo Make these user settable.
    m->_dp.epsilon = 5e-5;
    m->_dp.beta = 6e-11;
    m->_dp.rho_cp = 2600 * 1100;
    m->_dp.h_c = 1e-3;
    m->_dp.h = 1e-3;

    m->_dp.use_dilatancy = true;
    m->_dp.use_only_theta = false;
    m->_dp.use_membrane_diffusion = true;

    for (int i = 0; i < m->_mesh.nelem; i++) {
      if (m->_dp.t_f_p[i] >= 0) {
        if (m->_dp.use_dilatancy) m->_dp.use_p = true;
        if (m->_dp.t_f_T[i] >= 0 && m->_dp.Lambda > 0) {
          m->_dp.use_p = true;
          m->_dp.use_T = true;
          break;
        }
      }
    }
  }

  vs.SetString("stop_indicator", m->_ctrl.stop_indicator);
  if (!m->_ctrl.stop_indicator.empty())
    vs.SetScalar("stop_check_frequency", m->_ctrl.stop_check_frequency, 1);

  vs.SetScalar("disp_every", m->_io.disp_every, -1);
  int ow;
  vs.SetScalar("allow_overwrite", ow, 0);
  m->_io.allow_overwrite = (bool) ow;
    
  // Get stuff for listeners, but don't actually construct them.
  vs.SetString("save_filename", m->_io.save_filename);
  vs.SetScalar("save_v_every", m->_io.v_every, -1);
  vs.SetScalar("save_slip_every", m->_io.slip_every, -1);
  vs.SetScalar("save_state_every", m->_io.state_every, -1);

  vs.SetString("lineint_save_filename", m->_io.lineint_save_filename);
  vs.SetScalar("lineint_save_every", m->_io.lineint_save_every, 1);
  vector<real> w;
  mdl::SetVectorAllPid(kvf, "lineint_wts", w);
  m->_io.lineint_nx = w.size();
  mdl::TrimLixWts(&m->_ctrl.as, w, m->_io.lineint_wts);

  vs.SetScalar("sample_save_every", m->_io.sample_save_every, -1);
  vector<int> idxs;
  mdl::SetVectorAllPid(kvf, "sample_idxs", idxs);
  m->_io.sample_tot_size = idxs.size();
  mdl::TrimSampleIdxs(&m->_ctrl.as, idxs, m->_io.sample_idxs);

  return m;
}

void DeleteModel (Model* m) { delete m; }

void Go (const Model* m) { RunOde(m); }

}
