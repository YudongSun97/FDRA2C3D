#include "Elastostatics.hpp"

class OkadaRectGf : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd();
  virtual bool Call(const vector<uint>& rs, const vector<uint>& cs,
                    double* B);
  virtual void DoExtraTasksSerial() throw (Exception);
  
private:
  string _bc_filename;
  size_t _component;
  vector<hmmvp::es::dc3::Elem> _es, _bc_es;
  hmmvp::es::LameParms _lp;
  double _disl[3];
};

void OkadaRectGf::DoExtraTasksSerial() throw (Exception)
{
  // Compute the BC.
  Matrix<double> B(_es.size());
  double* pB = B.GetPtr();
  for (size_t ie = 0; ie < _es.size(); ie++) {
    pB[ie] = 0.0;
    for (size_t ib = 0; ib < _bc_es.size(); ib++)
      pB[ie] += hmmvp::es::GetTractionComp
        (_lp, _bc_es[ib], _disl, _es[ie], _component);
  }

  FILE* fid = fopen(_bc_filename.c_str(), "w");
  if (!fid) throw Exception("Can't write BC file.");
  write(pB, _es.size(), fid);
  fclose(fid);
}

Hd* OkadaRectGf::ComputeHd()
{
  size_t nes = _es.size();
  Matrix<double> ctrs(3, nes);
  double* pc = ctrs.GetPtr();
  for (size_t i = 0; i < nes; i++) {
    memcpy(pc, _es[i].Center(), 3*sizeof(double));
    pc += 3;
  }
  Hd* hd = NewHd(ctrs);
  errpr("nbr blocks = %ld\n", hd->NbrBlocks());
  return hd;
}

void OkadaRectGf::Init(const KeyValueFile* kvf) throw (Exception)
{
  using namespace hmmvp::es;
  using namespace dc3;
  const string* s;
  const Matrix<double>* m;
  double d;

  Matrix<double> x, eta;
  if (!kvf->GetMatd("x", m)) throw Exception("Missing x.");
  x = *m;
  if (!kvf->GetMatd("eta", m)) throw Exception("Missing eta.");
  eta = *m;

  if (!kvf->GetString("bc_filename", s))
    throw Exception("Missing bc_filename.");
  _bc_filename = *s;

#define get(a) if (!kvf->GetDouble(#a, a)) throw Exception("Missing " #a);
  double mu, nu, dipdeg, y_min, depth_min, disl_dip, disl_strike, disl_tensile;
  get(mu); get(nu);
  get(dipdeg); get(y_min); get(depth_min);
  get(disl_dip); get(disl_strike); get(disl_tensile);
#undef get

  _component = 1;
  if (kvf->GetDouble("component", d)) _component = (size_t) d;
  if (_component > 2) throw Exception("component must be 0, 2, or 2");

  _lp.Set(mu, nu);
  _disl[0] = disl_strike; _disl[1] = disl_dip; _disl[2] = disl_tensile;

  // Fault elements.
  PlanarTensorMeshToElems(x, eta, y_min, depth_min, dipdeg, _es);

  // Get these elements for the boundary conditions. The BCs are implemented by
  // huge slabs adjacent to the simulated fault.
  Matrix<double> xlim(2), etalim(2);
  vector<Elem> bc_es;
  xlim(1) = x(1); xlim(2) = x(x.Size());
  etalim(1) = eta(1); etalim(2) = eta(eta.Size());
  double L = 100.0*std::max(xlim(2) - xlim(1), etalim(2) - etalim(1));
  // -x
  xlim(1) = x(1) - L;
  xlim(2) = x(1);
  PlanarTensorMeshToElems(xlim, etalim, y_min, depth_min, dipdeg, bc_es);
  // +x
  xlim(1) = x(x.Size());
  xlim(2) = xlim(1) + L;
  PlanarTensorMeshToElems(xlim, etalim, y_min, depth_min, dipdeg, bc_es);
  // -y or y
  xlim(1) = x(1) - L;
  xlim(2) = x(x.Size()) + L;
  double dm, ym;
  if (dipdeg >= 0.0) {
    etalim(1) = eta(1) - L;
    etalim(2) = eta(2);
    dm = depth_min + (eta(eta.Size()) - eta(1))*sind(dipdeg);
    ym = y_min - L*cosd(dipdeg);
  } else {
    etalim(1) = eta(eta.Size());
    etalim(2) = etalim(1) + L;
    dm = depth_min - (eta(eta.Size()) - eta(1))*sind(dipdeg);
    ym = y_min + (eta(eta.Size()) - eta(1))*cosd(dipdeg);
  }
  PlanarTensorMeshToElems(xlim, etalim, ym, dm, dipdeg, bc_es);
  _bc_es.insert(_bc_es.end(), bc_es.begin(), bc_es.end());
}

bool OkadaRectGf::
Call(const vector<uint>& rs, const vector<uint>& cs, double* B)
{
  for (size_t ic = 0, k = 0; ic < cs.size(); ic++)
    for (size_t ir = 0; ir < rs.size(); ir++, k++)
      B[k] = hmmvp::es::GetTractionComp
        (_lp, _es[cs[ic] - 1], _disl, _es[rs[ir] - 1], _component);
  return true;
}
