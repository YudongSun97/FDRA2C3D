#include "Elastostatics.hpp"

class OssrpGf : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd();
  virtual bool Call(const vector<uint>& rs, const vector<uint>& cs,
                    double* B);
  virtual void DoExtraTasksSerial() throw (Exception);
  
private:
  string _bc_filename;
  size_t _component, _nx, _neta;
  vector<hmmvp::es::dc3::Elem> _es, _bc_es;
  hmmvp::es::LameParms _lp;
  double _disl[3];
  size_t _nrepeat; // Number of repeated images in each of -/+x directions
  double _L;       // Along-strike length
};

void OssrpGf::Init(const KeyValueFile* kvf) throw (Exception)
{
  using namespace hmmvp::es;
  using namespace dc3;
  const string* s;
  const Matrix<double>* m;

  Matrix<double> x, eta;
  if (!kvf->GetMatd("x", m)) throw Exception("Missing x.");
  x = *m;
  if (!kvf->GetMatd("eta", m)) throw Exception("Missing eta.");
  eta = *m;

  _nx = 2*(x.Size() - 1);
  _neta = eta.Size() - 1;
  _L = 2*x(x.Size());

  if (!kvf->GetString("bc_filename", s))
    throw Exception("Missing bc_filename.");
  _bc_filename = *s;

#define get(a) if (!kvf->GetDouble(#a, a)) throw Exception("Missing " #a);
  double mu, nu, dipdeg, y_min, depth_min, nrepeat;
  get(mu); get(nu);
  get(dipdeg); get(y_min); get(depth_min); get(nrepeat);
#undef get

  if (nrepeat < 0.0 || nrepeat > 100.0)
    throw Exception("nrepeat must be between 0 and 100.");
  _nrepeat = (size_t) nrepeat;
  _component = 1;
  _disl[0] = 0.0; _disl[1] = 1.0; _disl[2] = 0.0;
  _lp.Set(mu, nu);

  // Fault elements.
  PlanarTensorMeshToElems(x, eta, y_min, depth_min, dipdeg, _es);
  // Duplicate with -x for the other half of the fault.
  for (size_t i = 0, n = _es.size(); i < n; i++) {
    // Specify using the center so that we can flip the sign of the x component.
    const Elem& e = _es[i];
    double dx = 0.5*(e.al1() + e.al2());
    double de = 0.5*(e.aw1() + e.aw2());
    _es.push_back(Elem(-e.Center()[2], e.dipdeg(), dx, dx, de, de,
                       -e.Center()[0], e.Center()[1]));
  }

  // Get these elements for the boundary conditions. The BCs are implemented by
  // huge slabs adjacent to the simulated fault.
  Matrix<double> xlim(2), etalim(2);
  double L = 10.0*(1 + _nrepeat)*
    std::max(2*x(x.Size()), eta(eta.Size()) - eta(1));
  // -y or y
  xlim(2) = x(x.Size()) + L;
  xlim(1) = -xlim(2);
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
  PlanarTensorMeshToElems(xlim, etalim, ym, dm, dipdeg, _bc_es);
}

void OssrpGf::DoExtraTasksSerial() throw (Exception)
{
  // Compute the BC on the range elements.
  size_t n = _es.size() / 2;
  Matrix<double> B(n);
  double* pB = B.GetPtr();
  for (size_t ie = 0; ie < n; ie++) {
    pB[ie] = 0.0;
    for (size_t ib = 0; ib < _bc_es.size(); ib++)
      pB[ie] += hmmvp::es::GetTractionComp
        (_lp, _bc_es[ib], _disl, _es[ie], _component);
  }

  FILE* fid = fopen(_bc_filename.c_str(), "w");
  if (!fid) throw Exception("Can't write BC file.");
  write(pB, n, fid);
  fclose(fid);
}

Hd* OssrpGf::ComputeHd()
{
  size_t nes = _es.size();

  Matrix<double> D(3, nes);   // Domain is the full fault.
  Matrix<double> R(3, nes/2); // Range is the +x side of the fault.
  double* pD = D.GetPtr();
  double* pR = R.GetPtr();
  for (size_t i = 0; i < nes; i++) {
    memcpy(pD, _es[i].Center(), 3*sizeof(double));
    pD += 3;
    if (i < nes/2) {
      memcpy(pR, _es[i].Center(), 3*sizeof(double));
      pR += 3;
    }
  }

  Hd* hd = NewHd(D, R);
  errpr("nbr blocks = %ld\n", hd->NbrBlocks());
  return hd;
}

bool OssrpGf::
Call(const vector<uint>& rs, const vector<uint>& cs, double* B)
{
  using namespace hmmvp::es;
  using namespace dc3;
  for (size_t ic = 0, k = 0; ic < cs.size(); ic++)
    for (size_t ir = 0; ir < rs.size(); ir++, k++) {
      const Elem& es = _es[cs[ic] - 1];
      B[k] = GetTractionComp(_lp, es, _disl, _es[rs[ir] - 1], _component);
      // Periodic images
      for (size_t ip = 1; ip < _nrepeat + 1; ip++) {
        Elem esr1(es.depth(), es.dipdeg(), es.al1(), es.al2(), es.aw1(),
                  es.aw2(), es.gx() + ip*_L, es.gy());
        B[k] += GetTractionComp(_lp, esr1, _disl, _es[rs[ir] - 1], _component);
        Elem esr2(es.depth(), es.dipdeg(), es.al1(), es.al2(), es.aw1(),
                  es.aw2(), es.gx() - ip*_L, es.gy());
        B[k] += GetTractionComp(_lp, esr2, _disl, _es[rs[ir] - 1], _component);
      }
    }
  return true;
}
