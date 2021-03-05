#include <iostream>
#include <algorithm>
#include "util/include/CodeAnalysis.hpp"
#include "util/include/Elastostatics.hpp"
#include "util/src/BrickMeshBemBuilder_pri.hpp"

namespace bm {
using namespace std;
using namespace util;
using namespace rmesh;

class OkadaBrickMeshGf : public ImplGreensFn, public rmesh::GreensFn {
public:
  OkadaBrickMeshGf();
  virtual ~OkadaBrickMeshGf();
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd();
  virtual void AmEstimatingBfro() { _bmb->UsePoorQuality(); }
  virtual void AmFormingHmat() { _bmb->UseFullQuality(); }
  virtual bool Call(const vector<UInt>& rs, const vector<UInt>& cs, double* B);
  virtual double Call(const Rect& r, double x, double y) const;
  double Call(Dir::Enum bdy, double rx, double ry) const;
  virtual void DoExtraTasksMpi() throw (Exception);

protected:
  string _bc_filename;
  size_t _component;
  es::LameParms _lp;
  double _disl[3];

  RectMeshUnstruct* _rmu;
  RmuAnalyzer* _ra;
  MeshAnalyzer* _ma;
  LinearInterpolatorMatrix* _lim;
  BrickMeshBemBuilder* _bmb;

  struct {
    // Keep some data from the input.
    double y_min, depth_min, dipdeg;
    size_t n_per_layers;
    double dist_fac;
    // And some derived data.
    double depth0, depth1; // GetDomain().y and .y + .dy
    double cd, sd;
  } _in;

private:
  void ConstructBc(const KeyValueFile* kvf, Boundaries& b) throw (Exception);
  void PrintState();
};

OkadaBrickMeshGf::OkadaBrickMeshGf()
  : _rmu(NULL), _ra(NULL), _ma(NULL), _lim(NULL), _bmb(NULL)
{
  Ca::GetTimer()->Reset(40); Ca::GetTimer()->Reset(41);
}

OkadaBrickMeshGf::~OkadaBrickMeshGf()
{
  if (_bmb) {
    size_t mns[2];
    mns[0] = _bmb->GetCtrNcalls();
    mns[1] = _bmb->GetCtrNgfs();
    size_t ns[2];
    mpi::Reduce(mns, ns, 2, MPI_SUM);
    if (mpi::AmRoot()) {
      size_t max_nc = 0;
      for (size_t i = 0; i < _ma->Nc().size(); i++)
        max_nc = std::max(max_nc, _ma->Nc()[i]);
      const Boundaries& b = _ma->GetBoundaries();
      size_t pf = 0;
      if (b.GetBC(Dir::E) == Boundaries::bc_periodic) pf++;
      if (b.GetBC(Dir::N) == Boundaries::bc_periodic) pf++;
      switch (pf) {
      case 0: pf = 1; break;
      case 1: pf = 2 * _in.n_per_layers + 1; break;
      case 2: pf = 2 * _in.n_per_layers + 1; pf *= pf; break;
      default: assert(false); break;
      }
      size_t exact_tngfs = pf * ns[0] * max_nc*max_nc;
      printf("calls: std %lld IGA %lld (%1.2f) exact IGA %lld (%1.2f)\n",
             (long long int) (pf * ns[0]), (long long int) ns[1],
             (double) ns[1] / (pf * ns[0]), (long long int) exact_tngfs,
             (double) exact_tngfs / (pf * ns[0]));
    }
  }

  if (_bmb) DeleteBrickMeshBemBuilder(_bmb);
  if (_lim) delete _lim;
  if (_ma) delete _ma;
  if (_ra) DeleteRmuAnalyzer(_ra);
  if (_rmu) DeleteRectMeshUnstruct(_rmu);

  caprint("OkadaBrickMeshGf: tot call (%1.2e) 1.00 dc3d %1.3f\n",
          Ca::GetTimer()->TotEt(40),
          Ca::GetTimer()->TotEt(41)/Ca::GetTimer()->TotEt(40));
}

void OkadaBrickMeshGf::Init(const KeyValueFile* kvf) throw (Exception)
{
  using namespace util::es;
  using namespace dc3;
  const string* s;
  double d;

  if (!kvf->GetString("rmesh_filename", s))
    throw Exception("Missing rmesh_filename.");
  string rmesh_filename = *s + ".rmu";
    
  if (!kvf->GetString("bc_filename", s))
    throw Exception("Missing bc_filename.");
  _bc_filename = *s;

#define get(a) if (!kvf->GetDouble(#a, a)) throw Exception("Missing " #a);
  double mu, nu, dipdeg, y_min, depth_min, disl_dip, disl_strike,
    disl_tensile;
  get(mu); get(nu);
  get(dipdeg); get(depth_min); y_min = 0.0;
  get(disl_dip); get(disl_strike); get(disl_tensile);
#undef get
  if (depth_min < 0) throw Exception("depth_min must be >= 0.");
  _in.y_min = y_min;
  _in.depth_min = depth_min;
  _in.dipdeg = dipdeg;
  _in.cd = cosd(dipdeg);
  _in.sd = sind(dipdeg);

  _in.dist_fac = 4.0;
  if (kvf->GetDouble("neighborhood", d)) {
    if (d < 0.0) throw Exception("neighborhood must be >= 0");
    _in.dist_fac = d;
  }

  _in.n_per_layers = 1;
  if (kvf->GetDouble("bc_periodic_nlayers", d)) {
    if (d < 0 || (size_t) d > 5)
      throw Exception("Choose bc_periodic_nlayers to be between 0 and 5.");
    _in.n_per_layers = (size_t) d;
  }
    
  _component = 1;
  if (kvf->GetDouble("component", d)) {
    _component = (size_t) d;
    if (_component > 2) throw Exception("component must be 0, 1, or 2");
  }

  _lp.Set(mu, nu);
  _disl[0] = disl_strike; _disl[1] = disl_dip; _disl[2] = disl_tensile;

  _rmu = NewRectMeshUnstruct(rmesh_filename);
  Boundaries b;
  ConstructBc(kvf, b);
  _ra = NewRmuAnalyzer(_rmu, b);
  _ma = new MeshAnalyzer(_rmu, _ra);
  _lim = new LinearInterpolatorMatrix(_ma);
  _bmb = NewBrickMeshBemBuilder(this, _ma, _lim, _in.dist_fac,
                                _in.n_per_layers);

  _in.depth0 = _in.depth_min + fabs(_in.sd) * _rmu->GetDomain().dy;
  _in.depth1 = depth_min;
  if (_in.sd < 0.0) std::swap(_in.depth0, _in.depth1);

  PrintState();
}

namespace {
  Dir::Enum CharToDir(char d)
  {
    switch (d) {
    case 'e': return Dir::E;
    case 'w': return Dir::W;
    case 'n': return Dir::N;
    case 's': return Dir::S;
    default: return Dir::inside;
    }
  }

  struct Entry {
    double d;
    Dir::Enum dir;
    bool operator<(const Entry& e2) const { return d > e2.d; }
  };
}

void OkadaBrickMeshGf::
ConstructBc(const KeyValueFile* kvf, Boundaries& b)
  throw (Exception)
{
  //todo
  // - generate an error if _in.n_per_layers is too large for the depth and N
  //   bdy is periodic.
  // - free BC.
  // - for PID 1, output the setup.

  const char* bc_strs[] = {"e0vbc", "n0vbc", "w0vbc", "s0vbc",
                           "evbc", "nvbc", "wvbc", "svbc",
                           "ewpbc", "nspbc", "wepbc", "snpbc"};
  const size_t n_bc_strs = sizeof(bc_strs) / sizeof(char*);
  vector<Entry> dominance(4);
  for (size_t i = 0; i < 4; i++) dominance[i].dir = (Dir::Enum) i;

  // Process the BC specification.
  for (size_t i = 0; i < n_bc_strs; i++) {
    double d;
    if (kvf->GetDouble(bc_strs[i], d)) {
      if (bc_strs[i][1] == '0' || bc_strs[i][1] == 'v') {
        Dir::Enum dir = CharToDir(bc_strs[i][0]);
        if (bc_strs[i][1] == '0') b.SetZeroVelocityBC(dir);
        else {
          // Rect() is a stand-in for what comes later.
          b.SetVelocityBC(dir, Rect());
        }
        if (d < 0) d = 0;
        dominance[(int) dir].d = d;
        dominance[(int) dir].dir = dir;
      } else {
        Dir::Enum dir = bc_strs[i][0] == 'e' || bc_strs[i][1] == 'e'
          ? Dir::E : Dir::N;
        b.SetPeriodicBC(dir);
        dominance[(int) dir].d = -1;
        dominance[(int) dir].dir = dir;
        dir = Dir::Opposite(dir);
        dominance[(int) dir].d = -1;
        dominance[(int) dir].dir = dir;
      }
    }
  }
  std::sort(dominance.begin(), dominance.end());

  // Set up the v-BC slabs.
  const Rect& dmn = _rmu->GetDomain();
  const double L = std::max(100.0, 10.0*(1 + _in.n_per_layers)) *
    std::max(dmn.dx, dmn.dy);
  // Prevent a v-BC slab from going the the surface.
  Dir::Enum surf_dir = _in.dipdeg >= 0.0 ? Dir::N : Dir::S;
  double L_max = _in.depth_min / fabs(_in.sd);
  if (L_max > L) L_max = L;
  vector<bool> filled(9, false);
  for (size_t i = 0; i < 4; i++) {
    Dir::Enum dir = dominance[i].dir;
    Dir::Enum ccw = Dir::OneCcw(dir), cw = Dir::OneCw(dir);
    if (b.GetBC(dir) == Boundaries::bc_velocity) {
      const bool cw_filled = filled[(int) cw], ccw_filled = filled[(int) ccw];
      Rect r;
      switch (dir) {
      case Dir::E: {
        r.x = dmn.x + dmn.dx;
        r.dx = L;
        r.y = dmn.y - (cw_filled ? 0 : (surf_dir == Dir::S ? L_max : L));
        double y1 = dmn.y + dmn.dy +
          (ccw_filled ? 0 : (surf_dir == Dir::N ? L_max : L));
        r.dy = y1 - r.y;
      } break;
      case Dir::N: {
        r.y = dmn.y + dmn.dy;
        r.dy = (surf_dir == Dir::N ? L_max : L);
        r.x = dmn.x - (ccw_filled ? 0 : L);
        double x1 = dmn.x + dmn.dx + (cw_filled ? 0 : L);
        r.dx = x1 - r.x;
      } break;
      case Dir::W: {
        r.x = dmn.x - L;
        r.dx = L;
        r.y = dmn.y - (ccw_filled ? 0 : (surf_dir == Dir::S ? L_max : L));
        double y1 = dmn.y + dmn.dy +
          (cw_filled ? 0 : (surf_dir == Dir::N ? L_max : L));
        r.dy = y1 - r.y;
      } break;
      case Dir::S: {
        r.y = dmn.y - (surf_dir == Dir::S ? L_max : L);
        r.dy = surf_dir == Dir::S ? L_max : L;
        r.x = dmn.x - (cw_filled ? 0 : L);
        double x1 = dmn.x + dmn.dx + (ccw_filled ? 0 : L);
        r.dx = x1 - r.x;
      } break;
      default: assert(false);
      }
      b.SetVelocityBC(dir, r);
    }
    if (Boundaries::IsVBC(b.GetBC(dir)))
      filled[(int) dir] = filled[(int) cw] = filled[(int) ccw] = true;
  }
  // Check for a free boundary condition.
  if (_in.depth_min == 0) b.SetFreeBC(surf_dir);
}

namespace {
  ostream& operator<<(ostream& os, const Rect& r)
  {
    os << "[x: " << r.x << " to " << r.x + r.dx << ", eta: " << r.y << " to "
       << r.y + r.dy  << "]";
    return os;
  }

  string DirToString(Dir::Enum dir)
  {
    switch (dir) {
    case Dir::E: return " East";
    case Dir::N: return "North";
    case Dir::W: return " West";
    case Dir::S: return "South";
    default: assert(false); return "oops";
    }
  }

  string BCToString(Boundaries::BC bc)
  {
    switch (bc) {
    case Boundaries::bc_velocity:  return "    Velocity";
    case Boundaries::bc_0velocity: return "  0-Velocity";
    case Boundaries::bc_periodic:  return "    Periodic";
    case Boundaries::bc_free:      return "Free Surface";
    }
  }

  string ComponentToString(size_t c)
  {
    switch (c) {
    case 0: return "along-strike";
    case 1: return "along-dip";
    case 2: return "fault-normal";
    default: assert(false); return "oops";
    }
  }
}

void OkadaBrickMeshGf::PrintState()
{
  if (!mpi::AmRoot()) return;
  cout << "--------- start OkadaBrickMeshGf description ---------" << endl;
  cout << "Domain: " << _rmu->GetDomain() << "." << endl
       << "Boundary conditions:" << endl;
  bool have_periodic = false;
  for (size_t i = 0; i < 4; i++) {
    Dir::Enum dir = (Dir::Enum) i;
    Boundaries::BC bc = _ra->GetBoundaries().GetBC(dir);
    if (bc == Boundaries::bc_periodic) have_periodic = true;
    cout << "  " << DirToString(dir) << ": " << BCToString(bc);
    if (bc == Boundaries::bc_velocity)
      cout << " with slab " << _ra->GetBoundaries().GetRect(dir);
    cout << "." << endl;
  }
  if (have_periodic)
    
    cout << "Using " << _in.n_per_layers << " layer"
         << (_in.n_per_layers != 1 ? "s" : "") << " of periodic source images."
         << endl;
  Dir::Enum surf_dir = _in.dipdeg >= 0.0 ? Dir::N : Dir::S;
  cout << "The updip side is " << DirToString(surf_dir) << "; it is at depth "
    << _in.depth_min << "." << endl;
  cout << "The fault dips at " << fabs(_in.dipdeg) << " degrees." << endl;
  cout << "mu is " << _lp.mu() << " and nu is " << _lp.nu() << "." << endl;
  cout << "Computing component: " << ComponentToString(_component) << "."
       << endl;
  cout << "Tractions are: along-strike: " << _disl[0]
       << ", along-dip: " << _disl[1]
       << ", fault-normal: " << _disl[2] << "." << endl;
  cout << "The neighborhood factor is " << _in.dist_fac << "." << endl;
  const size_t ncpu = mpi::GetNproc();
  cout << "Using " << ncpu << " MPI process" << (ncpu > 1 ? "es." : ".")
       << endl;
  cout << "---------- end OkadaBrickMeshGf description ----------" << endl;
}

Hd* OkadaBrickMeshGf::ComputeHd()
{
  // Periodic boundaries for the hierarchical decomposition.
  Matrix<double> pbs(2, 3);
  pbs.Zero();
  const Rect& d = _ma->GetDomain();
  if (_ma->GetBoundaries().GetBC(Dir::E) == Boundaries::bc_periodic) {
    pbs(1, 1) = d.x;
    pbs(2, 1) = d.x + d.dx;
  }
  if (_ma->GetBoundaries().GetBC(Dir::N) == Boundaries::bc_periodic) {
    pbs(1, 2) = d.y;
    pbs(2, 2) = d.y + d.dy;
  }

  // Rect centers in (x, eta) space. Could do it in (x, y, z) space, but it's
  // easier to specifiy the periodic boundaries in (x, eta) space.
  const vector<Rect>& rs = _ma->GetRects();
  const Matrix<double>& xs = _ma->GetX();
  Matrix<double> ctrs(3, rs.size());
  for (size_t i = 1; i <= rs.size(); i++) {
    ctrs(1, i) = xs(1, i);
    ctrs(2, i) = xs(2, i);
    ctrs(3, i) = 0;
  }

  Hd* hd = NewHd(ctrs, &pbs);
  return hd;
}

void OkadaBrickMeshGf::DoExtraTasksMpi() throw (Exception)
{
  const size_t nr = _ma->GetRects().size();
  Matrix<double> B(nr, 4);
  B.Zero();
  Ca::GetTimer()->Tic(40);
  for (size_t j = 0; j < 4; j++)
    for (size_t i = mpi::Pid(); i < nr; i += mpi::GetNproc())
      B(i+1, j+1) = _bmb->Call(i, (Dir::Enum) j);
  Ca::GetTimer()->Toc(40);
  Matrix<double> Bg;
  if (mpi::AmRoot()) Bg.Resize(nr, 4);
  mpi::Reduce(B.GetPtr(), Bg.GetPtr(), 4*nr, MPI_SUM);
    
  if (mpi::AmRoot()) {
    FILE* fid = fopen(_bc_filename.c_str(), "w");
    if (!fid) throw Exception("Can't write BC file.");
    write(Bg.GetPtr(), 4*nr, fid);
    fclose(fid);
  }
}

bool OkadaBrickMeshGf::
Call(const vector<UInt>& rs, const vector<UInt>& cs, double* B)
{
  Ca::GetTimer()->Tic(40);
  for (size_t ic = 0, k = 0; ic < cs.size(); ic++)
    for (size_t ir = 0; ir < rs.size(); ir++, k++)
      B[k] = _bmb->Call(rs[ir] - 1, cs[ic] - 1);
  Ca::GetTimer()->Toc(40);
  return true;
}

double OkadaBrickMeshGf::Call(const Rect& src, double rx, double ry) const
{
  double eta0, alpha;
  eta0 = src.y - _ma->GetDomain().y;
  alpha = eta0 / _ma->GetDomain().dy;
  es::dc3::Elem es
    ((1.0 - alpha) * _in.depth0 + alpha * _in.depth1,
     _in.dipdeg, _in.cd, _in.sd, 0.0, src.dx, 0.0, src.dy, src.x,
     _in.y_min + _in.cd * eta0);
  double s[6];
  { double obs[3], u[3], du[9];
    eta0 = ry - _ma->GetDomain().y;
    alpha = eta0 / _ma->GetDomain().dy;
    obs[0] = rx;
    obs[1] = _in.y_min + eta0 * _in.cd;
    obs[2] = -((1.0 - alpha) * _in.depth0 + alpha * _in.depth1);
    Ca::GetTimer()->Tic(41);
    Dc3d(_lp, es, _disl, obs, u, du);
    Ca::GetTimer()->Toc(41);
    DuToS(_lp, du, s); }
  // Can just uses source's vectors since we're on a plane.
  return es::ProjectStress(s, es, _component);
}
}
