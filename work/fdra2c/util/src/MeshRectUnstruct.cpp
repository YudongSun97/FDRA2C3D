#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <queue>
#include "util/include/Util.hpp"
#include "util/include/MeshRectUnstruct.hpp"

namespace util {
  namespace rmesh {

    Rect::Rect(FILE* fid) { Deserialize(fid); }
    
    bool Rect::operator<(const Rect& r) const
    {
      if (dx < r.dx) return true;
      if (dx == r.dx) return dy < r.dy;
      return false;
    }

    void Rect::Serialize(FILE* fid) const { write(&x, 4, fid); }

    void Rect::Deserialize(FILE* fid) { read(&x, 4, fid); }

    RectOpts::RectOpts(FILE* fid) { Deserialize(fid); }

    void RectOpts::Serialize(FILE* fid) const { write(&min_len, 2, fid); }

    void RectOpts::Deserialize(FILE* fid) { read(&min_len, 2, fid); }
    
    struct Point {
      double x, y;
    };

    class QuadTree {
    private:
      class Node {
      private:
        Rect _r;
        Node* _kids[4]; // Ordered by quadrant number
        size_t _id;     // Used only in leaves

      public:
        Node(const Rect& r);
        Node(FILE* fid);
        ~Node();

        size_t Id() const { return _id; }
        void Split(const RectOpts& ro, ResolutionFn* rf, size_t* starting_id);
        void PushBackRects(vector<Rect>* rs);
        void Serialize(FILE* fid);
        void Deserialize(FILE* fid);
      };

    private:
      Node* _domain;

    public:
      QuadTree();
      QuadTree(const Rect& domain, const RectOpts& ro, ResolutionFn* rf,
               size_t* starting_id);
      QuadTree(FILE* fid);
      ~QuadTree();
      void PushBackRects(vector<Rect>* rs);
      void Serialize(FILE* fid);

    private:
      void Deserialize(FILE* fid);
    };

    QuadTree::QuadTree() : _domain(NULL) {}

    QuadTree::QuadTree(const Rect& domain, const RectOpts& ro, ResolutionFn* rf,
                       size_t* starting_id)
    {
      _domain = new Node(domain);
      _domain->Split(ro, rf, starting_id);
    }

    QuadTree::QuadTree(FILE* fid) { Deserialize(fid); }

    QuadTree::~QuadTree() { if (_domain) delete _domain; }

    QuadTree::Node::Node(const Rect& r)
      : _r(r), _id(0)
    { for (size_t i = 0; i < 4; i++) _kids[i] = NULL; }
    
    QuadTree::Node::Node(FILE* fid)
    {
      for (size_t i = 0; i < 4; i++) _kids[i] = NULL;
      Deserialize(fid);
    }

    QuadTree::Node::~Node()
    { for (size_t i = 0; i < 4; i++) if (_kids[i]) delete _kids[i]; }

    void QuadTree::Node::Serialize(FILE* fid)
    {
      _r.Serialize(fid);
      char have_kids = (char) (_kids[0] != NULL);
      write(&have_kids, 1, fid);
      if (have_kids)
        for (size_t i = 0; i < 4; i++)
          _kids[i]->Serialize(fid);
      else
        write(&_id, 1, fid);
    }

    void QuadTree::Node::Deserialize(FILE* fid)
    {
      _r.Deserialize(fid);
      char have_kids;
      read(&have_kids, 1, fid);
      if (have_kids) {
        for (size_t i = 0; i < 4; i++) _kids[i] = new Node(fid);
        _id = 0;
      } else
        read(&_id, 1, fid);
    }

    inline void QuadTree::PushBackRects(vector<Rect>* rs)
    { _domain->PushBackRects(rs); }

    void QuadTree::Serialize(FILE* fid) { _domain->Serialize(fid); }

    void QuadTree::Deserialize(FILE* fid) { _domain = new Node(fid); }

    static inline Point& GetRectCenter(const Rect& r, Point& p)
    {
      p.x = r.x + 0.5*r.dx;
      p.y = r.y + 0.5*r.dy;
      return p;
    }

    struct RfCallData {
      ResolutionFn* rf;
      vector<double> x, y, f;
      RfCallData(ResolutionFn* irf, size_t n) : rf(irf), x(1), y(1), f(1) {}
      void Resize(size_t n) { x.resize(n); y.resize(n); f.resize(n); }
    };

    static inline double CallOnSinglePoint(RfCallData& rd, const Point& p)
    {
      rd.Resize(1);
      rd.x[0] = p.x;
      rd.y[0] = p.y;
      rd.rf->Call(rd.x, rd.y, rd.f);
      return rd.f[0];
    }

    static inline double CallOnCenterPoint(RfCallData& rd, const Rect& r)
    {
      Point p;
      GetRectCenter(r, p);
      return CallOnSinglePoint(rd, p);
    }

    static inline bool
    TestShouldSplit(const Rect& r, const RectOpts& ro, double f)
    { return max(r.dx, r.dy) > max(f, ro.min_len); }

    // Breadth-first search on points inside r. The purpose of BFS is to test
    // points that are spread out rather than (by DFS) focusing on an
    // increasingly small region before proceeding to the next.
    static bool
    ShouldSplit(const Rect& r, const RectOpts& ro, RfCallData& rd)
    {
      priority_queue<Rect> pq;
      pq.push(r);
      while (!pq.empty()) {
        Rect cr = pq.top();
        pq.pop();
        if (TestShouldSplit(r, ro, CallOnCenterPoint(rd, cr))) return true;
        double dxh = 0.5*cr.dx, dyh = 0.5*cr.dy;
        if (max(dxh, dyh) > 0.5*ro.min_len) {
          pq.push(Rect(cr.x + dxh, cr.y + dyh, dxh, dyh));
          pq.push(Rect(cr.x      , cr.y + dyh, dxh, dyh));
          pq.push(Rect(cr.x      , cr.y      , dxh, dyh));
          pq.push(Rect(cr.x + dxh, cr.y      , dxh, dyh));
        }
      }
      return false;
    }

    void QuadTree::Node::Split(const RectOpts& ro, ResolutionFn* rf,
                               size_t* id)
    {
      RfCallData rd(rf, 1);
      if (ShouldSplit(_r, ro, rd)) {
        double dxh = 0.5*_r.dx, dyh = 0.5*_r.dy;
        _kids[0] = new Node(Rect(_r.x + dxh, _r.y + dyh, dxh, dyh));
        _kids[1] = new Node(Rect(_r.x      , _r.y + dyh, dxh, dyh));
        _kids[2] = new Node(Rect(_r.x      , _r.y      , dxh, dyh));
        _kids[3] = new Node(Rect(_r.x + dxh, _r.y      , dxh, dyh));
        for (size_t i = 0; i < 4; i++) _kids[i]->Split(ro, rf, id);
      } else {
        _id = *id;
        (*id)++;
      }
    }

    void QuadTree::Node::PushBackRects(vector<Rect>* rs)
    {
      if (_kids[0])
        for (size_t i = 0; i < 4; i++) _kids[i]->PushBackRects(rs);
      else
        rs->push_back(_r);
    }

    // i is a base-1 index.
    static inline void
    GetSurroundingPoints(const Matrix<double>& v, double x, size_t* i)
    {
      const double* pl = lower_bound(v.GetPtr(), v.GetPtr() + v.Size(), x);
      *i = (size_t) (pl - v.GetPtr());
    }

    void LinInterpRF::Call(const vector<double>& x, const vector<double>& y,
                           vector<double>& f)
    {
      for (size_t i = 0, n = x.size(); i < n; i++) {
        double alpha, f1, f2;
        size_t ix, iy;
        GetSurroundingPoints(_x, x[i], &ix);
        GetSurroundingPoints(_y, y[i], &iy);
        alpha = (x[i] - _x(ix)) / (_x(ix + 1) - _x(ix));
        f1 = (1.0 - alpha) * _f(ix, iy    ) + alpha * _f(ix + 1, iy    );
        f2 = (1.0 - alpha) * _f(ix, iy + 1) + alpha * _f(ix + 1, iy + 1);
        alpha = (y[i] - _y(iy)) / (_y(iy + 1) - _y(iy));
        f[i] = (1.0 - alpha) * f1 + alpha * f2;
      }
    }

    static inline void BreakUpDomain(double Dx, double max_len,
                                     double* dx, unsigned int* nx)
    {
      *nx = (unsigned int) ceil(Dx / max_len);
      *dx = Dx / *nx;
    }

    RectMeshUnstruct::
    RectMeshUnstruct(const Rect& domain, const RectOpts& ro, ResolutionFn* rf)
      : _domain(domain), _ro(ro)
    {
      // First break the domain into squares.
      double dx, dy;
      unsigned int nx, ny;
      BreakUpDomain(_domain.dx, ro.max_len, &dx, &nx);
      BreakUpDomain(_domain.dy, ro.max_len, &dy, &ny);
      // Now split each square.
      _qts.reserve(ny * nx);
      size_t id = 0;
      for (size_t iy = 0, k = 0; iy < ny; iy++)
        for (size_t ix = 0; ix < nx; ix++, k++) {
          QuadTree* qt = new QuadTree
            (Rect(_domain.x + ix*dx, _domain.y + iy*dy, dx, dy), ro, rf, &id);
          _qts.push_back(qt);
        }
    }

    RectMeshUnstruct::RectMeshUnstruct(const string& filename)
      throw (FileException)
      : _domain(0.0, 0.0, 0.0, 0.0), _ro(0.0, 0.0)
    {
      FILE* fid = fopen(filename.c_str(), "r");
      if (!fid) throw FileException("Can't read " + filename);
      Deserialize(fid);
      fclose(fid);
    }

    RectMeshUnstruct::~RectMeshUnstruct()
    { for (size_t i = 0; i < _qts.size(); i++) delete _qts[i]; }

    void RectMeshUnstruct::Serialize(const string& filename) const
      throw (FileException)
    {
      FILE* fid = fopen(filename.c_str(), "w");
      if (!fid) throw FileException("Can't read " + filename);
      Serialize(fid);
      fclose(fid);
    }
    
    void RectMeshUnstruct::Serialize(FILE* fid) const throw (FileException)
    {
      _domain.Serialize(fid);
      _ro.Serialize(fid);
      long long int n = _qts.size();
      write(&n, 1, fid);
      for (size_t i = 0; i < _qts.size(); i++) _qts[i]->Serialize(fid);
    }

    void RectMeshUnstruct::Deserialize(FILE* fid)
    {
      try {
        _domain.Deserialize(fid);
        _ro.Deserialize(fid);
        long long int n;
        read(&n, 1, fid);
        _qts.reserve(n);
        for (size_t i = 0; i < (size_t) n; i++) {
          QuadTree* qt = new QuadTree(fid);
          _qts.push_back(qt);
        }
      } catch (const FileException& e) {
        for (size_t i = 0; i < _qts.size(); i++) delete _qts[i];
        throw e;
      }
    }

    const vector<Rect>& RectMeshUnstruct::GetRects() const
    {
      if (_rs.empty()) // Lazy evaluation
        for (size_t i = 0; i < _qts.size(); i++)
          _qts[i]->PushBackRects(&_rs);
      return _rs;
    }

    int GetRectId(double x, double y) const
    {
      //todo
    }

  }
}
