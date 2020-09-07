#ifndef INCLUDE_UTIL_MESHRECTUNSTRUCT
#define INCLUDE_UTIL_MESHRECTUNSTRUCT

#include <vector>
#include "util/include/Matrix.hpp"
#include "util/include/Exception.hpp"

namespace util {
  namespace rmesh {
    using namespace std;
    
    struct Rect {
      double x, y;   // lower-left corner
      double dx, dy; // x,y-direction lengths

      Rect() : x(0.0), y(0.0), dx(0.0), dy(0.0) {}
      Rect(FILE* fid);
      Rect(double ix, double iy, double idx, double idy)
        : x(ix), y(iy), dx(idx), dy(idy) {}
      bool operator<(const Rect& r) const;
      void Serialize(FILE* fid) const;
      void Deserialize(FILE* fid);
    };

    struct RectOpts {
      double min_len, max_len;

      RectOpts() : min_len(0.0), max_len(0.0) {}
      RectOpts(FILE* fid);
      RectOpts(double imin_len, double imax_len)
        : min_len(imin_len), max_len(imax_len) {}
      void Serialize(FILE* fid) const;
      void Deserialize(FILE* fid);
    };

    class ResolutionFn {
    public:
      virtual ~ResolutionFn() {}
      // Calculate f[i] = f(x[i], y[i]) for each i.
      virtual void Call(const vector<double>& x, const vector<double>& y,
                        vector<double>& f) = 0;
    };

    class LinInterpRF : public ResolutionFn {
    public:
      LinInterpRF(const Matrix<double>& x, const Matrix<double>& y,
                  const Matrix<double>& f)
        : _x(x), _y(y), _f(f)
      {
        assert(_y.Size() == _f.Size(1));
        assert(_x.Size() == _f.Size(2));
      }

      virtual void Call(const vector<double>& x, const vector<double>& y,
                        vector<double>& f);

    private:
      Matrix<double> _x, _y, _f;
    };

    class QuadTree;

    class RectMeshUnstruct {
    public:
      RectMeshUnstruct(const Rect& domain, const RectOpts& ro,
                       ResolutionFn* rf);
      RectMeshUnstruct(const string& filename) throw (FileException);
      ~RectMeshUnstruct();

      void Serialize(const string& filename) const throw (FileException);
      void Serialize(FILE* fid) const throw (FileException);

      const Rect& GetDomain() const { return _domain; }
      const RectOpts& GetRectOpts() const { return _ro; }
      const vector<Rect>& GetRects() const;

      // Return index i such that GetRects()[i] is the rectangle that contains
      // the point (x, y). Result is -1 if (x, y) is outside of the domain.
      int GetRectId(double x, double y) const;

    private:
      Rect _domain;
      RectOpts _ro;
      vector<QuadTree*> _qts;
      mutable vector<Rect> _rs;

    private:
      void Deserialize(FILE* fid);
    };

  }
}

#endif
