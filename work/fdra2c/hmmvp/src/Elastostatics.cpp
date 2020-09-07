#include <stdio.h>
#include <string.h>
#include "Elastostatics.hpp"

namespace hmmvp {
  namespace es {
    namespace dc3 {

      Elem::Elem(double depth, double dipdeg,
                 double al1, double al2, double aw1, double aw2,
                 double gx, double gy)
        : _depth(depth), _dipdeg(dipdeg),
          _al1(al1), _al2(al2), _aw1(aw1), _aw2(aw2),
          _gx(gx), _gy(gy)
      {
        double cd = cosd(dipdeg), sd = sind(dipdeg);

        _ctr[0] = gx    + 0.5*   (al2 - al1);
        _ctr[1] =  gy    + 0.5*cd*(aw2 - aw1);
        _ctr[2] = -depth + 0.5*sd*(aw2 - aw1);

        _adip[0]    = 0.0; _adip[1]    =  cd;  _adip[2]    = sd;
        _astrike[0] = 1.0; _astrike[1] =  0.0; _astrike[2] = 0.0;
        _normal[0]  = 0.0; _normal[1]  = -sd;  _normal[2]  = cd;
      }

      static bool IsAscending(const Matrix<double>& v)
      {
        const double* pv = v.GetPtr();
        for (int i = 1; i < v.Size(); i++) if (pv[i] <= pv[i-1]) return false;
        return true;
      }

      void PlanarTensorMeshToElems
      (const Matrix<double>& mx, const Matrix<double>& meta, double y_min,
       double depth_min, double dipdeg, vector<Elem>& es)
        throw (Exception)
      {
        if (!IsAscending(mx))
          throw Exception("x must be in ascending order.");
        if (!IsAscending(meta)) {
          throw Exception("eta must be in ascending order.");
        }
        if (depth_min < 0.0) throw Exception("depth_min must be >= 0");

        size_t nx = mx.Size() - 1, ny = meta.Size() - 1;
        const double* x = mx.GetPtr();
        const double* eta = meta.GetPtr();
        es.reserve(es.size() + nx*ny);
        double al1 = 0.0;
        double aw1 = 0.0;
        double sd = sind(dipdeg), cd = cosd(dipdeg);
        double depth_start = (sd >= 0.0) ?
          depth_min + sd*(eta[ny] - eta[0]) :
          depth_min;
        for (size_t iy = 0; iy < ny; iy++) {
          double gy = y_min + cd*(eta[iy] - eta[0]);
          double aw2 = eta[iy+1] - eta[iy];
          double depth = depth_start - sd*(eta[iy] - eta[0]);
          for (size_t ix = 0; ix < nx; ix++) {
            double al2 = x[ix+1] - x[ix];
            double gx = x[ix];
            es.push_back(Elem(depth, dipdeg, al1, al2, aw1, aw2, gx, gy));
          }
        }
      }
    }

    void DuToS(const LameParms& lp, const double* d, double* s)
    {
      double
        theta = d[0] + d[4] + d[8],
        mu = lp.mu(), lambda = lp.lambda();
      s[0] = lambda*theta + 2.0*mu*d[0];
      s[1] = mu*(d[1] + d[3]);
      s[2] = mu*(d[2] + d[6]);
      s[3] = lambda*theta + 2.0*mu*d[4];
      s[4] = mu*(d[5] + d[7]);
      s[5] = lambda*theta + 2.0*mu*d[8];      
    }

    void ProjectStress
    (const double* s, const double* normal,
     const double* along1, const double* along2,
     double* ss1, double* ss2, double* ns)
    {
      /* Equivalent to the Matlab code
           sigma = [Sxx(i) Sxy(i) Sxz(i)
                    Sxy(i) Syy(i) Syz(i)
                    Sxz(i) Syz(i) Szz(i)];
           ss(i) = along' *sigma*normal;
           ns(i) = normal'*sigma*normal; */

      // s = [sxx sxy sxz syy syz szz]. Compute Sigma*normal:
      double sn[3];
      sn[0] = s[0]*normal[0] + s[1]*normal[1] + s[2]*normal[2];
      sn[1] = s[1]*normal[0] + s[3]*normal[1] + s[4]*normal[2];
      sn[2] = s[2]*normal[0] + s[4]*normal[1] + s[5]*normal[2];
      if (along1) *ss1 = along1[0]*sn[0] + along1[1]*sn[1] + along1[2]*sn[2];
      if (along2) *ss2 = along2[0]*sn[0] + along2[1]*sn[1] + along2[2]*sn[2];
      if (ns)     *ns  = normal[0]*sn[0] + normal[1]*sn[1] + normal[2]*sn[2];
    }

  }
}
