#include <stdio.h>
#include <string.h>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include "util/include/Mpi.hpp"
#include "NLA.hpp"
#include "Ode.hpp"
#include "DDouble.hpp"
using namespace std;
using namespace util;

namespace fdra {

  static const double __meps__ = numeric_limits<double>::epsilon();
  // I'm using eps^2 because t is DDouble.
  const double OdeMethod::_eps = __meps__ * __meps__;
  const double Ode23::_pow = 1.0/3.0;
  const double Ode23::_A[3] = { 0.5, 0.75, 1.0 };
  const double Ode23::_B[12] = { 0.5, 0.0, 0.0, 0.0,
                                 0.0, 0.75, 0.0, 0.0,
                                 2.0/9.0, 1.0/3.0, 4.0/9.0, 0.0 };
  const double Ode23::_E[4] = { -5.0/72.0, 1.0/12.0, 1.0/9.0, -1.0/8.0 };

  OdeMethod::OdeMethod ()
    : _abs_tol(1.0e-6), _rel_tol(1.0e-3),
      _initial_step(1.0e-6), _ti(0), _tf(0), _have_ts(false)
  {}

  bool OdeMethod::SetAbsTol (double abs_tol) {
    if (abs_tol < 0.0) return false;
    _abs_tol = abs_tol;
    return true;
  }

  bool OdeMethod::SetRelTol (double rel_tol) {
    if (rel_tol < _eps) return false;
    _rel_tol = rel_tol;
    return true;
  }

  bool OdeMethod::SetInitialStep (double initial_step) {
    if (initial_step <= 0.0) return false;
    _initial_step = initial_step;
    return true;
  }

  bool OdeMethod::SetTspan (double ti, double tf) {
    _ti = ti;
    _tf = tf;
    return true;
  }

  struct Ode23Data {
    real *f0, *f1, *f2, *f3, *y, *yt, *w;

    Ode23Data (int n) throw (OdeException, OutOfMemoryException) {
      if (n <= 0) throw OdeException("n <= 0.");
      _data.resize(7*n);
      real* d = &_data[0];
      f0 = d;
      f1 = d + n;
      f2 = f1 + n;
      f3 = f2 + n;
      y = f3 + n;
      yt = y + n;
      w = yt + n;
    }

  private:
    vector<real> _data;
  };

  /* I'm representing t internally by a double-double number. This allows dt to
     be accurate, which I want for use in diffusion calculations. double(t) is
     passed to the user, however. */
  void Ode23::Run (OdeFunction& ofn) throw (OdeException) {
    const double threshold = _abs_tol / _rel_tol;
    const double tdir = (_tf > _ti) ? 1.0 : -1.0;

    if (mpi::AmRoot()) printf("rel_tol: %1.2e\n", _rel_tol);

    const int n = ofn.GetN();
    Ode23Data d(n);
    ofn.GetIc(d.y);

    DDouble t = _ti;
    double absh = _initial_step;
    bool is_error = false;
    ofn.Call(t, 0, d.y, d.f0, is_error);
    if (is_error) throw OdeException("is_error set at t = t0.");
    if (!ofn.ViewOutput(t, 0, d.y, OdeFunction::m_init)) return;

    while (t < _tf) {
      // 16 is a magic number in matlab's ode23. I have no reason to change it.
      const double hmin = 16.0*_eps*t;
      bool nofailed = true;
      double tnew, err, h;
      for (;;) { // Loop for one step
        h = tdir*absh;
        if (t + h > _tf) {
          h = _tf - t;
          absh = abs(h);
        }
        tnew = t + h;
        h = tnew - t;

        err = numeric_limits<double>::infinity();
        while (true) { // Not really a loop; just allow break-out.
          // yt = y + f h B(:,1)
          memcpy(d.yt, d.y, n*sizeof(real));
          axpy<real>(n, h*_B[0], d.f0, 1, d.yt, 1);
          double dt = h*_A[0];
          ofn.Call(t + dt, dt, d.yt, d.f1, is_error);
          if (is_error) break;

          // yt = y + f h B(:,2)
          memcpy(d.yt, d.y, n*sizeof(real));
          axpy<real>(n, h*_B[5], d.f1, 1, d.yt, 1);
          dt = h*_A[1];
          ofn.Call(t + dt, dt, d.yt, d.f2, is_error);
          if (is_error) break;

          // yt = y + f h B(:,3)
          memcpy(d.yt, d.y, n*sizeof(real));
          axpy<real>(n, h*_B[8 ], d.f0, 1, d.yt, 1);
          axpy<real>(n, h*_B[9 ], d.f1, 1, d.yt, 1);
          axpy<real>(n, h*_B[10], d.f2, 1, d.yt, 1);
          ofn.Call(tnew, h, d.yt, d.f3, is_error);
          if (is_error) break;

          // err = absh*norm((f*E)./max(max(abs(y),abs(yt)),threshold),inf)
          
          memset(d.w, 0, n*sizeof(real));
          axpy<real>(n, _E[0], d.f0, 1, d.w, 1);
          axpy<real>(n, _E[1], d.f1, 1, d.w, 1);
          axpy<real>(n, _E[2], d.f2, 1, d.w, 1);
          axpy<real>(n, _E[3], d.f3, 1, d.w, 1);
          double my_err = 0.0;
          for (int i = 0; i < n; i++)
            my_err = max
              (my_err,
               abs(d.w[i]) / max((double) max(abs(d.y[i]), abs(d.yt[i])),
                                 threshold));
          my_err *= absh;
          mpi::Allreduce(&my_err, &err, 1, MPI_MAX);
          break;
        }
        if (is_error) nofailed = false;

        if (err > _rel_tol) {
          if (absh <= hmin) {
            ofn.ViewOutput(t, h, d.y, OdeFunction::m_tolfail);
            return;
          }
          if (nofailed) {
            nofailed = false;
            absh = max(hmin, absh*max(0.5, 0.8*pow(_rel_tol / err, _pow)));
          } else {
            absh = max(hmin, 0.5*absh);
          }
        } else {
          break;
        }
      } // Loop for one step

      if (nofailed)
        absh *= min(5.0, 0.8*pow(_rel_tol / err, _pow));

      t = tnew;
      swap(d.y, d.yt);
      swap(d.f0, d.f3);
      if (!ofn.ViewOutput(t, h, d.y, OdeFunction::m_none)) return;
    } // t < _tf
  }

};
