#ifndef INCLUDE_FDRA_ODE
#define INCLUDE_FDRA_ODE

// Clone of Matlab's ode23 with a few different options. Cite
//   P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
//   Appl. Math Lett. 2(4), 321-325, 1989.

#include <vector>
#include "util/include/Exception.hpp"
#include "util/include/Mpi.hpp"
#include "Defs.hpp"

namespace fdra {
  using namespace std;
  using namespace util;

  class OdeException : public Exception {
  public:
    OdeException(const std::string& msg = "odee") : Exception(msg) {}
  };

  // All data are for this MPI process.
  class OdeFunction {
  public:
    enum Message { m_none, m_init, m_final, m_tolfail };

    virtual ~OdeFunction() {}

    // Problem size.
    virtual int GetN() = 0;
    virtual void GetIc(real* y0) = 0;

    // The ODE function. Set yd to the time derivative of y at time t. Set
    // is_error = true if there is a problem. If MPI is being used, every
    // process must set is_error to the same value.
    virtual void Call
      (double t, double dt, const real* y, real* yd, bool& is_error) = 0;

    // Called at the end of each successful step. Use this method to save
    // data. Return true to continue, false to stop.
    virtual bool ViewOutput(double t, double dt, const real* y, Message m)
    { return true; }

    // Called at the end of a step. Optionally set certain elements to values
    // you want. This is useful when solving an index-1 DAE. You can include
    // equations the REC monitors, but then set the values at the end of each
    // step (equivalently, beginning of the next step) to the values you obtain
    // by solving algebraic equations.
    virtual void ResetY(real* y) {};
  };

  class OdeMethod {
  public:
    OdeMethod();
    virtual ~OdeMethod() {}

    bool SetAbsTol(double abs_tol);
    bool SetRelTol(double rel_tol);
    bool SetInitialStep(double initial_step);

    double GetAbsTol() { return _abs_tol; }
    double GetRelTol() { return _rel_tol; }
    double GetInitialStep() { return _initial_step; }

    bool SetTspan(double ti, double tf);
    bool SetTspan(const vector<double>& ts) { return false; }

    // y0 is needed only until before the first call to an ofn method, so you
    // can use temporary memory for it.
    virtual void Run(OdeFunction& ofn) throw (OdeException) = 0;

  protected:
    double _abs_tol;
    double _rel_tol;
    double _initial_step;

    double _ti, _tf;
    bool _have_ts;
    vector<double> _ts;

    static const double _eps;
  };

  class Ode23 : public OdeMethod {
  public:

    virtual void Run(OdeFunction& ofn) throw (OdeException);

  private:
    
    static const double _pow, _A[3], _B[12], _E[4];
  };

};

#endif
