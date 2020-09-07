#ifndef INCLUDE_FDRA_ODEFN
#define INCLUDE_FDRA_ODEFN

#include <vector>
#include <string>
#include "util/include/Defs.hpp"
#include "util/include/Mpi.hpp"
#include "util/include/WorkArray.hpp"
#include "Ode.hpp"

namespace fdra {
  using namespace std;
  using namespace util;

  class Model;
  class StressFn;

  class OdeStateVector {
  public:
    const real *slip, *v_o_v0, *theta_v0_o_dc, *dlte_v0_o_dc;

    OdeStateVector(const Model* m);

    bool Decode(const real* y);

  private:
    int _nelem, _n;
    bool _use_only_theta, _use_dilatancy;
  };

  class MembraneDiffusionIntegrator;

  class OdeFn : public OdeFunction {
  public:
    OdeFn(const Model* m);
    ~OdeFn();

    virtual int GetN();
    virtual void GetIc(real* y0);
    virtual void Call(double t, double dt, const real* y, real* yd,
                      bool& is_error);
    virtual bool ViewOutput(double t, double dt, const real* y, Message m);
    virtual void ResetY(real* y);

  private:
    const Model* _m;
    OdeStateVector _osv;
    MembraneDiffusionIntegrator* _mdi;
    WorkArray<real> _rwrk;
    WorkArray<tau_real> _trwrk;
    int _ctr, _disp_every;
    struct StopIndicator {
      bool check;
      int frequency;
      string indicator;
    } _si;
  };

  void RunOde(const Model* m);
}

#endif
