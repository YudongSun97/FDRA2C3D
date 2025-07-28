#ifndef INCLUDE_FDRA_HMATRIXSTRESSFN
#define INCLUDE_FDRA_HMATRIXSTRESSFN

#include <vector>
#include "MyMpiHmat.hpp"
#include "Fdra.hpp"

namespace fdra {
  using namespace util;

  typedef hmmvp::MpiHmat<tau_real> MpiHmat;

  class HmatrixStressFn : public StressFn {
  public:
    HmatrixStressFn(const mpi::ArraySegmenter* as, const vector<MpiHmat*>& hm,
                    int ncomp);
    virtual ~HmatrixStressFn();
    //Yudong May 6 2025
    //virtual bool IncludeNormalComponent() { return false; }
    virtual bool IncludeNormalComponent () { return _IncludeNormal; }

    virtual void Call(int deriv, double t, const real* x,
                      tau_real* tau, tau_real* taun) = 0;

  protected:
    virtual int NeedNrwrk() = 0;
    //Yudong May 6 2025
    bool _IncludeNormal;
  };

  StressFn* NewHmatrixStressFn(const Model* m, KeyValueFile* kvf)
    throw (FileException, Exception);

}

#endif
