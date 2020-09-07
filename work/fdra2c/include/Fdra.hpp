#ifndef INCLUDE_FDRA_FDRA
#define INCLUDE_FDRA_FDRA

#include "util/include/Defs.hpp"
#include "util/include/Mpi.hpp"
#include "util/include/KeyValueFile.hpp"
#include "StreamDataFile.hpp"
#include "Defs.hpp"

namespace fdra {
  using namespace util;

  class StressFn {
  public:
    virtual ~StressFn() {}

    // Does not always have to return the same value.
    virtual bool IncludeNormalComponent() = 0;

    // tau is NULL if only the normal component is wanted.
    virtual void Call(int deriv, double t, const real* x,
                      tau_real* tau, tau_real* taun) = 0;
  };

  class StateVector {
  public:
    int Ncomp() const;
    const real* slip() const;
    const real* v() const;
    const real* theta() const;
    const real* dlte() const;
    const real* p() const;
    const real* T() const;
  };

  class OutputListener {
  public:
    OutputListener(const mpi::ArraySegmenter* as) : _as(as) {}

    virtual ~OutputListener() {}

    virtual bool Call(double t, double dt, const StateVector* sv) = 0;

  protected:
    const mpi::ArraySegmenter* _as;
  };

  class Model
  {
  public:
    int GetN() const;
    int GetNcomp() const;
    int GetNelem() const;

    // Report whether the current state of the model is consistent and
    // sufficient to use in the simulation. It is not safe to call any other
    // method until IsOk returns true.
    bool IsOk(bool disp_messages = true) const;
    // Initialize the listeners only after IsOk() returns true.
    bool InitListeners(bool disp_messages = true);
    // Put contents in a KeyValueFile. Only the root process needs a non-NULL
    // kvf.
    void ToKvf(KeyValueFile* kvf) const;

    const mpi::ArraySegmenter* GetArraySegmenter() const;
    void SetStressFn(StressFn* sf);
    StressFn* GetStressFn() const;

  private:
    explicit Model();
    Model(const Model&);
    Model& operator=(const Model&);
  };

  // Only the root process needs a non-NULL kvf.
  Model* BuildModelFromKeyValueFile(const KeyValueFile* kvf,
                                    bool disp_messages = true);
  void DeleteModel(Model* m);

  void Go(const Model* m);

}

#endif
