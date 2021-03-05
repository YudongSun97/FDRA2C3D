#ifndef INCLUDE_UTIL_WORKARRAY
#define INCLUDE_UTIL_WORKARRAY

#include <vector>

namespace util {
  using namespace std;

  template<typename T> class WorkArray {
  public:
    WorkArray() : _ptr(0), _cap(0) {}
    WorkArray(size_t n) : _w(n), _ptr(0), _cap(n) {}
    WorkArray(size_t n, const T& t) : _w(n, t), _ptr(0), _cap(n) {}

    void Reset(size_t n)
    {
      if (n > _cap) {
        _w.resize(n);
        _cap = n;
      }
      _ptr = 0;
    }

    void Reset() { _ptr = 0; }

    T* GetPtr() { return &_w[0]; }

    T* AllocWork(size_t n)
    {
      size_t sz = _ptr + n;
      if (sz > _cap) return NULL;
      T* w = &_w[0] + _ptr;
      _ptr = sz;
      return w;
    }

  private:
    vector<T> _w;
    size_t _ptr;
    unsigned int _cap;
  };

}

#endif
