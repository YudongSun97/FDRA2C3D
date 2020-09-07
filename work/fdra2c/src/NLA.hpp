#ifndef INCLUDE_FDRA_NLA
#define INCLUDE_FDRA_NLA

namespace fdra {

  // LAPACK/BLAS interface
  typedef long long blas_int;

  extern "C" void sgemm_
  (char*, char*, blas_int*, blas_int*, blas_int*, float*, const float*,
   blas_int*, const float*, blas_int*, float*, float*, blas_int*);
  extern "C" void sgemv_
  (char*, blas_int*, blas_int*, float*, const float*, blas_int*, const float*,
   blas_int*, float*, float*, blas_int*);
  extern "C" float sdot_
  (blas_int*, const float*, blas_int*, const float*, blas_int*);
  extern "C" float saxpy_
  (blas_int*, float*, const float*, blas_int*, float*, blas_int*);
  extern "C" void dgemm_
  (char*, char*, blas_int*, blas_int*, blas_int*, double*, const double*,
   blas_int*, const double*, blas_int*, double*, double*, blas_int*);
  extern "C" void dgemv_
  (char*, blas_int*, blas_int*, double*, const double*, blas_int*, const double*,
   blas_int*, double*, double*, blas_int*);
  extern "C" double ddot_
  (blas_int*, const double*, blas_int*, const double*, blas_int*);
  extern "C" double daxpy_
  (blas_int*, double*, const double*, blas_int*, double*, blas_int*);

  template<typename T> void gemm
  (char transa, char transb, blas_int m, blas_int n, blas_int k,
   T alpha, const T* A, blas_int lda, const T* B, blas_int ldb, T beta,
   T* C, blas_int ldc);
  template<typename T> void gemv
  (char trans, blas_int m, blas_int n, T alpha, const T* a, blas_int lda,
   const T* x, blas_int incx, T beta, T* y, blas_int incy);
  template<typename T> T dot
  (blas_int n, const T* x, blas_int incx, const T* y, blas_int incy);
  template<typename T> T axpy
  (blas_int n, T a, const T* x, blas_int incx, T* y, blas_int incy);

  template<> inline void gemm<double>
  (char transa, char transb, blas_int m, blas_int n, blas_int k,
   double alpha, const double* A, blas_int lda,
   const double* B, blas_int ldb, double beta,
   double* C, blas_int ldc)
  {
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
           C, &ldc);
  }

  template<> inline void gemm<float>
  (char transa, char transb, blas_int m, blas_int n, blas_int k,
   float alpha, const float* A, blas_int lda,
   const float* B, blas_int ldb, float beta,
   float* C, blas_int ldc)
  {
    sgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
           C, &ldc);
  }

  template<> inline void gemv<double>
  (char trans, blas_int m, blas_int n, double alpha,
   const double* a, blas_int lda, const double* x, blas_int incx, double beta,
   double* y, blas_int incy)
  {
    dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
  }

  template<> inline void gemv<float>
  (char trans, blas_int m, blas_int n, float alpha,
   const float* a, blas_int lda, const float* x, blas_int incx, float beta,
   float* y, blas_int incy)
  {
    sgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
  }

  template<> inline double dot<double>
  (blas_int n, const double* x, blas_int incx, const double* y, blas_int incy)
  {
    return ddot_(&n, x, &incx, y, &incy);
  }

  template<> inline float dot<float>
  (blas_int n, const float* x, blas_int incx, const float* y, blas_int incy)
  {
    return sdot_(&n, x, &incx, y, &incy);
  }

  template<> inline double axpy<double>
  (blas_int n, double a, const double* x, blas_int incx, double* y, blas_int incy)
  {
    return daxpy_(&n, &a, x, &incx, y, &incy);
  }

  template<> inline float axpy<float>
  (blas_int n, float a, const float* x, blas_int incx, float* y, blas_int incy)
  {
    return saxpy_(&n, &a, x, &incx, y, &incy);
  }

};

#endif
