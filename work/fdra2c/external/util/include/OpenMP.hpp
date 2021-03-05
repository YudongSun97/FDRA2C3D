#ifndef INCLUDE_UTIL_OPENMP
#define INCLUDE_UTIL_OPENMP

#ifdef UTIL_OMP
#  include <omp.h>
#else
#  define omp_set_num_threads(nthreads) do {} while (0)
#  define omp_get_max_threads() 1
#  define omp_get_thread_num() 0
#  define omp_init_lock(a) do {} while (0)
#  define omp_destroy_lock(a) do {} while (0)
#  define omp_set_lock(a) do {} while (0)
#  define omp_unset_lock(a) do {} while (0)
typedef char omp_lock_t;
#endif

#endif
