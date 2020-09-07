#include <sys/time.h>

typedef double real;

static inline double difftime(struct timeval& t1, struct timeval& t2)
{
  static const double us = 1.0e6;
  return (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
}

static inline real DoSomeSillyStuff(real v, int i)
{
  // I'm not even joking.
  v = (v + 0.1);
  v *= v;
  v = pow(v, 1.5);
  v = pow(v, 0.68);
  v *= i / (i + 0.1);
  return v;
}
