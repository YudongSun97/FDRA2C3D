#ifndef INCLUDE_FDRA_DDOUBLE
#define INCLUDE_FDRA_DDOUBLE

namespace fdra {

  class DDouble {
  public:
    DDouble(double x_) : x(x_), y(0.0) {}
    DDouble(double x_, double y_) : x(x_), y(y_) {}
    
    DDouble& operator+=(const DDouble& b)
    {
      DDouble& a = *this;
      // Adapted from ddadd in David H. Bailey's DDFUN. This implementation is
      // based on older work. In particular:
      double t1 = a.x + b.x;
#if 1
      //   Two-Sum(a.x, b.x) in Shewchuk's "Robust Geometric Predicates" paper,
      // attributed to one of D.E. Knuth's books. Two-Sum is a branch-less
      // alternative to a branching call to Fast-Two-Sum(a.x, a.b) based on
      // |a.x| >= |b.x|.
      double e = t1 - a.x;
      double t2 = ((b.x - e) + (a.x - (t1 - e)));
#else
      //   Branch-based Two-Sum.
      double t2;
      // |a.x| > |b.x|; an idea in Shewchuck's paper
      if ((a.x > b.x) == (a.x > -b.x)) {
        // Fast-Two-Sum(a.x , b.x)
        t2 = b.x - (t1 - a.x);
      } else {
        // Fast-Two-Sum(b.x , a.x)
        t2 = a.x - (t1 - b.x);
      }
#endif
      // Now accumulate the low-order parts.
      t2 += a.y + b.y;
      //   Fast-Two-Sum(t1, t2) in Shewchuk's paper, attributed to
      //     T.J. Dekker 1971.
      // Applicable because |t1| >= |t2|.
      a.x = t1 + t2;
      a.y = t2 - (a.x - t1);
      return a;
    }

    DDouble& operator+=(const double bx)
    {
      DDouble& a = *this;
      // This is the previous method specialized to the case b.y = 0.
      double t1 = a.x + bx;
      double e = t1 - a.x;
      double t2 = (((bx - e) + (a.x - (t1 - e)))) + a.y;
      a.x = t1 + t2;
      a.y = t2 - (a.x - t1);
      return a;
    }

    DDouble& operator*=(const double s)
    {
      x *= s;
      y *= s;
      return *this;
    }

    DDouble& operator-=(const DDouble& b_)
    {
      DDouble b = b_;
      b *= -1.0;
      return operator+=(b);
    }

    double ToDouble() const
    {
      return x;
    }

    operator double() const
    {
      return x;
    }

  private:
    double x, y;
  };

  DDouble operator+(const DDouble& a, const DDouble& b)
  {
    DDouble c = a;
    c += b;
    return c;
  }

  DDouble operator-(const DDouble& a, const DDouble& b)
  {
    DDouble c = a;
    c -= b;
    return c;
  }

  DDouble operator+(const DDouble& a, const double b)
  {
    DDouble c = a;
    c += b;
    return c;
  }

  DDouble operator-(const double a, const DDouble& b)
  {
    DDouble c = b;
    c *= -1.0;
    c += a;
    return c;
  }

}

#endif
