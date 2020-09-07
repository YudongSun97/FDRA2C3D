#ifndef INCLUDE_HMMVP_HD
#define INCLUDE_HMMVP_HD

#include <vector>
#include "util/include/Defs.hpp"
#include "util/include/Exception.hpp"
#include "util/include/Matrix.hpp"

namespace hmmvp {
  using namespace util;

  // Row and column global indices for a matrix block. Indexing is 1-based.
  struct Block {
    uint r0, m, c0, n;
  };

  class Hd {
  private:
    virtual ~Hd() {}

  public:
    // Permutations use base-1 indexing.
    virtual void Permutations(std::vector<uint>& p, std::vector<uint>& q)
      const = 0;

    // Methods to access matrix blocks
    typedef std::vector<Block>::const_iterator iterator;
    iterator Begin() const;
    iterator End() const;
    std::vector<Block>::size_type NbrBlocks() const;

  private:
    Hd();
    Hd(const Hd&);
    Hd& operator=(const Hd&);
  };

  // D and D are 3xN matrices. R can be empty or not provided if R = D.
  Hd* NewHd(const Matrix<double>& D);
  Hd* NewHd(const Matrix<double>& D, const Matrix<double>& R);

  void WriteHd(const Hd* hd, const std::string& hd_filename)
    throw (FileException);

  Hd* NewHd(const std::string& hd_filename)
    throw (Exception, FileException);

  void DeleteHd(Hd* hd);

}

#endif
