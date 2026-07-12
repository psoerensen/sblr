#ifndef ST_LD_OPERATOR_H
#define ST_LD_OPERATOR_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include <RcppArmadillo.h>

#include "st_csr_common.h"

struct CsrOperator {
  STLDCSR ld;
  arma::rowvec xx;

  CsrOperator(const STLDCSR& ld_, const arma::rowvec& xx_) :
    ld(ld_), xx(xx_) {}

  inline const arma::rowvec& diag() const { return xx; }

  inline void apply_offdiag(int i, double diff, arma::rowvec& r) const {
    const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
    const uint64_t end = ld.ptr[static_cast<std::size_t>(i + 1)];

    for (uint64_t p = start; p < end; ++p) {
      const int j = ld.idx[static_cast<std::size_t>(p)];
      r(static_cast<arma::uword>(j)) -=
        static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
    }
  }

  inline void rebuild(
      const arma::rowvec& wy,
      const arma::rowvec& b,
      arma::rowvec& r) const {
    const int m = static_cast<int>(b.n_elem);
    r = wy;

    for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      const double bi = b(iu);
      if (bi == 0.0) continue;

      r(iu) -= xx(iu) * bi;
      apply_offdiag(i, bi, r);
    }
  }
};

struct EigenBlock {
  int start = 0;
  int size = 0;
  std::vector<float> tri;

  inline std::size_t pidx(int i, int j) const {
    const std::size_t ii = static_cast<std::size_t>(i);
    const std::size_t jj = static_cast<std::size_t>(j);
    const std::size_t ss = static_cast<std::size_t>(size);
    return ii * (2 * ss - ii + 1) / 2 + (jj - ii);
  }

  inline double sym_at(int a, int b) const {
    const int i = (a <= b) ? a : b;
    const int j = (a <= b) ? b : a;
    return static_cast<double>(tri[pidx(i, j)]);
  }
};

struct BlockEigenOperator {
  std::vector<EigenBlock> blocks;
  std::vector<int> block_of;
  std::vector<int> local_of;
  arma::rowvec diag_;

  inline const arma::rowvec& diag() const { return diag_; }

  inline void apply_offdiag(int i, double diff, arma::rowvec& r) const {
    const int g = block_of[static_cast<std::size_t>(i)];
    const int li = local_of[static_cast<std::size_t>(i)];
    const EigenBlock& B = blocks[static_cast<std::size_t>(g)];

    for (int lj = 0; lj < B.size; ++lj) {
      if (lj == li) continue;

      const double aij = B.sym_at(li, lj);
      if (aij != 0.0) {
        r(static_cast<arma::uword>(B.start + lj)) -= aij * diff;
      }
    }
  }

  inline void rebuild(
      const arma::rowvec& wy,
      const arma::rowvec& b,
      arma::rowvec& r) const {
    r = wy;

    for (const EigenBlock& B : blocks) {
      for (int li = 0; li < B.size; ++li) {
        const double bi = b(static_cast<arma::uword>(B.start + li));
        if (bi == 0.0) continue;

        for (int lj = 0; lj < B.size; ++lj) {
          r(static_cast<arma::uword>(B.start + lj)) -= B.sym_at(li, lj) * bi;
        }
      }
    }
  }
};

#endif
