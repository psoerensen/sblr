#ifndef ST_LD_OPERATOR_H
#define ST_LD_OPERATOR_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include <RcppArmadillo.h>

#include "st_csr_common.h"
#include "blr_block_eigen.h"

struct CsrOperator {
  STLDCSR ld;
  arma::rowvec xx;

  CsrOperator(const STLDCSR& ld_, const arma::rowvec& xx_) :
    ld(ld_), xx(xx_) {}

  inline sblr::core::SparseLdCsrView view() const {
    sblr::core::SparseLdCsrView result;
    result.marker_count = static_cast<std::size_t>(xx.n_elem);
    result.row_ptr = ld.ptr.data();
    result.row_ptr_size = ld.ptr.size();
    result.column_index = ld.idx.empty() ? nullptr : ld.idx.data();
    result.offdiag_xij = ld.xij.empty() ? nullptr : ld.xij.data();
    result.nonzero_count = ld.idx.size();
    result.diagonal = &xx;
    return result;
  }

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

using EigenBlock = sblr::core::BlockEigenBlockStorage;
using BlockEigenOperator = sblr::core::BlockEigenStorage;

#endif
