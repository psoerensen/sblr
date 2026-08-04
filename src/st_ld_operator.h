#ifndef ST_LD_OPERATOR_H
#define ST_LD_OPERATOR_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include <RcppArmadillo.h>

#include "st_csr_common.h"
#include "blr_block_eigen.h"
#include "blr_block_low_rank.h"

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

  inline arma::uword residual_size() const { return xx.n_elem; }

  inline bool uses_retained_low_rank() const { return false; }

  inline const char* diagnostic_name() const { return "csr"; }

  inline double rebuild_and_measure_drift(int, const arma::rowvec&,
                                          const arma::rowvec&,
                                          arma::rowvec&) const {
    return 0.0;
  }

  inline double corrected_rhs(int i, double beta_old, const arma::rowvec& r) const {
    const arma::uword iu = static_cast<arma::uword>(i);
    return r(iu) + xx(iu) * beta_old;
  }

  inline void apply_difference(int i, double diff, arma::rowvec& r) const {
    const arma::uword iu = static_cast<arma::uword>(i);
    r(iu) -= xx(iu) * diff;
    apply_offdiag(i, diff, r);
  }

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

  inline void rebuild(int, const arma::rowvec& wy, const arma::rowvec& b,
                      arma::rowvec& r) const { rebuild(wy, b, r); }

  inline double projected_score_dot(int, const arma::rowvec& b,
                                    const arma::rowvec& wy) const {
    double bwy = 0.0;
    for (arma::uword i = 0; i < b.n_elem; ++i) bwy += b(i) * wy(i);
    return bwy;
  }

  inline double quadratic_form(const arma::rowvec& b) const {
    // Validation/reference operation. Ordinary operator-aware MCMC iterations
    // use fitted_quadratic() from the maintained residual state.
    arma::rowvec zero(xx.n_elem, arma::fill::zeros);
    arma::rowvec rb;
    rebuild(zero, b, rb);
    return -arma::dot(b, rb);
  }

  inline double fitted_quadratic(int, const arma::rowvec& b,
                                 const arma::rowvec& wy,
                                 const arma::rowvec& r) const {
    double result = 0.0;
    const double* b_ptr = b.memptr();
    const double* wy_ptr = wy.memptr();
    const double* r_ptr = r.memptr();
    for (arma::uword i = 0; i < b.n_elem; ++i)
      result += b_ptr[i] * (wy_ptr[i] - r_ptr[i]);
    return result;
  }

  inline double residual_sse(int, double yy, const arma::rowvec& b,
                             const arma::rowvec& wy, const arma::rowvec& r) const {
    double b_dot_r_plus_wy = 0.0;
    for (arma::uword i = 0; i < b.n_elem; ++i)
      b_dot_r_plus_wy += b(i) * (r(i) + wy(i));
    return yy - b_dot_r_plus_wy;
  }

  inline double genetic_variance(int trait, const arma::rowvec& b,
                                 const arma::rowvec& wy, const arma::rowvec& r,
                                 double n) const {
    return fitted_quadratic(trait, b, wy, r) / n;
  }

  inline void materialize_residual(int, const arma::rowvec& r,
                                   arma::rowvec& marker_residual) const {
    marker_residual = r;
  }
};

using EigenBlock = sblr::core::BlockEigenBlockStorage;
using BlockEigenOperator = sblr::core::BlockEigenStorage;
using BlockLowRankOperator = sblr::core::BlockLowRankStorage;

struct BlockEigenDispatchOperator {
  bool low_rank = false;
  BlockEigenOperator dense;
  BlockLowRankOperator retained;

  inline const arma::rowvec& diag() const {
    return low_rank ? retained.diag() : dense.diag();
  }
  inline arma::uword residual_size() const {
    return low_rank ? retained.residual_size() : dense.residual_size();
  }
  inline bool uses_retained_low_rank() const { return low_rank; }
  inline const char* diagnostic_name() const {
    return low_rank ? "retained_block_eigen" : "dense_block_eigen";
  }
  inline double corrected_rhs(int i, double beta, const arma::rowvec& r) const {
    return low_rank ? retained.corrected_rhs(i, beta, r) : dense.corrected_rhs(i, beta, r);
  }
  inline void apply_difference(int i, double difference, arma::rowvec& r) const {
    if (low_rank) retained.apply_difference(i, difference, r);
    else dense.apply_difference(i, difference, r);
  }
  inline void rebuild(int trait, const arma::rowvec& wy, const arma::rowvec& b,
                      arma::rowvec& r) const {
    if (low_rank) retained.rebuild(trait, wy, b, r);
    else dense.rebuild(trait, wy, b, r);
  }
  inline double rebuild_and_measure_drift(int trait, const arma::rowvec& wy,
                                          const arma::rowvec& b,
                                          arma::rowvec& r) const {
    if (!low_rank) return 0.0;
    return retained.rebuild_and_measure_drift(trait, wy, b, r);
  }
  inline double projected_score_dot(int trait, const arma::rowvec& b,
                                    const arma::rowvec& wy) const {
    return low_rank ? retained.projected_score_dot(trait, b, wy) :
      dense.projected_score_dot(trait, b, wy);
  }
  inline double quadratic_form(const arma::rowvec& b) const {
    return low_rank ? retained.quadratic_form(b) : dense.quadratic_form(b);
  }
  inline double fitted_quadratic(int trait, const arma::rowvec& b,
                                 const arma::rowvec& wy,
                                 const arma::rowvec& r) const {
    return low_rank ? retained.fitted_quadratic(trait, b, wy, r) :
      dense.fitted_quadratic(trait, b, wy, r);
  }
  inline double residual_sse(int trait, double yy, const arma::rowvec& b,
                             const arma::rowvec& wy, const arma::rowvec& r) const {
    return low_rank ? retained.residual_sse(trait, yy, b, wy, r) :
      dense.residual_sse(trait, yy, b, wy, r);
  }
  inline double genetic_variance(int trait, const arma::rowvec& b,
                                 const arma::rowvec& wy, const arma::rowvec& r,
                                 double n) const {
    return low_rank ? retained.genetic_variance(trait, b, wy, r, n) :
      dense.genetic_variance(trait, b, wy, r, n);
  }
  inline void materialize_residual(int trait, const arma::rowvec& r,
                                   arma::rowvec& marker_residual) const {
    if (low_rank) retained.materialize_residual(trait, r, marker_residual);
    else dense.materialize_residual(trait, r, marker_residual);
  }
};

inline void sampleE_ST_operator_from_scale(
    double nue, double& ve, double residual_scale, int n, std::mt19937& gen) {
  if (!std::isfinite(residual_scale) || residual_scale <= 0.0)
    throw std::runtime_error("sampleE_ST_operator: invalid projected residual scale.");
  std::chi_squared_distribution<double> rchisq(n + nue);
  const double chi2 = std::max(rchisq(gen), 1e-300);
  const double ve_new = residual_scale / chi2;
  if (!std::isfinite(ve_new) || ve_new <= 0.0)
    throw std::runtime_error("sampleE_ST_operator: sampled ve is invalid.");
  ve = std::max(ve_new, 1e-12);
}

template <typename OpT>
inline void sampleE_ST_operator(
    const OpT& op, int trait, double nue, double& ve,
    const arma::rowvec& b, const arma::rowvec& wy, const arma::rowvec& r,
    double sse_prior, double yy, int n, std::mt19937& gen) {
  const double sse = op.residual_sse(trait, yy, b, wy, r);
  const double scale = sse + nue * sse_prior;
  sampleE_ST_operator_from_scale(nue, ve, scale, n, gen);
}

template <typename OpT>
inline double computeG_ST_operator(
    const OpT& op, int trait, const arma::rowvec& b,
    const arma::rowvec& wy, const arma::rowvec& r, int n) {
  return op.genetic_variance(trait, b, wy, r, static_cast<double>(n));
}

#endif
