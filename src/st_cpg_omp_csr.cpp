// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "st_chain_utils.h"
#include "st_csr_common.h"
#include "st_block_eigen.h"
#include "st_block_eigen_execution.h"
#include "st_block_eigen_rcpp.h"
#include "st_ld_operator.h"
#include "blr_csr_bayesc_operator_adapter.h"
#include "blr_csr_bayesc_policy.h"
#include "blr_csr_bayesc_rcpp_adapter.h"
#include "blr_phase3_execution.h"
#define SBLR_CSR_BAYESC_CORE_IMPL_TRANSLATION_UNIT 1
// Emit the ordinary no-op-policy entry point from the reusable engine.
#include "blr_csr_bayesc_core_impl.h"
#undef SBLR_CSR_BAYESC_CORE_IMPL_TRANSLATION_UNIT


#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

static std::vector<std::string> stblr_copy_character_vector(
  Rcpp::CharacterVector x,
  const char* label
) {
 std::vector<std::string> out(static_cast<std::size_t>(x.size()));
 for (int i = 0; i < x.size(); ++i) {
  if (x[i] == NA_STRING) {
   throw std::runtime_error(std::string(label) + " contains NA.");
  }
  out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(x[i]);
  if (out[static_cast<std::size_t>(i)].empty()) {
   throw std::runtime_error(std::string(label) + " contains an empty string.");
  }
 }
 return out;
}

static std::vector<std::vector<int>> stblr_copy_int_list(
  Rcpp::List x,
  const char* label
) {
 std::vector<std::vector<int>> out(static_cast<std::size_t>(x.size()));
 for (int k = 0; k < x.size(); ++k) {
  Rcpp::IntegerVector v = x[k];
  out[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(v.size()));
  for (int i = 0; i < v.size(); ++i) {
   if (v[i] == NA_INTEGER) {
    throw std::runtime_error(std::string(label) + " contains NA.");
   }
   if (v[i] <= 0) {
    throw std::runtime_error(std::string(label) + " must contain positive 1-based marker indices.");
   }
   out[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = v[i];
  }
 }
 return out;
}

static std::vector<int> stblr_copy_rows0_or_empty(
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  int n_bed
) {
 std::vector<int> out;
 if (rows.isNull()) return out;

 Rcpp::IntegerVector r(rows);
 out.resize(static_cast<std::size_t>(r.size()));
 bool identity_rows = (r.size() == n_bed);
 for (int i = 0; i < r.size(); ++i) {
  if (r[i] == NA_INTEGER) {
   throw std::runtime_error("rows contains NA.");
  }
  if (r[i] < 1 || r[i] > n_bed) {
   throw std::runtime_error("rows contains an index outside [1, n].");
  }
  if (r[i] != i + 1) identity_rows = false;
  out[static_cast<std::size_t>(i)] = r[i] - 1;
 }
 if (identity_rows) out.clear();
 return out;
}

static std::vector<double> stblr_copy_numeric_vector(
  Rcpp::NumericVector x,
  const char* label
) {
 std::vector<double> out(static_cast<std::size_t>(x.size()));
 for (int i = 0; i < x.size(); ++i) {
  const double val = x[i];
  if (!std::isfinite(val)) {
   throw std::runtime_error(std::string(label) + " must contain finite values.");
  }
  out[static_cast<std::size_t>(i)] = val;
 }
 return out;
}

static std::vector<int> stblr_copy_integer_vector(
  Rcpp::IntegerVector x,
  const char* label
) {
 std::vector<int> out(static_cast<std::size_t>(x.size()));
 for (int i = 0; i < x.size(); ++i) {
  if (x[i] == NA_INTEGER) {
   throw std::runtime_error(std::string(label) + " contains NA.");
  }
  out[static_cast<std::size_t>(i)] = x[i];
 }
 return out;
}


// -----------------------------------------------------------------------------
// ST-specific LD structure: flat symmetric CSR
// Stores pre-scaled X_i'X_j, not raw LD correlation.
// Disk input is expected to be upper-triangular or otherwise non-symmetric CSR
// with raw correlations r_ij. This builder symmetrizes it.
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Single-trait BayesC marker update
// -----------------------------------------------------------------------------

template <class OpT>
inline void sampleBetaC_ST_csr(
  int i,
  const std::vector<double>& pi,
  double vb,
  const arma::rowvec& prior_scale,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const OpT& op,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);

 const double wi = ww(iu);

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm01(0.0, 1.0);

 const double pi0 = std::max(pi[0], 1e-300);
 const double pi1 = std::max(pi[1], 1e-300);

 const double vei_safe = std::max(vei_i, 1e-300);
 const double vbi = std::max(vb * prior_scale(iu), 1e-300);

 // score = x_i' residual_without_i
 const double score = op.corrected_rhs(i, b(iu), r);

 // Same BayesC scalar marginal likelihood as old sbayes(),
 // but using scalar adjusted residual variance vei_i.
 const double denom = std::max(vei_safe + wi * vbi, 1e-300);

 const double logBF =
  0.5 * std::log(vei_safe / denom)
  + 0.5 * score * score * vbi / (vei_safe * denom);

 const double logp1 = std::log(pi1) + logBF;
 const double logp0 = std::log(pi0);

 const double delta_log = logp0 - logp1;

 double p1 = 0.0;

 if (delta_log > 35.0) {
  p1 = 0.0;
 } else if (delta_log < -35.0) {
  p1 = 1.0;
 } else {
  p1 = 1.0 / (1.0 + std::exp(delta_log));
 }

 const int di = (runif(gen) < p1) ? 1 : 0;

 double b_new = 0.0;

 if (di == 1) {
  const double lhs = wi + vei_safe / vbi;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);

  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  op.apply_difference(i, diff, r);
 }

 b(iu) = b_new;
 d(iu) = di;
}

template <class OpT>
inline void sampleBetaC_ST_csr_unscaled(
  int i,
  const std::vector<double>& pi,
  double vb,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const OpT& op,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);

 const double wi = ww(iu);

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm01(0.0, 1.0);

 const double pi0 = std::max(pi[0], 1e-300);
 const double pi1 = std::max(pi[1], 1e-300);

 const double vei_safe = std::max(vei_i, 1e-300);
 const double vbi = std::max(vb, 1e-300);

 // score = x_i' residual_without_i
 const double score = op.corrected_rhs(i, b(iu), r);

 // Same BayesC scalar marginal likelihood as old sbayes(),
 // but using scalar adjusted residual variance vei_i.
 const double denom = std::max(vei_safe + wi * vbi, 1e-300);

 const double logBF =
  0.5 * std::log(vei_safe / denom)
  + 0.5 * score * score * vbi / (vei_safe * denom);

 const double logp1 = std::log(pi1) + logBF;
 const double logp0 = std::log(pi0);

 const double delta_log = logp0 - logp1;

 double p1 = 0.0;

 if (delta_log > 35.0) {
  p1 = 0.0;
 } else if (delta_log < -35.0) {
  p1 = 1.0;
 } else {
  p1 = 1.0 / (1.0 + std::exp(delta_log));
 }

 const int di = (runif(gen) < p1) ? 1 : 0;

 double b_new = 0.0;

 if (di == 1) {
  const double lhs = wi + vei_safe / vbi;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);

  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  op.apply_difference(i, diff, r);
 }

 b(iu) = b_new;
 d(iu) = di;
}

inline double bayesc_maf_effect_scale(double s, double log_h) {
 return std::exp((s + 1.0) * log_h);
}

inline void fill_maf_effect_s_prior_scale(
  int m,
  double s,
  const arma::rowvec& log_h,
  arma::rowvec& prior_scale
) {
 prior_scale.set_size(static_cast<arma::uword>(m));
 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double scale_i = bayesc_maf_effect_scale(s, log_h(iu));
  if (!std::isfinite(scale_i) || scale_i <= 0.0) {
   throw std::runtime_error("dynamic maf_effect_s prior scale became invalid.");
  }
  prior_scale(iu) = scale_i;
 }
}

// -----------------------------------------------------------------------------
// Single-trait variance updates
// -----------------------------------------------------------------------------

inline void sampleB_ST_csr(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  const arma::rowvec& prior_scale,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (d(iu) > 0) {
   ssb += b(iu) * b(iu) / prior_scale(iu);
   dfb += 1.0;
  }
 }

 const double scale = ssb + nub * ssb_prior;

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);

 vb = std::max(scale / chi2, 1e-12);
}

inline void sampleB_ST_csr_unscaled(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (d(iu) > 0) {
   ssb += b(iu) * b(iu);
   dfb += 1.0;
  }
 }

 const double scale = ssb + nub * ssb_prior;

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);

 vb = std::max(scale / chi2, 1e-12);
}

inline double logpost_maf_effect_s_bayesc(
  double s,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  double vb,
  const arma::rowvec& log_h,
  double prior_lower,
  double prior_upper
) {
 if (!std::isfinite(s) || s < prior_lower || s > prior_upper) {
  return -std::numeric_limits<double>::infinity();
 }
 const double vb_safe = std::max(vb, 1e-300);
 double lp = 0.0;
 const arma::uword m = b.n_elem;
 for (arma::uword j = 0; j < m; ++j) {
  if (d(j) <= 0) continue;
  const double log_q = (s + 1.0) * log_h(j);
  const double q = std::exp(log_q);
  if (!std::isfinite(q) || q <= 0.0) {
   return -std::numeric_limits<double>::infinity();
  }
  lp += -0.5 * (log_q + b(j) * b(j) / (vb_safe * q));
 }
 return lp;
}

inline bool update_maf_effect_s_bayesc(
  double& maf_effect_s_current,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  double vb,
  const arma::rowvec& log_h,
  double prior_lower,
  double prior_upper,
  double proposal_sd,
  std::mt19937& gen
) {
 std::normal_distribution<double> proposal(0.0, proposal_sd);
 const double maf_effect_s_prop = maf_effect_s_current + proposal(gen);
 if (maf_effect_s_prop < prior_lower || maf_effect_s_prop > prior_upper ||
     !std::isfinite(maf_effect_s_prop)) {
  return false;
 }

 const double lp_current = logpost_maf_effect_s_bayesc(
  maf_effect_s_current, b, d, vb, log_h, prior_lower, prior_upper
 );
 const double lp_prop = logpost_maf_effect_s_bayesc(
  maf_effect_s_prop, b, d, vb, log_h, prior_lower, prior_upper
 );
 const double log_alpha = lp_prop - lp_current;
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 if (std::log(std::max(runif(gen), 1e-300)) < log_alpha) {
  maf_effect_s_current = maf_effect_s_prop;
  return true;
 }
 return false;
}


inline void samplePi_ST(
  const arma::Row<int>& d,
  std::vector<double>& pi,
  double pi_prior_a,
  double pi_prior_b,
  std::mt19937& gen
) {
 // pi[1] is inclusion probability; pi[0] is exclusion probability.
 // Prior: pi[1] ~ Beta(pi_prior_a, pi_prior_b).
 double c1 = pi_prior_a;
 double c0 = pi_prior_b;

 for (arma::uword i = 0; i < d.n_elem; ++i) {
  if (d(i) > 0) c1 += 1.0;
  else c0 += 1.0;
 }

 std::gamma_distribution<double> rg0(c0, 1.0);
 std::gamma_distribution<double> rg1(c1, 1.0);

 const double g0 = std::max(rg0(gen), 1e-300);
 const double g1 = std::max(rg1(gen), 1e-300);
 const double s = g0 + g1;

 pi[0] = g0 / s;
 pi[1] = g1 / s;
}

// -----------------------------------------------------------------------------
// Linkage-equilibrium and linkage-disequilibrium variance components
// -----------------------------------------------------------------------------

inline double computeLE_ST_csr(
  int m,
  const arma::rowvec& b,
  const arma::rowvec& ww,
  int n
) {
 double vle = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double bi = b(iu);

  if (bi != 0.0) {
   vle += ww(iu) * bi * bi;
  }
 }

 return vle / static_cast<double>(n);
}

template <class OpT>
inline void set_marker_effect_st_csr(
  int i,
  double b_new,
  int d_new,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const OpT& op
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  op.apply_difference(i, diff, r);
 }

 b(iu) = b_new;
 d(iu) = d_new;
}

inline double residual_sse_st_csr(
  int m,
  const arma::rowvec& b,
  const arma::rowvec& wy,
  const arma::rowvec& r,
  double yy
) {
 double b_dot_r_plus_wy = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  b_dot_r_plus_wy += b(iu) * (r(iu) + wy(iu));
 }

 return yy - b_dot_r_plus_wy;
}

inline int count_excluded_ld_friends(
  int i,
  const arma::Row<int>& d,
  const LDLDFriends& friends
) {
 int n = 0;
 const uint64_t start = friends.ptr[static_cast<std::size_t>(i)];
 const uint64_t end   = friends.ptr[static_cast<std::size_t>(i + 1)];

 for (uint64_t p = start; p < end; ++p) {
  const int j = friends.idx[static_cast<std::size_t>(p)];
  if (d(static_cast<arma::uword>(j)) == 0) ++n;
 }

 return n;
}

inline int collect_ld_swap_candidates(
  int m,
  const arma::Row<int>& d,
  const LDLDFriends& friends,
  std::vector<int>& candidates,
  std::vector<int>& n_excluded
) {
 candidates.clear();
 n_excluded.clear();

 for (int i = 0; i < m; ++i) {
  if (d(static_cast<arma::uword>(i)) <= 0) continue;

  const int nf = count_excluded_ld_friends(i, d, friends);
  if (nf > 0) {
   candidates.push_back(i);
   n_excluded.push_back(nf);
  }
 }

 return static_cast<int>(candidates.size());
}

template <class OpT>
inline bool attempt_ld_swap_st_csr(
  int m,
  double vei,
  double vb,
  double yy,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  const arma::rowvec& prior_scale,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const OpT& op,
  const LDLDFriends& friends,
  std::mt19937& gen
) {
 if (!std::isfinite(vei) || vei <= 0.0) return false;

 std::vector<int> candidates;
 std::vector<int> n_excluded;
 const int n_candidates = collect_ld_swap_candidates(m, d, friends, candidates, n_excluded);
 if (n_candidates <= 0) return false;

 std::uniform_int_distribution<int> pick_candidate(0, n_candidates - 1);
 const int cand_pos = pick_candidate(gen);
 const int j = candidates[static_cast<std::size_t>(cand_pos)];
 const int n_forward_friends = n_excluded[static_cast<std::size_t>(cand_pos)];
 if (n_forward_friends <= 0) return false;

 std::uniform_int_distribution<int> pick_friend(0, n_forward_friends - 1);
 const int friend_pos = pick_friend(gen);

 int k = -1;
 int seen = 0;
 const uint64_t start = friends.ptr[static_cast<std::size_t>(j)];
 const uint64_t end   = friends.ptr[static_cast<std::size_t>(j + 1)];

 for (uint64_t p = start; p < end; ++p) {
  const int jj = friends.idx[static_cast<std::size_t>(p)];
  if (d(static_cast<arma::uword>(jj)) != 0) continue;

  if (seen == friend_pos) {
   k = jj;
   break;
  }
  ++seen;
 }

 if (k < 0) return false;

 const arma::uword ju = static_cast<arma::uword>(j);
 const arma::uword ku = static_cast<arma::uword>(k);
 const double b_j_old = b(ju);
 const double b_k_old = b(ku);
 const int d_j_old = d(ju);
 const int d_k_old = d(ku);
 if (d_j_old <= 0 || d_k_old != 0 || b_j_old == 0.0) return false;
 const double vb_j = std::max(vb * prior_scale(ju), 1e-300);
 const double vb_k = std::max(vb * prior_scale(ku), 1e-300);

 const double sse_old = residual_sse_st_csr(m, b, wy, r, yy);
 if (!std::isfinite(sse_old)) return false;

 const arma::rowvec r_old = r;

 set_marker_effect_st_csr(j, 0.0, 0, ww, r, b, d, op);
 set_marker_effect_st_csr(k, b_j_old, 1, ww, r, b, d, op);

 const double sse_new = residual_sse_st_csr(m, b, wy, r, yy);
 bool accept = false;

 if (std::isfinite(sse_new)) {
  std::vector<int> reverse_candidates;
  std::vector<int> reverse_n_excluded;
  const int n_reverse_candidates =
   collect_ld_swap_candidates(m, d, friends, reverse_candidates, reverse_n_excluded);

  int n_reverse_friends = 0;
  for (std::size_t pos = 0; pos < reverse_candidates.size(); ++pos) {
   if (reverse_candidates[pos] == k) {
    n_reverse_friends = reverse_n_excluded[pos];
    break;
   }
  }

  if (n_reverse_candidates > 0 && n_reverse_friends > 0) {
   const double log_q_forward =
    -std::log(static_cast<double>(n_candidates)) -
    std::log(static_cast<double>(n_forward_friends));
   const double log_q_reverse =
    -std::log(static_cast<double>(n_reverse_candidates)) -
    std::log(static_cast<double>(n_reverse_friends));
   // With constant BayesC prior variance, active/null LD-swap relocations
   // cancel the prior density. Fixed maf_effect_s makes the included-marker
   // prior variance marker-specific, so the N(0, vb * prior_scale_j) density
   // must enter the MH ratio for the moved effect.
   const double log_prior_ratio =
    -0.5 * (std::log(vb_k / vb_j) +
            b_j_old * b_j_old * (1.0 / vb_k - 1.0 / vb_j));
   const double log_alpha =
    -0.5 * (sse_new - sse_old) / vei +
    log_prior_ratio + log_q_reverse - log_q_forward;

   std::uniform_real_distribution<double> runif(0.0, 1.0);
   accept = std::log(std::max(runif(gen), 1e-300)) < log_alpha;
  }
 }

 if (!accept) {
  r = r_old;
  b(ju) = b_j_old;
  b(ku) = b_k_old;
  d(ju) = d_j_old;
  d(ku) = d_k_old;
 }

 return accept;
}

template <class OpT>
inline bool attempt_ld_swap_st_csr_unscaled(
  int m,
  double vei,
  double yy,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const OpT& op,
  const LDLDFriends& friends,
  std::mt19937& gen
) {
 if (!std::isfinite(vei) || vei <= 0.0) return false;

 std::vector<int> candidates;
 std::vector<int> n_excluded;
 const int n_candidates = collect_ld_swap_candidates(m, d, friends, candidates, n_excluded);
 if (n_candidates <= 0) return false;

 std::uniform_int_distribution<int> pick_candidate(0, n_candidates - 1);
 const int cand_pos = pick_candidate(gen);
 const int j = candidates[static_cast<std::size_t>(cand_pos)];
 const int n_forward_friends = n_excluded[static_cast<std::size_t>(cand_pos)];
 if (n_forward_friends <= 0) return false;

 std::uniform_int_distribution<int> pick_friend(0, n_forward_friends - 1);
 const int friend_pos = pick_friend(gen);

 int k = -1;
 int seen = 0;
 const uint64_t start = friends.ptr[static_cast<std::size_t>(j)];
 const uint64_t end   = friends.ptr[static_cast<std::size_t>(j + 1)];

 for (uint64_t p = start; p < end; ++p) {
  const int jj = friends.idx[static_cast<std::size_t>(p)];
  if (d(static_cast<arma::uword>(jj)) != 0) continue;

  if (seen == friend_pos) {
   k = jj;
   break;
  }
  ++seen;
 }

 if (k < 0) return false;

 const arma::uword ju = static_cast<arma::uword>(j);
 const arma::uword ku = static_cast<arma::uword>(k);
 const double b_j_old = b(ju);
 const double b_k_old = b(ku);
 const int d_j_old = d(ju);
 const int d_k_old = d(ku);
 if (d_j_old <= 0 || d_k_old != 0 || b_j_old == 0.0) return false;

 const double sse_old = residual_sse_st_csr(m, b, wy, r, yy);
 if (!std::isfinite(sse_old)) return false;

 const arma::rowvec r_old = r;

 set_marker_effect_st_csr(j, 0.0, 0, ww, r, b, d, op);
 set_marker_effect_st_csr(k, b_j_old, 1, ww, r, b, d, op);

 const double sse_new = residual_sse_st_csr(m, b, wy, r, yy);
 bool accept = false;

 if (std::isfinite(sse_new)) {
  std::vector<int> reverse_candidates;
  std::vector<int> reverse_n_excluded;
  const int n_reverse_candidates =
   collect_ld_swap_candidates(m, d, friends, reverse_candidates, reverse_n_excluded);

  int n_reverse_friends = 0;
  for (std::size_t pos = 0; pos < reverse_candidates.size(); ++pos) {
   if (reverse_candidates[pos] == k) {
    n_reverse_friends = reverse_n_excluded[pos];
    break;
   }
  }

  if (n_reverse_candidates > 0 && n_reverse_friends > 0) {
   const double log_q_forward =
    -std::log(static_cast<double>(n_candidates)) -
    std::log(static_cast<double>(n_forward_friends));
   const double log_q_reverse =
    -std::log(static_cast<double>(n_reverse_candidates)) -
    std::log(static_cast<double>(n_reverse_friends));
   const double log_alpha =
    -0.5 * (sse_new - sse_old) / vei +
    log_q_reverse - log_q_forward;

   std::uniform_real_distribution<double> runif(0.0, 1.0);
   accept = std::log(std::max(runif(gen), 1e-300)) < log_alpha;
  }
 }

 if (!accept) {
  r = r_old;
  b(ju) = b_j_old;
  b(ku) = b_k_old;
  d(ju) = d_j_old;
  d(ku) = d_k_old;
 }

 return accept;
}

// This is the sole typed-result-to-stblr_raw_v1 conversion for ordinary CSR
// BayesC. It is binding-specific and executes only after the core returns.
Rcpp::List stblr_csr_bayesc_result_to_raw(
  const sblr::core::CsrBayesCResult& result,
  const CsrBayesCRawConversionContext& context
) {
 const int m = context.marker_count;
 const int nt = context.trait_count;
 const int nit = context.nit;
 const int nburn = context.nburn;
 const int nthin = context.nthin;
 const int ncores = context.ncores;
 const int nchains = context.nchains;
 const bool keep_chains = context.keep_chains;
 const double pi_prior_a = context.pi_prior_a;
 const double pi_prior_b = context.pi_prior_b;
 const bool updateLDswap = context.update_ld_swap;
 const bool use_maf_effect_s_prior_scale =
  context.use_fixed_maf_effect_scale;
 const bool estimate_maf_effect_s = context.estimate_maf_effect_s;
 if (context.sample_size == nullptr) {
  throw std::runtime_error(
   "stblr_csr_bayesc_result_to_raw: sample-size metadata is missing."
  );
 }
 const std::vector<int>& n = *context.sample_size;
 if (context.convergence_markers == nullptr) {
  throw std::runtime_error(
   "stblr_csr_bayesc_result_to_raw: convergence-marker metadata is missing."
  );
 }
 const std::vector<int>& convergence_markers = *context.convergence_markers;
 const BlrPhase3ExecutionContract legacy_execution;
 const BlrPhase3ExecutionContract& execution_contract =
  context.execution_contract == nullptr
   ? legacy_execution : *context.execution_contract;

#ifdef _OPENMP
 const int nthreads = stblr_num_threads_for_tasks(
  ncores, stblr_num_chain_tasks(nt, nchains)
 );
 for (int task = 0; task < stblr_num_chain_tasks(nt, nchains); ++task) {
  const sblr::core::CsrBayesCChainResult& chain_result =
   result.chains[static_cast<std::size_t>(task)];
 }
#endif

 const int n_trace = nit + nburn;
 const bool return_chain_summaries = (nchains > 1) || keep_chains;
 auto marker_matrix = [&](const arma::mat& values) {
  Rcpp::NumericMatrix output(m, nt);
  for (int trait = 0; trait < nt; ++trait) {
   for (int marker = 0; marker < m; ++marker) {
    output(marker, trait) = values(trait, marker);
   }
  }
  return output;
 };
 auto trace_matrix = [&](const arma::mat& values) {
  Rcpp::NumericMatrix output(n_trace, nt);
  for (int trait = 0; trait < nt; ++trait) {
   for (int iteration = 0; iteration < n_trace; ++iteration) {
    output(iteration, trait) = values(trait, iteration);
   }
  }
  return output;
 };
 auto diagonal_matrix = [&](const arma::vec& values) {
  Rcpp::NumericMatrix output(nt, nt);
  for (int trait = 0; trait < nt; ++trait) output(trait, trait) = values(trait);
  return output;
 };

 Rcpp::NumericMatrix pi_final(nt, 2);
 Rcpp::NumericMatrix pi_mean(nt, 2);
 Rcpp::NumericVector selection_mean(nt);
 Rcpp::NumericVector maf_effect_sd(nt);
 Rcpp::NumericVector selection_min(nt);
 Rcpp::NumericVector selection_max(nt);
 Rcpp::NumericVector selection_acceptance(nt);
 for (int trait = 0; trait < nt; ++trait) {
  const arma::uword trait_u = static_cast<arma::uword>(trait);
  pi_final(trait, 0) = 1.0 - result.final_inclusion_probability(trait_u);
  pi_final(trait, 1) = result.final_inclusion_probability(trait_u);
  double mean_pi = 0.0;
  int pi_samples = 0;
  for (int iteration = nburn; iteration < n_trace; ++iteration) {
   mean_pi += result.inclusion_trace(trait_u, iteration);
   ++pi_samples;
  }
  if (pi_samples > 0) mean_pi /= static_cast<double>(pi_samples);
  else mean_pi = result.final_inclusion_probability(trait_u);
  pi_mean(trait, 0) = 1.0 - mean_pi;
  pi_mean(trait, 1) = mean_pi;
  selection_acceptance[trait] = result.maf_effect_s_attempted(trait_u) > 0.0
   ? result.maf_effect_s_accepted(trait_u) /
     result.maf_effect_s_attempted(trait_u) : 0.0;
  if (estimate_maf_effect_s) {
   double mean = 0.0;
   double minimum = std::numeric_limits<double>::infinity();
   double maximum = -std::numeric_limits<double>::infinity();
   int count = 0;
   for (int iteration = nburn; iteration < n_trace; ++iteration) {
    const double value = result.maf_effect_s_trace(trait_u, iteration);
    mean += value;
    minimum = std::min(minimum, value);
    maximum = std::max(maximum, value);
    ++count;
   }
   if (count > 0) mean /= static_cast<double>(count);
   double sum_squares = 0.0;
   if (count > 1) {
    for (int iteration = nburn; iteration < n_trace; ++iteration) {
     const double difference =
      result.maf_effect_s_trace(trait_u, iteration) - mean;
     sum_squares += difference * difference;
    }
    maf_effect_sd[trait] = std::sqrt(
     sum_squares / static_cast<double>(count - 1)
    );
   } else maf_effect_sd[trait] = NA_REAL;
   selection_mean[trait] = mean;
   selection_min[trait] = count > 0 ? minimum : NA_REAL;
   selection_max[trait] = count > 0 ? maximum : NA_REAL;
  } else {
   selection_mean[trait] = NA_REAL;
   maf_effect_sd[trait] = NA_REAL;
   selection_min[trait] = NA_REAL;
   selection_max[trait] = NA_REAL;
  }
 }

 Rcpp::NumericMatrix ld_swap(nt, 3);
 Rcpp::NumericVector nsamples(nt);
 Rcpp::IntegerVector n_used(nt);
 Rcpp::NumericVector seconds_mean(nt);
 Rcpp::NumericVector seconds_max(nt);
 for (int trait = 0; trait < nt; ++trait) {
  const arma::uword trait_u = static_cast<arma::uword>(trait);
  const double attempted = result.ld_swap_attempted(trait_u);
  const double accepted = result.ld_swap_accepted(trait_u);
  ld_swap(trait, 0) = attempted;
  ld_swap(trait, 1) = accepted;
  ld_swap(trait, 2) = attempted > 0.0 ? accepted / attempted : 0.0;
  nsamples[trait] = result.retained_samples(trait_u);
  n_used[trait] = n[static_cast<std::size_t>(trait)];
  double seconds_sum = 0.0;
  double maximum = 0.0;
  for (int chain = 0; chain < nchains; ++chain) {
   const double seconds = result.chains[
    static_cast<std::size_t>(trait * nchains + chain)
   ].seconds;
   seconds_sum += seconds;
   maximum = std::max(maximum, seconds);
  }
  seconds_mean[trait] = seconds_sum / static_cast<double>(nchains);
  seconds_max[trait] = maximum;
 }

 Rcpp::List marker = Rcpp::List::create(
  Rcpp::Named("bm") = marker_matrix(result.marker_mean),
  Rcpp::Named("dm") = marker_matrix(result.marker_pip),
  Rcpp::Named("wy") = marker_matrix(result.marker_score),
  Rcpp::Named("r") = marker_matrix(result.final_residual),
  Rcpp::Named("b") = marker_matrix(result.final_effect),
  Rcpp::Named("state") = marker_matrix(result.final_state)
 );
 if (return_chain_summaries) {
  marker["bm_sd"] = marker_matrix(result.marker_mean_sd);
  marker["bm_min"] = marker_matrix(result.marker_mean_min);
  marker["bm_max"] = marker_matrix(result.marker_mean_max);
  marker["dm_sd"] = marker_matrix(result.marker_pip_sd);
  marker["dm_min"] = marker_matrix(result.marker_pip_min);
  marker["dm_max"] = marker_matrix(result.marker_pip_max);
 }
 Rcpp::List trace = Rcpp::List::create(
  Rcpp::Named("vbs") = trace_matrix(result.marker_variance_trace),
  Rcpp::Named("vgs") = trace_matrix(result.genetic_variance_trace),
  Rcpp::Named("ves") = trace_matrix(result.residual_variance_trace),
  Rcpp::Named("vle") = trace_matrix(result.le_variance_trace),
  Rcpp::Named("vld") = trace_matrix(result.ld_variance_trace),
  Rcpp::Named("pis") = trace_matrix(result.inclusion_trace)
 );
 Rcpp::List variance = Rcpp::List::create(
  Rcpp::Named("covb") = diagonal_matrix(result.final_marker_variance),
  Rcpp::Named("covg") = diagonal_matrix(result.final_genetic_variance),
  Rcpp::Named("cove") = diagonal_matrix(result.final_residual_variance),
  Rcpp::Named("vb") = diagonal_matrix(result.final_marker_variance),
  Rcpp::Named("vg") = diagonal_matrix(result.final_genetic_variance),
  Rcpp::Named("ve") = diagonal_matrix(result.final_residual_variance)
 );
 Rcpp::List diagnostics = Rcpp::List::create(
  Rcpp::Named("nsamples") = nsamples,
  Rcpp::Named("n_used") = n_used,
  Rcpp::Named("log_cpo") = Rcpp::NumericVector(nt),
  Rcpp::Named("mean_log_cpo") = Rcpp::NumericVector(nt),
  Rcpp::Named("seconds_mean") = seconds_mean,
  Rcpp::Named("seconds_max") = seconds_max,
  Rcpp::Named("ld_swap") = updateLDswap ? Rcpp::wrap(ld_swap) : R_NilValue
 );
 if (execution_contract.active()) {
  diagnostics.push_back([&]() {
   std::vector<int> worker_ids(result.chains.size(), 0);
   std::vector<int> team_sizes(result.chains.size(), 1);
   for (std::size_t task = 0; task < result.chains.size(); ++task) {
    worker_ids[task] = result.chains[task].thread_used;
    team_sizes[task] = result.chains[task].team_size_used;
   }
   return blr_phase3_worker_diagnostics(
    execution_contract, ncores,
    stblr_num_threads_for_tasks(ncores, static_cast<int>(result.chains.size())),
    worker_ids, team_sizes
   );
  }(), "workers");
 }

 Rcpp::List chains = R_NilValue;
 if (keep_chains) {
  chains = Rcpp::List(nt);
  Rcpp::CharacterVector trait_names(nt);
  for (int trait = 0; trait < nt; ++trait) {
   trait_names[trait] = "trait" + std::to_string(trait + 1);
   Rcpp::List trait_chains(nchains);
   Rcpp::CharacterVector chain_names(nchains);
   for (int chain = 0; chain < nchains; ++chain) {
    chain_names[chain] = "chain" + std::to_string(chain + 1);
    const sblr::core::CsrBayesCChainResult& current = result.chains[
     static_cast<std::size_t>(trait * nchains + chain)
    ];
    Rcpp::NumericVector chain_bm(m);
    Rcpp::NumericVector chain_dm(m);
    Rcpp::NumericVector chain_state(m);
    for (int marker_index = 0; marker_index < m; ++marker_index) {
     chain_bm[marker_index] = current.marker_mean(marker_index);
     chain_dm[marker_index] = current.marker_pip(marker_index);
     chain_state[marker_index] = current.final_state(marker_index);
    }
    Rcpp::NumericMatrix chain_ld(1, 3);
    chain_ld(0, 0) = current.ld_swap_attempted;
    chain_ld(0, 1) = current.ld_swap_accepted;
    chain_ld(0, 2) = current.ld_swap_attempted > 0.0
     ? current.ld_swap_accepted / current.ld_swap_attempted : 0.0;
    Rcpp::List chain_selection = Rcpp::List::create(
     Rcpp::Named("trace") = R_NilValue,
     Rcpp::Named("acceptance") = R_NilValue
    );
    if (estimate_maf_effect_s) {
     Rcpp::NumericVector chain_selection_trace(n_trace);
     for (int iteration = 0; iteration < n_trace; ++iteration) {
      chain_selection_trace[iteration] = current.maf_effect_s_trace(iteration);
     }
     chain_selection["trace"] = chain_selection_trace;
     chain_selection["acceptance"] = current.maf_effect_s_attempted > 0.0
      ? current.maf_effect_s_accepted / current.maf_effect_s_attempted : 0.0;
    }
    trait_chains[chain] = Rcpp::List::create(
     Rcpp::Named("marker") = Rcpp::List::create(
      Rcpp::Named("bm") = chain_bm,
      Rcpp::Named("dm") = chain_dm,
      Rcpp::Named("state") = chain_state
     ),
     Rcpp::Named("trace") = Rcpp::List::create(
      Rcpp::Named("vbs") = current.marker_variance_trace,
      Rcpp::Named("vgs") = current.genetic_variance_trace,
      Rcpp::Named("ves") = current.residual_variance_trace,
      Rcpp::Named("vle") = current.le_variance_trace,
      Rcpp::Named("vld") = current.ld_variance_trace,
      Rcpp::Named("pis") = current.inclusion_trace
     ),
     Rcpp::Named("pi") = Rcpp::List::create(
      Rcpp::Named("final") = Rcpp::NumericVector::create(
       1.0 - current.final_inclusion_probability,
       current.final_inclusion_probability
      ),
      Rcpp::Named("mean") = R_NilValue
     ),
     Rcpp::Named("selection") = chain_selection,
     Rcpp::Named("convergence_trace") = Rcpp::List::create(
      Rcpp::Named("b") = current.convergence_b,
      Rcpp::Named("d") = current.convergence_d,
      Rcpp::Named("component") = R_NilValue,
      Rcpp::Named("marker_index") = Rcpp::wrap(convergence_markers)),
     Rcpp::Named("diagnostics") = Rcpp::List::create(
      Rcpp::Named("ld_swap") = updateLDswap
       ? Rcpp::wrap(chain_ld) : R_NilValue
     )
    );
   }
   trait_chains.attr("names") = chain_names;
   chains[trait] = trait_chains;
  }
  chains.attr("names") = trait_names;
 }

 Rcpp::List selection = Rcpp::List::create(
  Rcpp::Named("enabled") = estimate_maf_effect_s || use_maf_effect_s_prior_scale,
  Rcpp::Named("fixed") = use_maf_effect_s_prior_scale,
  Rcpp::Named("scale") = "standardized_genotype_effect",
  Rcpp::Named("trace") = estimate_maf_effect_s
   ? Rcpp::wrap(trace_matrix(result.maf_effect_s_trace)) : R_NilValue,
  Rcpp::Named("mean") = estimate_maf_effect_s
   ? Rcpp::wrap(selection_mean) : R_NilValue,
  Rcpp::Named("sd") = estimate_maf_effect_s
   ? Rcpp::wrap(maf_effect_sd) : R_NilValue,
  Rcpp::Named("min") = estimate_maf_effect_s
   ? Rcpp::wrap(selection_min) : R_NilValue,
  Rcpp::Named("max") = estimate_maf_effect_s
   ? Rcpp::wrap(selection_max) : R_NilValue,
  Rcpp::Named("acceptance") = estimate_maf_effect_s
   ? Rcpp::wrap(selection_acceptance) : R_NilValue
 );

 Rcpp::List raw = Rcpp::List::create(
  Rcpp::Named("schema") = Rcpp::List::create(
   Rcpp::Named("class") = "stblr_raw",
   Rcpp::Named("version") = 1
  ),
  Rcpp::Named("meta") = Rcpp::List::create(
   Rcpp::Named("model") = "bayesc",
   Rcpp::Named("backend") = "csr_bayesc",
   Rcpp::Named("data_level") = "summary",
   Rcpp::Named("prior_type") = "global",
   Rcpp::Named("m") = m,
   Rcpp::Named("nt") = nt,
   Rcpp::Named("n_trace") = n_trace,
   Rcpp::Named("nit") = nit,
   Rcpp::Named("nburn") = nburn,
   Rcpp::Named("nthin") = nthin,
   Rcpp::Named("nchains") = nchains,
   Rcpp::Named("keep_chains") = keep_chains,
   Rcpp::Named("n_components") = 2,
   Rcpp::Named("n_annotations") = 0,
   Rcpp::Named("n_groups") = 0
  ),
  Rcpp::Named("marker") = marker,
  Rcpp::Named("trace") = trace,
  Rcpp::Named("variance") = variance,
  Rcpp::Named("pi") = Rcpp::List::create(
   Rcpp::Named("final") = pi_final,
   Rcpp::Named("mean") = pi_mean,
   Rcpp::Named("names") = Rcpp::CharacterVector::create("pi0", "pi1")
  ),
  Rcpp::Named("diagnostics") = diagnostics,
  Rcpp::Named("chains") = chains,
  Rcpp::Named("prior") = Rcpp::List::create(),
  Rcpp::Named("group") = Rcpp::List::create(),
  Rcpp::Named("annotation") = Rcpp::List::create(),
  Rcpp::Named("component") = Rcpp::List::create(),
  Rcpp::Named("selection") = selection
 );
 raw.attr("class") = Rcpp::CharacterVector::create(
  "stblr_raw_v1", "stblr_raw", "list"
 );
 return raw;
}

// The ordinary CSR native adapter constructs borrowed typed views, invokes
// the one canonical core, and delegates all R object construction above.
static Rcpp::List stblr_csr_bayesc_run_canonical(
  int m,
  int nt,
  const std::string& ld_prefix,
  const CsrOperator& op,
  const LDLDFriends& ld_swap_friends,
  const arma::mat& wy_mat,
  const arma::mat& b_mat,
  const std::vector<std::vector<double>>& d_init,
  bool use_d_init,
  const std::vector<std::vector<double>>& r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  const arma::vec& yy_vec,
  const arma::mat& B,
  const arma::mat& E,
  const arma::mat& ssb_prior_mat,
  const arma::mat& sse_prior_mat,
  const std::vector<double>& pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  double adjE,
  const std::vector<int>& n,
  int nit,
  int nburn,
  int nthin,
  double pi_prior_a,
  double pi_prior_b,
  int ncores,
  int seed,
  int nchains,
  bool keep_chains,
  const std::vector<int>& chain_seeds,
  bool updateLDswap,
  double ld_swap_prob,
  int ld_swap_moves,
  bool use_maf_effect_s_prior_scale,
  const arma::rowvec& prior_scale,
  bool estimate_maf_effect_s,
  double maf_effect_s_init,
  const Rcpp::NumericVector& maf_effect_s_prior,
  double maf_effect_s_proposal_sd,
  const arma::rowvec& maf_effect_s_log_h_row,
  const std::vector<int>& convergence_markers,
  bool convergence_b,
  bool convergence_d,
  const std::vector<int>& order,
  const BlrPhase3ExecutionContract& execution_contract
) {
 sblr::core::ResolvedSpec specification;
 specification.data.marker_count = static_cast<std::size_t>(m);
 specification.data.trait_count = static_cast<std::size_t>(nt);
 specification.data.marker_ids.reserve(static_cast<std::size_t>(m));
 specification.data.trait_ids.reserve(static_cast<std::size_t>(nt));
 for (int marker = 0; marker < m; ++marker) {
  specification.data.marker_ids.push_back(
   "marker_" + std::to_string(marker + 1)
  );
 }
 for (int trait = 0; trait < nt; ++trait) {
  specification.data.trait_ids.push_back(
   "trait_" + std::to_string(trait + 1)
  );
 }
 specification.data.sample_size = n;
 specification.data.csr.resource_id = ld_prefix;
 specification.data.csr.marker_count = static_cast<std::size_t>(m);
 specification.mcmc.nit = nit;
 specification.mcmc.nburn = nburn;
 specification.mcmc.nthin = nthin;
 specification.mcmc.nchains = nchains;
 specification.mcmc.ncores = ncores;
 specification.mcmc.seed = seed;
 specification.mcmc.has_explicit_chain_seeds = !chain_seeds.empty();
 specification.mcmc.chain_seeds = chain_seeds;
 specification.output.keep_chain_summaries = keep_chains;

 sblr::core::CsrBayesCExecutionInput input;
 input.specification = std::move(specification);
 input.data.marker_count = static_cast<std::size_t>(m);
 input.data.trait_count = static_cast<std::size_t>(nt);
 input.data.ld = op.view();
 input.data.wy = &wy_mat;
 input.data.yy = &yy_vec;
 input.data.sample_size = &n;
 input.priors.marker_variance = &B;
 input.priors.residual_variance = &E;
 input.priors.marker_scale_prior = &ssb_prior_mat;
 input.priors.residual_scale_prior = &sse_prior_mat;
 input.priors.inclusion_probability = pi;
 input.priors.marker_degrees_freedom = nub;
 input.priors.residual_degrees_freedom = nue;
 input.priors.inclusion_prior_active = pi_prior_a;
 input.priors.inclusion_prior_null = pi_prior_b;
 input.initial.effects = &b_mat;
 input.initial.inclusion = &d_init;
 input.initial.residual = &r_init;
 input.initial.use_inclusion = use_d_init;
 input.initial.use_residual = use_r_init;
 input.controls.nit = nit;
 input.controls.nburn = nburn;
 input.controls.nthin = nthin;
 input.controls.nchains = nchains;
 input.controls.ncores = ncores;
 input.controls.seed = seed;
 input.controls.chain_seeds = chain_seeds;
 input.controls.seed_contract_version = execution_contract.seed_contract_version;
 input.controls.retention_contract_version =
  execution_contract.retention_contract_version;
 input.controls.scheduler_version = execution_contract.scheduler_version;
 input.controls.task_seeds = execution_contract.task_seeds;
 input.controls.retained_transition_indices =
  execution_contract.retained_transition_indices;
 input.controls.keep_chains = keep_chains;
 input.controls.update_marker_variance = updateB;
 input.controls.update_residual_variance = updateE;
 input.controls.update_inclusion_probability = updatePi;
 input.controls.rebuild_residual_before_update = rebuild_r_before_updateE;
 input.controls.residual_adjustment = adjE;
 input.controls.update_ld_swap = updateLDswap;
 input.controls.ld_swap_probability = ld_swap_prob;
 input.controls.ld_swap_moves = ld_swap_moves;
 input.controls.use_fixed_maf_effect_scale = use_maf_effect_s_prior_scale;
 input.controls.fixed_maf_effect_scale =
  use_maf_effect_s_prior_scale ? &prior_scale : nullptr;
 input.controls.estimate_maf_effect_s = estimate_maf_effect_s;
 input.controls.maf_effect_s_initial = maf_effect_s_init;
 input.controls.maf_effect_s_prior_lower = maf_effect_s_prior[0];
 input.controls.maf_effect_s_prior_upper = maf_effect_s_prior[1];
 input.controls.maf_effect_s_proposal_sd = maf_effect_s_proposal_sd;
 input.controls.maf_effect_s_log_h =
  estimate_maf_effect_s ? &maf_effect_s_log_h_row : nullptr;
 input.controls.convergence_markers = convergence_markers;
 input.controls.convergence_b = convergence_b;
 input.controls.convergence_d = convergence_d;
 input.output.keep_chains = keep_chains;
 input.ld_friends.row_ptr = ld_swap_friends.ptr.empty()
  ? nullptr : ld_swap_friends.ptr.data();
 input.ld_friends.row_ptr_size = ld_swap_friends.ptr.size();
 input.ld_friends.index = ld_swap_friends.idx.empty()
  ? nullptr : ld_swap_friends.idx.data();
 input.ld_friends.friend_count = ld_swap_friends.idx.size();
 input.marker_order = &order;

 const sblr::core::CsrBayesCResult result =
  sblr::core::run_csr_bayesc(input);
 const CsrBayesCRawConversionContext conversion = {
  m, nt, nit, nburn, nthin, ncores, nchains, keep_chains,
  pi_prior_a, pi_prior_b, updateLDswap, use_maf_effect_s_prior_scale,
  estimate_maf_effect_s, &n, &convergence_markers, &execution_contract
 };
 return stblr_csr_bayesc_result_to_raw(result, conversion);
}

// -----------------------------------------------------------------------------
// Main implementation: parallel single-trait BayesC over traits
// -----------------------------------------------------------------------------

template <class OperatorFactory>
Rcpp::List stblr_cpg_omp_csr_impl(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> d_init,
  bool use_d_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  std::string ld_prefix,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<double> pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  double adjE,
  std::vector<int> n,
  int nit,
  int nburn,
  int nthin,
  double pi_prior_a,
  double pi_prior_b,
  int ncores,
  int seed,
  int nchains,
  bool keep_chains,
  std::vector<int> chain_seeds,
  bool updateLDswap,
  double ld_swap_prob,
  double ld_swap_r2,
  int ld_swap_max_friends,
  int ld_swap_moves,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale,
  bool estimate_maf_effect_s,
  double maf_effect_s_init,
  Rcpp::NumericVector maf_effect_s_prior,
  double maf_effect_s_proposal_sd,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h,
  std::vector<int> convergence_markers,
  bool convergence_b,
  bool convergence_d,
  int low_rank_residual_rebuild_every,
  CsrBayesCPolicyFactory* policy_factory,
  const BlrPhase3ExecutionContract& execution_contract,
  OperatorFactory make_operator
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr: nt must be positive.");
 }

 const int m = static_cast<int>(wy[0].size());

 if (m <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr: m must be positive.");
 }

 if (nit <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr: nit must be positive.");
 }

 if (nburn < 0) {
  throw std::runtime_error("stblr_cpg_omp_csr: nburn must be non-negative.");
 }

 if (nthin <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr: nthin must be positive.");
 }

 if (nchains <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr: nchains must be positive.");
 }

 if (!chain_seeds.empty() && static_cast<int>(chain_seeds.size()) != nchains) {
  throw std::runtime_error("stblr_cpg_omp_csr: chain_seeds must have length nchains.");
 }

 if (!std::isfinite(ld_swap_prob) || ld_swap_prob < 0.0 || ld_swap_prob > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr: ld_swap_prob must be in [0, 1].");
 }

 if (!std::isfinite(ld_swap_r2) || ld_swap_r2 < 0.0 || ld_swap_r2 > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr: ld_swap_r2 must be in [0, 1].");
 }

 if (ld_swap_max_friends <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr: ld_swap_max_friends must be positive.");
 }

 if (ld_swap_moves < 0) {
  throw std::runtime_error("stblr_cpg_omp_csr: ld_swap_moves must be non-negative.");
 }

 bool use_maf_effect_s_prior_scale = maf_effect_s_prior_scale.isNotNull();
 Rcpp::NumericVector maf_effect_s_prior_scale_vec;
 if (use_maf_effect_s_prior_scale) {
  maf_effect_s_prior_scale_vec = Rcpp::NumericVector(maf_effect_s_prior_scale);
  use_maf_effect_s_prior_scale = maf_effect_s_prior_scale_vec.size() > 0;
 }

 if (use_maf_effect_s_prior_scale &&
     static_cast<int>(maf_effect_s_prior_scale_vec.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_prior_scale must have length m.");
 }

 if (estimate_maf_effect_s && use_maf_effect_s_prior_scale) {
  throw std::runtime_error(
    "stblr_cpg_omp_csr: fixed maf_effect_s_prior_scale and estimate_maf_effect_s cannot both be used."
  );
 }

 if (estimate_maf_effect_s) {
  if (!std::isfinite(maf_effect_s_init)) {
   throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_init must be finite.");
  }
  if (maf_effect_s_prior.size() != 2) {
   throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_prior must have length 2.");
  }
  if (!std::isfinite(maf_effect_s_prior[0]) || !std::isfinite(maf_effect_s_prior[1]) ||
      maf_effect_s_prior[0] >= maf_effect_s_prior[1]) {
   throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_prior must have finite lower < upper.");
  }
  if (maf_effect_s_init < maf_effect_s_prior[0] || maf_effect_s_init > maf_effect_s_prior[1]) {
   throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_init must lie within maf_effect_s_prior.");
  }
  if (!std::isfinite(maf_effect_s_proposal_sd) || maf_effect_s_proposal_sd <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_proposal_sd must be positive.");
  }
  if (maf_effect_s_log_h.isNull()) {
   throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_log_h is required when estimate_maf_effect_s = TRUE.");
  }
 }

 if ((int)ww.size() != nt || (int)b_init.size() != nt ||
     (int)yy.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr: inconsistent trait dimensions.");
 }

 if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr: priors must be nt x nt.");
 }

 if (pi.size() != 2) {
  throw std::runtime_error("stblr_cpg_omp_csr: pi must have length 2, c(pi0, pi1).");
 }

 if (!std::isfinite(pi_prior_a) || pi_prior_a <= 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr: pi_prior_a must be finite and positive.");
 }

 if (!std::isfinite(pi_prior_b) || pi_prior_b <= 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr: pi_prior_b must be finite and positive.");
 }

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m ||
      (int)ww[t].size() != m ||
      (int)b_init[t].size() != m) {
   throw std::runtime_error("stblr_cpg_omp_csr: inconsistent marker dimensions.");
  }
 }

 if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr: B must be nt x nt.");
 }

 if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr: E must be nt x nt.");
 }

 if (use_r_init) {
  if (static_cast<int>(r_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_state: r_init must have length nt when use_r_init = true.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(r_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_state: each r_init[t] must have length m.");
   }
  }
 }

 if (use_d_init) {
  if (static_cast<int>(d_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_state: d_init must have length nt when use_d_init = true.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(d_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_state: each d_init[t] must have length m.");
   }
  }
 }

 // --------------------------------------------------------------------------
 // Convert inputs to Armadillo
 // --------------------------------------------------------------------------

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);
 arma::rowvec prior_scale;
 arma::rowvec maf_effect_s_log_h_row;

 arma::vec yy_vec(nt, arma::fill::zeros);
 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  yy_vec(static_cast<arma::uword>(t)) = yy[t];

  for (int i = 0; i < m; ++i) {
   wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
   ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
   b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i))  = b_init[t][i];
  }

  if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr: priors must be nt x nt.");
  }

  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
   sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
  }
 }

 if (use_maf_effect_s_prior_scale) {
  prior_scale.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double scale_i = maf_effect_s_prior_scale_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(scale_i) || scale_i <= 0.0) {
    throw std::runtime_error(
      "stblr_cpg_omp_csr: maf_effect_s_prior_scale must contain positive finite values."
    );
   }
   prior_scale(static_cast<arma::uword>(i)) = scale_i;
  }
 }

 if (estimate_maf_effect_s) {
  Rcpp::NumericVector log_h_vec(maf_effect_s_log_h);
  if (static_cast<int>(log_h_vec.size()) != m) {
   throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_log_h must have length m.");
  }
  maf_effect_s_log_h_row.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double log_h_i = log_h_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(log_h_i)) {
    throw std::runtime_error("stblr_cpg_omp_csr: maf_effect_s_log_h must contain finite values.");
   }
   maf_effect_s_log_h_row(static_cast<arma::uword>(i)) = log_h_i;
  }
 }

 // --------------------------------------------------------------------------
 // Validate shared scaling and build shared flat LD object
 // --------------------------------------------------------------------------

 for (int t = 1; t < nt; ++t) {
  if (n[t] != n[0]) {
   throw std::runtime_error(
     "stblr_cpg_omp_csr: current shared-LD scaling assumes equal n across traits."
   );
  }
 }

 for (int t = 1; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   const double w0 = ww_mat(0, static_cast<arma::uword>(i));
   const double wt = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   const double tol = 1e-8 * std::max(1.0, std::abs(w0));

   if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
    throw std::runtime_error(
      "stblr_cpg_omp_csr: ww contains invalid value before LD pre-scaling."
    );
   }

   if (std::abs(w0 - wt) > tol) {
    throw std::runtime_error(
      "stblr_cpg_omp_csr: ww differs across traits; pre-scaled shared ST LD is invalid."
    );
   }
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);

 for (int i = 0; i < m; ++i) {
  const double wi = ww_mat(0, static_cast<arma::uword>(i));
  if (!std::isfinite(wi) || wi <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr: ww contains invalid value in trait 0.");
  }
  xx[static_cast<std::size_t>(i)] = wi;
 }

 arma::rowvec xx_row(static_cast<arma::uword>(m));
 for (int i = 0; i < m; ++i) {
  xx_row(static_cast<arma::uword>(i)) = xx[static_cast<std::size_t>(i)];
 }
 auto operator_context = make_operator(
  m,
  xx,
  xx_row,
  wy_mat,
  updateLDswap,
  ld_swap_r2,
  ld_swap_max_friends
 );
 const auto& op = operator_context.op;
 const LDLDFriends& ld_swap_friends = operator_context.ld_swap_friends;
 if (low_rank_residual_rebuild_every < 0) {
  throw std::runtime_error("low_rank_residual_rebuild_every must be non-negative.");
 }
 if (op.uses_retained_low_rank() && use_r_init) {
  throw std::runtime_error(
   "r_init is not supported for representation = \"retained_low_rank\"; "
   "a reduced-residual restart contract has not been implemented."
  );
 }
 if (op.uses_retained_low_rank() && rebuild_r_before_updateE) {
  throw std::runtime_error(
   "rebuild_r_before_updateE is incompatible with retained low rank; use "
   "low_rank_residual_rebuild_every instead."
  );
 }
 if (!op.uses_retained_low_rank() && low_rank_residual_rebuild_every != 0) {
  throw std::runtime_error(
   "low_rank_residual_rebuild_every is only supported by retained low rank."
  );
 }

 // --------------------------------------------------------------------------
 // Marker update order based on max single-trait marginal effect
 // --------------------------------------------------------------------------

 std::vector<double> x2(static_cast<std::size_t>(m), 0.0);
 std::vector<int> order(static_cast<std::size_t>(m));

 for (int i = 0; i < m; ++i) {
  double best = 0.0;

  for (int t = 0; t < nt; ++t) {
   const double wi = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   if (wi > 0.0) {
    const double bhat = wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) / wi;
    best = std::max(best, bhat * bhat);
   }
  }

  x2[static_cast<std::size_t>(i)] = best;
  order[static_cast<std::size_t>(i)] = i;
 }

 std::sort(order.begin(), order.end(),
           [&](int a, int b) { return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)]; });

 // Ordinary CSR has one compile-time route: the canonical typed core. The
 // other instantiation is the separate protected block-eigen backend; this is
 // not an old/new runtime selector or an ordinary-CSR fallback.
 if constexpr (std::is_same<
                 typename std::decay<decltype(op)>::type,
                 CsrOperator
               >::value) {
  return stblr_csr_bayesc_run_canonical(
   m, nt, ld_prefix, op, ld_swap_friends, wy_mat, b_mat, d_init,
   use_d_init, r_init, use_r_init, rebuild_r_before_updateE, yy_vec,
   B, E, ssb_prior_mat, sse_prior_mat, pi, nub, nue, updateB, updateE,
   updatePi, adjE, n, nit, nburn, nthin, pi_prior_a, pi_prior_b,
   ncores, seed, nchains, keep_chains, chain_seeds, updateLDswap,
   ld_swap_prob, ld_swap_moves, use_maf_effect_s_prior_scale, prior_scale,
   estimate_maf_effect_s, maf_effect_s_init, maf_effect_s_prior,
   maf_effect_s_proposal_sd, maf_effect_s_log_h_row, convergence_markers,
   convergence_b, convergence_d, order, execution_contract
  );
 } else {

 // --------------------------------------------------------------------------
 // Output storage
 // --------------------------------------------------------------------------

 const int ntasks = stblr_num_chain_tasks(nt, nchains);

 arma::mat bm_task(ntasks, m, arma::fill::zeros);
 arma::mat dm_task(ntasks, m, arma::fill::zeros);
 arma::mat b_task(ntasks, m, arma::fill::zeros);
 arma::mat r_task(ntasks, m, arma::fill::zeros);
 arma::mat d_task_double(ntasks, m, arma::fill::zeros);

 arma::mat vbs_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vgs_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat ves_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat pis_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vles_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vlds_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat maf_effect_s_task(ntasks, nit + nburn, arma::fill::zeros);

 arma::vec final_vb_task(ntasks, arma::fill::zeros);
 arma::vec final_vg_task(ntasks, arma::fill::zeros);
 arma::vec final_ve_task(ntasks, arma::fill::zeros);
 arma::vec final_pi_task(ntasks, arma::fill::zeros);
 arma::vec final_vle_task(ntasks, arma::fill::zeros);
 arma::vec final_vld_task(ntasks, arma::fill::zeros);
 arma::vec nsamples_task(ntasks, arma::fill::zeros);
 arma::vec ld_swap_attempted_task(ntasks, arma::fill::zeros);
 arma::vec ld_swap_accepted_task(ntasks, arma::fill::zeros);
 arma::vec maf_effect_s_attempted_task(ntasks, arma::fill::zeros);
 arma::vec maf_effect_s_accepted_task(ntasks, arma::fill::zeros);
 arma::ivec low_rank_residual_rebuild_count_task(ntasks, arma::fill::zeros);
 arma::vec low_rank_residual_max_abs_drift_task(ntasks, arma::fill::zeros);
 std::vector<arma::mat> convergence_b_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::imat> convergence_d_task(static_cast<std::size_t>(ntasks));
 const arma::uword convergence_marker_count=static_cast<arma::uword>(convergence_markers.size());
 for (int task=0;task<ntasks;++task) {
  if (convergence_b && convergence_marker_count>0) convergence_b_task[static_cast<std::size_t>(task)].zeros(nit,convergence_marker_count);
  if (convergence_d && convergence_marker_count>0) convergence_d_task[static_cast<std::size_t>(task)].zeros(nit,convergence_marker_count);
 }

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);
 arma::mat d_mat_double(nt, m, arma::fill::zeros);
 arma::mat bm_sd_mat(nt, m, arma::fill::zeros);
 arma::mat dm_sd_mat(nt, m, arma::fill::zeros);
 arma::mat bm_min_mat(nt, m, arma::fill::zeros);
 arma::mat dm_min_mat(nt, m, arma::fill::zeros);
 arma::mat bm_max_mat(nt, m, arma::fill::zeros);
 arma::mat dm_max_mat(nt, m, arma::fill::zeros);

 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vles_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vlds_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat maf_effect_s_mat(nt, nit + nburn, arma::fill::zeros);

 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi(nt, arma::fill::zeros);
 arma::vec final_vle(nt, arma::fill::zeros);
 arma::vec final_vld(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);
 arma::vec ld_swap_attempted_vec(nt, arma::fill::zeros);
 arma::vec ld_swap_accepted_vec(nt, arma::fill::zeros);
 arma::vec maf_effect_s_attempted_vec(nt, arma::fill::zeros);
 arma::vec maf_effect_s_accepted_vec(nt, arma::fill::zeros);

 // --------------------------------------------------------------------------
 // Parallel over traits
 // --------------------------------------------------------------------------

 std::vector<int> failed(static_cast<std::size_t>(ntasks), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(ntasks));
 std::vector<int> thread_used(static_cast<std::size_t>(ntasks), 0);
 std::vector<int> team_size_used(static_cast<std::size_t>(ntasks), 1);
 std::vector<double> task_seconds(static_cast<std::size_t>(ntasks), 0.0);

 int nthreads = 1;

#ifdef _OPENMP
 omp_set_dynamic(0);

 nthreads = stblr_num_threads_for_tasks(ncores, ntasks);

 // Explicitly set the OpenMP thread count for this runtime.
 omp_set_num_threads(nthreads);

#endif

 // #ifdef _OPENMP
 //  omp_set_dynamic(0);
 //  nthreads = std::max(1, std::min(ncores, nt));
 //  nthreads = std::min(nthreads, omp_get_max_threads());
 //
 //  Rcpp::Rcout
 //  << "STBLR OpenMP threads = "
 //  << nthreads
 //  << ", max threads = "
 //  << omp_get_max_threads()
 //  << ", num procs = "
 //  << omp_get_num_procs()
 //  << "\n";
 // #endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int task = 0; task < ntasks; ++task) {
  const int t = stblr_task_trait(task, nchains);
  const int chain = stblr_task_chain(task, nchains);

#ifdef _OPENMP
  const double wall_start = omp_get_wtime();
  thread_used[static_cast<std::size_t>(task)] = omp_get_thread_num();
  team_size_used[static_cast<std::size_t>(task)] = omp_get_num_threads();
#else
  const double wall_start = 0.0;
  thread_used[static_cast<std::size_t>(task)] = 0;
#endif

  try {
   unsigned int task_seed = 0u;
   if (!chain_seeds.empty()) {
    task_seed = stblr_seed_with_chain_base(
     chain_seeds[static_cast<std::size_t>(chain)],
     t
    );
   } else if (nchains == 1) {
    task_seed = stblr_trait_seed(seed, t);
   } else {
    task_seed = stblr_chain_seed(seed, t, chain);
   }
   task_seed = blr_phase3_task_seed(execution_contract, task, task_seed);
   std::mt19937 gen_t(task_seed);
   CsrBayesCPolicyHandle policy = policy_factory
    ? policy_factory->make(task, t, chain, m)
    : make_bayesc_noop_policy();

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   const arma::rowvec& ww_t = op.diag();

   arma::rowvec b_t(m, arma::fill::zeros);
   for (int i = 0; i < m; ++i) {
    b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   }

   arma::rowvec r_t(op.residual_size(), arma::fill::zeros);
   arma::Row<int> d_t(m, arma::fill::zeros);

   if (use_d_init) {
    for (int i = 0; i < m; ++i) {
     d_t(static_cast<arma::uword>(i)) = d_init[t][i] > 0 ? 1 : 0;
    }
   } else {
    for (int i = 0; i < m; ++i) {
     d_t(static_cast<arma::uword>(i)) =
      (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
    }
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) {
     r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
    }

    if (!r_t.is_finite()) {
     throw std::runtime_error("stblr_cpg_omp_csr_state: r_init contains NaN/Inf.");
    }
   } else {
    op.rebuild(t, wy_t, b_t, r_t);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = 0.0;
   double vle_t = computeLE_ST_csr(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   std::vector<double> pi_t = pi;

   if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
       pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
    throw std::runtime_error(
      "invalid initial pi: pi0=" + std::to_string(pi_t[0]) +
       ", pi1=" + std::to_string(pi_t[1])
    );
   }

   {
    const double psum = pi_t[0] + pi_t[1];

    if (!std::isfinite(psum) || psum <= 0.0) {
     throw std::runtime_error("invalid initial pi sum.");
    }

    pi_t[0] /= psum;
    pi_t[1] /= psum;
   }

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);
   arma::rowvec maf_effect_s_t(nit + nburn, arma::fill::zeros);
   arma::rowvec dynamic_prior_scale;

   double nsamples_t = 0.0;
   double ld_swap_attempted_t = 0.0;
   double ld_swap_accepted_t = 0.0;
   double maf_effect_s_current = maf_effect_s_init;
   double maf_effect_s_attempted_t = 0.0;
   double maf_effect_s_accepted_t = 0.0;
   int low_rank_residual_rebuild_count_t = 0;
   double low_rank_residual_max_abs_drift_t = 0.0;

   for (int it = 0; it < nit + nburn; ++it) {
    if (estimate_maf_effect_s) {
     fill_maf_effect_s_prior_scale(
      m,
      maf_effect_s_current,
      maf_effect_s_log_h_row,
      dynamic_prior_scale
     );
    }

    // -------------------------------------------------------
    // Marker updates
    // -------------------------------------------------------
    if (estimate_maf_effect_s || use_maf_effect_s_prior_scale ||
        policy.provides_prior_scale()) {
     const arma::rowvec& active_prior_scale =
      estimate_maf_effect_s ? dynamic_prior_scale :
       (use_maf_effect_s_prior_scale ? prior_scale : policy.prior_scale());
     for (int isort = 0; isort < m; ++isort) {
      const int i = order[static_cast<std::size_t>(isort)];

      sampleBetaC_ST_csr(
       i,
       pi_t,
       vb_t,
       active_prior_scale,
       vei_t,
       ww_t,
       r_t,
       b_t,
       d_t,
       op,
       gen_t
      );
     }
    } else {
     for (int isort = 0; isort < m; ++isort) {
      const int i = order[static_cast<std::size_t>(isort)];

      sampleBetaC_ST_csr_unscaled(
       i,
       pi_t,
       vb_t,
       vei_t,
       ww_t,
       r_t,
       b_t,
       d_t,
       op,
       gen_t
      );
     }
    }

    if (updateLDswap && ld_swap_moves > 0 && ld_swap_prob > 0.0) {
     std::uniform_real_distribution<double> runif(0.0, 1.0);

     if (runif(gen_t) < ld_swap_prob) {
      for (int move = 0; move < ld_swap_moves; ++move) {
       ld_swap_attempted_t += 1.0;
       const bool accepted = (estimate_maf_effect_s || use_maf_effect_s_prior_scale) ?
        attempt_ld_swap_st_csr(
            m,
            vei_t,
            vb_t,
            yy_vec(static_cast<arma::uword>(t)),
            ww_t,
            wy_t,
            estimate_maf_effect_s ? dynamic_prior_scale : prior_scale,
            r_t,
            b_t,
            d_t,
            op,
            ld_swap_friends,
            gen_t
           ) :
        attempt_ld_swap_st_csr_unscaled(
            m,
            vei_t,
            yy_vec(static_cast<arma::uword>(t)),
            ww_t,
            wy_t,
            r_t,
            b_t,
            d_t,
            op,
            ld_swap_friends,
            gen_t
           );
       if (accepted) {
        ld_swap_accepted_t += 1.0;
       }
      }
     }
    }

    // -------------------------------------------------------
    // Variance updates
    // -------------------------------------------------------
    if (updateB) {
     if (estimate_maf_effect_s || use_maf_effect_s_prior_scale ||
         policy.provides_prior_scale()) {
      sampleB_ST_csr(
       m,
       nub,
       vb_t,
       b_t,
       d_t,
       estimate_maf_effect_s ? dynamic_prior_scale :
        (use_maf_effect_s_prior_scale ? prior_scale : policy.prior_scale()),
       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
       gen_t
      );
     } else {
      sampleB_ST_csr_unscaled(
       m,
       nub,
       vb_t,
       b_t,
       d_t,
       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
       gen_t
      );
     }

     if (!std::isfinite(vb_t) || vb_t <= 0.0) {
      throw std::runtime_error(
        "vb became invalid after sampleB. iter=" +
         std::to_string(it) +
       ", vb=" + std::to_string(vb_t)
      );
     }
    }

    if (estimate_maf_effect_s) {
     maf_effect_s_attempted_t += 1.0;
     const bool accepted_s = update_maf_effect_s_bayesc(
      maf_effect_s_current,
      b_t,
      d_t,
      vb_t,
      maf_effect_s_log_h_row,
      maf_effect_s_prior[0],
      maf_effect_s_prior[1],
      maf_effect_s_proposal_sd,
      gen_t
     );
     if (accepted_s) maf_effect_s_accepted_t += 1.0;
    }

    policy.after_vb_update(b_t, d_t, vb_t, gen_t, it);

    if (op.uses_retained_low_rank() &&
        low_rank_residual_rebuild_every > 0 &&
        ((it + 1) % low_rank_residual_rebuild_every == 0)) {
     const double drift = op.rebuild_and_measure_drift(t, wy_t, b_t, r_t);
     low_rank_residual_max_abs_drift_t = std::max(
      low_rank_residual_max_abs_drift_t, drift
     );
     ++low_rank_residual_rebuild_count_t;
    }

    if (updateE) {
     if (rebuild_r_before_updateE) {
      op.rebuild(t, wy_t, b_t, r_t);
     }

     sampleE_ST_operator(
      op,
      t,
      nue,
      ve_t,
      b_t,
      wy_t,
      r_t,
      sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      yy_vec(static_cast<arma::uword>(t)),
      n[t],
       gen_t
     );
    }

    if (updatePi) {
     samplePi_ST(d_t, pi_t, pi_prior_a, pi_prior_b, gen_t);

     if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
         pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
      throw std::runtime_error(
        "pi became invalid after samplePi. iter=" +
         std::to_string(it) +
         ", pi0=" + std::to_string(pi_t[0]) +
         ", pi1=" + std::to_string(pi_t[1])
      );
     }
    }

    vg_t = computeG_ST_operator(
     op,
     t,
     b_t,
     wy_t,
     r_t,
     n[t]
    );

    vle_t = computeLE_ST_csr(
     m,
     b_t,
     ww_t,
     n[t]
    );

    vld_t = vg_t - vle_t;

    if (!std::isfinite(vg_t)) {
     throw std::runtime_error(
       "vg became NaN/Inf after computeG. iter=" +
        std::to_string(it)
     );
    }

    if (!std::isfinite(vle_t)) {
     throw std::runtime_error(
       "vle became NaN/Inf after computeLE. iter=" +
        std::to_string(it)
     );
    }

    if (!std::isfinite(vld_t)) {
     throw std::runtime_error(
       "vld became NaN/Inf after computeLE. iter=" +
        std::to_string(it)
     );
    }

    vei_t = ve_t + adjE * vg_t;

    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error(
       "adjusted residual variance vei became invalid. iter=" +
        std::to_string(it) +
        ", vei=" + std::to_string(vei_t)
     );
    }

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    pis_t(static_cast<arma::uword>(it)) = pi_t[1];
    vles_t(static_cast<arma::uword>(it)) = vle_t;
    vlds_t(static_cast<arma::uword>(it)) = vld_t;
    maf_effect_s_t(static_cast<arma::uword>(it)) = maf_effect_s_current;
    policy.capture(it);
    if (it>=nburn) {
     const arma::uword draw=static_cast<arma::uword>(it-nburn);
     for (arma::uword s=0;s<convergence_marker_count;++s) {
      const arma::uword marker=static_cast<arma::uword>(convergence_markers[static_cast<std::size_t>(s)]);
      if (convergence_b) convergence_b_task[static_cast<std::size_t>(task)](draw,s)=b_t(marker);
      if (convergence_d) convergence_d_task[static_cast<std::size_t>(task)](draw,s)=d_t(marker);
     }
    }

    // -------------------------------------------------------
    // Store posterior summaries
    // -------------------------------------------------------
    if (blr_phase3_iteration_is_retained(
        execution_contract, it, nburn, nthin)) {
     policy.retain(it);
     nsamples_t += 1.0;

     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += static_cast<double>(d_t(iu));
     }
    }
   }

   policy.finish();

   if (op.uses_retained_low_rank()) {
    const double drift = op.rebuild_and_measure_drift(t, wy_t, b_t, r_t);
    low_rank_residual_max_abs_drift_t = std::max(
     low_rank_residual_max_abs_drift_t, drift
    );
    ++low_rank_residual_rebuild_count_t;
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;

   bm_t /= nsamples_t;
   dm_t /= nsamples_t;

   if (!bm_t.is_finite()) {
    throw std::runtime_error("posterior mean bm contains NaN/Inf.");
   }

   if (!dm_t.is_finite()) {
    throw std::runtime_error("posterior mean dm contains NaN/Inf.");
   }

   const arma::uword task_u = static_cast<arma::uword>(task);
   bm_task.row(task_u) = bm_t;
   dm_task.row(task_u) = dm_t;
   b_task.row(task_u)  = b_t;
   arma::rowvec marker_residual;
   op.materialize_residual(t, r_t, marker_residual);
   r_task.row(task_u)  = marker_residual;
   for (int i = 0; i < m; ++i) {
    d_task_double(task_u, static_cast<arma::uword>(i)) =
     static_cast<double>(d_t(static_cast<arma::uword>(i)));
   }

   vbs_task.row(task_u) = vbs_t;
   vgs_task.row(task_u) = vgs_t;
   ves_task.row(task_u) = ves_t;
   pis_task.row(task_u) = pis_t;
   vles_task.row(task_u) = vles_t;
   vlds_task.row(task_u) = vlds_t;
   maf_effect_s_task.row(task_u) = maf_effect_s_t;

   final_vb_task(task_u) = vb_t;
   final_ve_task(task_u) = ve_t;
   final_vg_task(task_u) = vg_t;
   final_vle_task(task_u) = vle_t;
   final_vld_task(task_u) = vld_t;
   final_pi_task(task_u) = pi_t[1];
   nsamples_task(task_u) = nsamples_t;
   ld_swap_attempted_task(task_u) = ld_swap_attempted_t;
   ld_swap_accepted_task(task_u) = ld_swap_accepted_t;
   maf_effect_s_attempted_task(task_u) = maf_effect_s_attempted_t;
   maf_effect_s_accepted_task(task_u) = maf_effect_s_accepted_t;
   low_rank_residual_rebuild_count_task(task_u) =
    low_rank_residual_rebuild_count_t;
   low_rank_residual_max_abs_drift_task(task_u) =
    low_rank_residual_max_abs_drift_t;

#ifdef _OPENMP
   task_seconds[static_cast<std::size_t>(task)] = omp_get_wtime() - wall_start;
#endif

  } catch (const std::exception& e) {
   failed[static_cast<std::size_t>(task)] = 1;
   errors[static_cast<std::size_t>(task)] = e.what();
#ifdef _OPENMP
   task_seconds[static_cast<std::size_t>(task)] = omp_get_wtime() - wall_start;
#endif
  } catch (...) {
   failed[static_cast<std::size_t>(task)] = 1;
   errors[static_cast<std::size_t>(task)] = "unknown error";
#ifdef _OPENMP
   task_seconds[static_cast<std::size_t>(task)] = omp_get_wtime() - wall_start;
#endif
  }
 }

#ifdef _OPENMP
 for (int task = 0; task < ntasks; ++task) {
  const int t = stblr_task_trait(task, nchains);
  const int chain = stblr_task_chain(task, nchains);
 }
#endif

 for (int task = 0; task < ntasks; ++task) {
  if (failed[static_cast<std::size_t>(task)]) {
   const int t = stblr_task_trait(task, nchains);
   const int chain = stblr_task_chain(task, nchains);
   throw std::runtime_error(
     "stblr_cpg_omp_csr failed for trait " +
      std::to_string(t) +
      ", chain " +
      std::to_string(chain) +
      ": " +
      errors[static_cast<std::size_t>(task)]
   );
  }
 }

 const double inv_chains = 1.0 / static_cast<double>(nchains);

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);

  bm_min_mat.row(tu).fill(std::numeric_limits<double>::infinity());
  dm_min_mat.row(tu).fill(std::numeric_limits<double>::infinity());
  bm_max_mat.row(tu).fill(-std::numeric_limits<double>::infinity());
  dm_max_mat.row(tu).fill(-std::numeric_limits<double>::infinity());

  for (int chain = 0; chain < nchains; ++chain) {
   const int task = t * nchains + chain;
   const arma::uword task_u = static_cast<arma::uword>(task);

   bm_mat.row(tu) += bm_task.row(task_u);
   dm_mat.row(tu) += dm_task.row(task_u);
   b_mat.row(tu) += b_task.row(task_u);
   r_mat.row(tu) += r_task.row(task_u);
   d_mat_double.row(tu) += d_task_double.row(task_u);

   vbs_mat.row(tu) += vbs_task.row(task_u);
   vgs_mat.row(tu) += vgs_task.row(task_u);
   ves_mat.row(tu) += ves_task.row(task_u);
   pis_mat.row(tu) += pis_task.row(task_u);
   vles_mat.row(tu) += vles_task.row(task_u);
   vlds_mat.row(tu) += vlds_task.row(task_u);
   maf_effect_s_mat.row(tu) += maf_effect_s_task.row(task_u);

   final_vb(tu) += final_vb_task(task_u);
   final_vg(tu) += final_vg_task(task_u);
   final_ve(tu) += final_ve_task(task_u);
   final_vle(tu) += final_vle_task(task_u);
   final_vld(tu) += final_vld_task(task_u);
   final_pi(tu) += final_pi_task(task_u);
   nsamples_vec(tu) += nsamples_task(task_u);
   ld_swap_attempted_vec(tu) += ld_swap_attempted_task(task_u);
   ld_swap_accepted_vec(tu) += ld_swap_accepted_task(task_u);
   maf_effect_s_attempted_vec(tu) += maf_effect_s_attempted_task(task_u);
   maf_effect_s_accepted_vec(tu) += maf_effect_s_accepted_task(task_u);

   for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    bm_min_mat(tu, iu) = std::min(bm_min_mat(tu, iu), bm_task(task_u, iu));
    dm_min_mat(tu, iu) = std::min(dm_min_mat(tu, iu), dm_task(task_u, iu));
    bm_max_mat(tu, iu) = std::max(bm_max_mat(tu, iu), bm_task(task_u, iu));
    dm_max_mat(tu, iu) = std::max(dm_max_mat(tu, iu), dm_task(task_u, iu));
   }
  }

  bm_mat.row(tu) *= inv_chains;
  dm_mat.row(tu) *= inv_chains;
  b_mat.row(tu) *= inv_chains;
  r_mat.row(tu) *= inv_chains;
  d_mat_double.row(tu) *= inv_chains;
  vbs_mat.row(tu) *= inv_chains;
  vgs_mat.row(tu) *= inv_chains;
  ves_mat.row(tu) *= inv_chains;
  pis_mat.row(tu) *= inv_chains;
  vles_mat.row(tu) *= inv_chains;
  vlds_mat.row(tu) *= inv_chains;
  maf_effect_s_mat.row(tu) *= inv_chains;
  final_vb(tu) *= inv_chains;
  final_vg(tu) *= inv_chains;
  final_ve(tu) *= inv_chains;
  final_vle(tu) *= inv_chains;
  final_vld(tu) *= inv_chains;
  final_pi(tu) *= inv_chains;
  nsamples_vec(tu) *= inv_chains;

  if (nchains > 1) {
   for (int chain = 0; chain < nchains; ++chain) {
    const int task = t * nchains + chain;
    const arma::uword task_u = static_cast<arma::uword>(task);
    arma::rowvec bm_diff = bm_task.row(task_u) - bm_mat.row(tu);
    arma::rowvec dm_diff = dm_task.row(task_u) - dm_mat.row(tu);
    bm_sd_mat.row(tu) += bm_diff % bm_diff;
    dm_sd_mat.row(tu) += dm_diff % dm_diff;
   }
   bm_sd_mat.row(tu) = arma::sqrt(bm_sd_mat.row(tu) / static_cast<double>(nchains - 1));
   dm_sd_mat.row(tu) = arma::sqrt(dm_sd_mat.row(tu) / static_cast<double>(nchains - 1));
  }
 }

 // --------------------------------------------------------------------------
 // Build named raw schema v1
 // --------------------------------------------------------------------------

 const int n_trace = nit + nburn;
 const bool return_chain_summaries = (nchains > 1) || keep_chains;

 auto marker_matrix = [&](const arma::mat& x) {
  Rcpp::NumericMatrix out(m, nt);
  for (int t = 0; t < nt; ++t) {
   for (int i = 0; i < m; ++i) {
    out(i, t) = x(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   }
  }
  return out;
 };

 auto trace_matrix = [&](const arma::mat& x) {
  Rcpp::NumericMatrix out(n_trace, nt);
  for (int t = 0; t < nt; ++t) {
   for (int it = 0; it < n_trace; ++it) {
    out(it, t) = x(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
   }
  }
  return out;
 };

 auto diagonal_matrix = [&](const arma::vec& x) {
  Rcpp::NumericMatrix out(nt, nt);
  for (int t = 0; t < nt; ++t) {
   out(t, t) = x(static_cast<arma::uword>(t));
  }
  return out;
 };

 Rcpp::NumericMatrix pi_final(nt, 2);
 Rcpp::NumericMatrix pi_mean(nt, 2);
 Rcpp::NumericVector selection_mean(nt);
 Rcpp::NumericVector maf_effect_sd(nt);
 Rcpp::NumericVector selection_min(nt);
 Rcpp::NumericVector selection_max(nt);
 Rcpp::NumericVector selection_acceptance(nt);

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  pi_final(t, 0) = 1.0 - final_pi(tu);
  pi_final(t, 1) = final_pi(tu);

  double mean_pi = 0.0;
  int npi = 0;
  for (int it = nburn; it < n_trace; ++it) {
   mean_pi += pis_mat(tu, static_cast<arma::uword>(it));
   ++npi;
  }
  if (npi > 0) mean_pi /= static_cast<double>(npi);
  else mean_pi = final_pi(tu);
  pi_mean(t, 0) = 1.0 - mean_pi;
  pi_mean(t, 1) = mean_pi;

  selection_acceptance[t] =
   maf_effect_s_attempted_vec(tu) > 0.0
   ? maf_effect_s_accepted_vec(tu) / maf_effect_s_attempted_vec(tu)
   : 0.0;

  if (estimate_maf_effect_s) {
   double mean_s = 0.0;
   double min_s = std::numeric_limits<double>::infinity();
   double max_s = -std::numeric_limits<double>::infinity();
   int ns = 0;
   for (int it = nburn; it < n_trace; ++it) {
    const double val = maf_effect_s_mat(tu, static_cast<arma::uword>(it));
    mean_s += val;
    min_s = std::min(min_s, val);
    max_s = std::max(max_s, val);
    ++ns;
   }
   if (ns > 0) mean_s /= static_cast<double>(ns);
   double ss = 0.0;
   if (ns > 1) {
    for (int it = nburn; it < n_trace; ++it) {
     const double diff = maf_effect_s_mat(tu, static_cast<arma::uword>(it)) - mean_s;
     ss += diff * diff;
    }
    maf_effect_sd[t] = std::sqrt(ss / static_cast<double>(ns - 1));
   } else {
    maf_effect_sd[t] = NA_REAL;
   }
   selection_mean[t] = mean_s;
   selection_min[t] = ns > 0 ? min_s : NA_REAL;
   selection_max[t] = ns > 0 ? max_s : NA_REAL;
  } else {
   selection_mean[t] = NA_REAL;
   maf_effect_sd[t] = NA_REAL;
   selection_min[t] = NA_REAL;
   selection_max[t] = NA_REAL;
  }
 }

 Rcpp::NumericMatrix ld_swap(nt, 3);
 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  const double attempted = ld_swap_attempted_vec(tu);
  const double accepted = ld_swap_accepted_vec(tu);
  ld_swap(t, 0) = attempted;
  ld_swap(t, 1) = accepted;
  ld_swap(t, 2) = attempted > 0.0 ? accepted / attempted : 0.0;
 }

 Rcpp::NumericVector nsamples(nt);
 Rcpp::IntegerVector n_used(nt);
 Rcpp::NumericVector seconds_mean(nt);
 Rcpp::NumericVector seconds_max(nt);
 Rcpp::NumericVector low_rank_residual_rebuild_count(nt);
 Rcpp::NumericVector low_rank_residual_max_abs_drift(nt);
 for (int t = 0; t < nt; ++t) {
  nsamples[t] = nsamples_vec(static_cast<arma::uword>(t));
  n_used[t] = n[static_cast<std::size_t>(t)];
  double sec_sum = 0.0;
  double sec_max = 0.0;
  for (int chain = 0; chain < nchains; ++chain) {
   const int task = t * nchains + chain;
   const double sec = task_seconds[static_cast<std::size_t>(task)];
   sec_sum += sec;
   sec_max = std::max(sec_max, sec);
  }
  seconds_mean[t] = sec_sum / static_cast<double>(nchains);
   seconds_max[t] = sec_max;
  for (int chain = 0; chain < nchains; ++chain) {
   const arma::uword task_u = static_cast<arma::uword>(t * nchains + chain);
   low_rank_residual_rebuild_count[t] = std::max(
    low_rank_residual_rebuild_count[t],
    static_cast<double>(low_rank_residual_rebuild_count_task(task_u))
   );
   low_rank_residual_max_abs_drift[t] = std::max(
    low_rank_residual_max_abs_drift[t],
    low_rank_residual_max_abs_drift_task(task_u)
   );
  }
 }

 Rcpp::List marker = Rcpp::List::create(
  Rcpp::Named("bm") = marker_matrix(bm_mat),
  Rcpp::Named("dm") = marker_matrix(dm_mat),
  Rcpp::Named("wy") = marker_matrix(wy_mat),
  Rcpp::Named("r") = marker_matrix(r_mat),
  Rcpp::Named("b") = marker_matrix(b_mat),
  Rcpp::Named("state") = marker_matrix(d_mat_double)
 );
 if (return_chain_summaries) {
  marker["bm_sd"] = marker_matrix(bm_sd_mat);
  marker["bm_min"] = marker_matrix(bm_min_mat);
  marker["bm_max"] = marker_matrix(bm_max_mat);
  marker["dm_sd"] = marker_matrix(dm_sd_mat);
  marker["dm_min"] = marker_matrix(dm_min_mat);
  marker["dm_max"] = marker_matrix(dm_max_mat);
 }

 Rcpp::List trace = Rcpp::List::create(
  Rcpp::Named("vbs") = trace_matrix(vbs_mat),
  Rcpp::Named("vgs") = trace_matrix(vgs_mat),
  Rcpp::Named("ves") = trace_matrix(ves_mat),
  Rcpp::Named("vle") = trace_matrix(vles_mat),
  Rcpp::Named("vld") = trace_matrix(vlds_mat),
  Rcpp::Named("pis") = trace_matrix(pis_mat)
 );

 Rcpp::List variance = Rcpp::List::create(
  Rcpp::Named("covb") = diagonal_matrix(final_vb),
  Rcpp::Named("covg") = diagonal_matrix(final_vg),
  Rcpp::Named("cove") = diagonal_matrix(final_ve),
  Rcpp::Named("vb") = diagonal_matrix(final_vb),
  Rcpp::Named("vg") = diagonal_matrix(final_vg),
  Rcpp::Named("ve") = diagonal_matrix(final_ve)
 );

 Rcpp::List diagnostics = Rcpp::List::create(
  Rcpp::Named("nsamples") = nsamples,
  Rcpp::Named("n_used") = n_used,
  Rcpp::Named("log_cpo") = Rcpp::NumericVector(nt),
  Rcpp::Named("mean_log_cpo") = Rcpp::NumericVector(nt),
  Rcpp::Named("seconds_mean") = seconds_mean,
  Rcpp::Named("seconds_max") = seconds_max,
  Rcpp::Named("ld_swap") = updateLDswap ? Rcpp::wrap(ld_swap) : R_NilValue
 );
 if (execution_contract.active()) {
  diagnostics.push_back(blr_phase3_worker_diagnostics(
   execution_contract, ncores, nthreads, thread_used, team_size_used),
   "workers");
 }
 if (operator_context.diagnostics.size() > 0) {
  diagnostics["block_eigen"] = operator_context.diagnostics;
 }
 if (op.uses_retained_low_rank()) {
  diagnostics["low_rank_residual"] = Rcpp::List::create(
   Rcpp::Named("low_rank_residual_rebuild_every") =
    Rcpp::IntegerVector(nt, low_rank_residual_rebuild_every),
   Rcpp::Named("low_rank_residual_rebuild_count") =
    low_rank_residual_rebuild_count,
   Rcpp::Named("low_rank_residual_max_abs_drift") =
    low_rank_residual_max_abs_drift
  );
 }

 Rcpp::List chains = R_NilValue;
 if (keep_chains) {
  chains = Rcpp::List(nt);
  Rcpp::CharacterVector trait_names(nt);
  for (int t = 0; t < nt; ++t) {
   trait_names[t] = "trait" + std::to_string(t + 1);
   Rcpp::List trait_chains(nchains);
   Rcpp::CharacterVector chain_names(nchains);
   for (int chain = 0; chain < nchains; ++chain) {
    chain_names[chain] = "chain" + std::to_string(chain + 1);
    const int task = t * nchains + chain;
    const arma::uword task_u = static_cast<arma::uword>(task);
    Rcpp::NumericVector chain_bm(m);
    Rcpp::NumericVector chain_dm(m);
    Rcpp::NumericVector chain_state(m);
    for (int i = 0; i < m; ++i) {
     const arma::uword iu = static_cast<arma::uword>(i);
     chain_bm[i] = bm_task(task_u, iu);
     chain_dm[i] = dm_task(task_u, iu);
     chain_state[i] = d_task_double(task_u, iu);
    }
    Rcpp::NumericMatrix chain_ld(1, 3);
    const double attempted = ld_swap_attempted_task(task_u);
    const double accepted = ld_swap_accepted_task(task_u);
    chain_ld(0, 0) = attempted;
    chain_ld(0, 1) = accepted;
    chain_ld(0, 2) = attempted > 0.0 ? accepted / attempted : 0.0;
    Rcpp::List chain_low_rank_residual = R_NilValue;
    if (op.uses_retained_low_rank()) {
     chain_low_rank_residual = Rcpp::List::create(
      Rcpp::Named("low_rank_residual_rebuild_every") =
       low_rank_residual_rebuild_every,
      Rcpp::Named("low_rank_residual_rebuild_count") =
       low_rank_residual_rebuild_count_task(task_u),
      Rcpp::Named("low_rank_residual_max_abs_drift") =
       low_rank_residual_max_abs_drift_task(task_u)
     );
    }
    Rcpp::List chain_selection = Rcpp::List::create(
     Rcpp::Named("trace") = R_NilValue,
     Rcpp::Named("acceptance") = R_NilValue
    );
    if (estimate_maf_effect_s) {
     Rcpp::NumericVector chain_s_trace(n_trace);
     for (int it = 0; it < n_trace; ++it) {
      chain_s_trace[it] = maf_effect_s_task(task_u, static_cast<arma::uword>(it));
     }
     const double s_attempted = maf_effect_s_attempted_task(task_u);
     const double s_accepted = maf_effect_s_accepted_task(task_u);
     chain_selection["trace"] = chain_s_trace;
     chain_selection["acceptance"] =
      s_attempted > 0.0 ? s_accepted / s_attempted : 0.0;
    }
    trait_chains[chain] = Rcpp::List::create(
     Rcpp::Named("marker") = Rcpp::List::create(
      Rcpp::Named("bm") = chain_bm,
      Rcpp::Named("dm") = chain_dm,
      Rcpp::Named("state") = chain_state
     ),
     Rcpp::Named("trace") = Rcpp::List::create(
      Rcpp::Named("vbs") = Rcpp::wrap(vbs_task.row(task_u)),
      Rcpp::Named("vgs") = Rcpp::wrap(vgs_task.row(task_u)),
      Rcpp::Named("ves") = Rcpp::wrap(ves_task.row(task_u)),
      Rcpp::Named("vle") = Rcpp::wrap(vles_task.row(task_u)),
      Rcpp::Named("vld") = Rcpp::wrap(vlds_task.row(task_u)),
      Rcpp::Named("pis") = Rcpp::wrap(pis_task.row(task_u))
     ),
     Rcpp::Named("pi") = Rcpp::List::create(
      Rcpp::Named("final") = Rcpp::NumericVector::create(
       1.0 - final_pi_task(task_u), final_pi_task(task_u)
      ),
      Rcpp::Named("mean") = R_NilValue
     ),
     Rcpp::Named("selection") = chain_selection,
     Rcpp::Named("convergence_trace") = Rcpp::List::create(
      Rcpp::Named("b") = convergence_b_task[static_cast<std::size_t>(task)],
      Rcpp::Named("d") = convergence_d_task[static_cast<std::size_t>(task)],
      Rcpp::Named("component") = R_NilValue,
      Rcpp::Named("marker_index") = Rcpp::wrap(convergence_markers)),
     Rcpp::Named("diagnostics") = Rcpp::List::create(
      Rcpp::Named("ld_swap") = updateLDswap ? Rcpp::wrap(chain_ld) : R_NilValue,
      Rcpp::Named("low_rank_residual") = chain_low_rank_residual
     )
    );
   }
   trait_chains.attr("names") = chain_names;
   chains[t] = trait_chains;
  }
  chains.attr("names") = trait_names;
 }

 Rcpp::List selection = Rcpp::List::create(
  Rcpp::Named("enabled") = estimate_maf_effect_s || use_maf_effect_s_prior_scale,
  Rcpp::Named("fixed") = use_maf_effect_s_prior_scale,
  Rcpp::Named("scale") = "standardized_genotype_effect",
  Rcpp::Named("trace") = estimate_maf_effect_s ? Rcpp::wrap(trace_matrix(maf_effect_s_mat)) : R_NilValue,
  Rcpp::Named("mean") = estimate_maf_effect_s ? Rcpp::wrap(selection_mean) : R_NilValue,
  Rcpp::Named("sd") = estimate_maf_effect_s ? Rcpp::wrap(maf_effect_sd) : R_NilValue,
  Rcpp::Named("min") = estimate_maf_effect_s ? Rcpp::wrap(selection_min) : R_NilValue,
  Rcpp::Named("max") = estimate_maf_effect_s ? Rcpp::wrap(selection_max) : R_NilValue,
  Rcpp::Named("acceptance") = estimate_maf_effect_s ? Rcpp::wrap(selection_acceptance) : R_NilValue
 );

 Rcpp::List raw = Rcpp::List::create(
  Rcpp::Named("schema") = Rcpp::List::create(
   Rcpp::Named("class") = "stblr_raw",
   Rcpp::Named("version") = 1
  ),
  Rcpp::Named("meta") = Rcpp::List::create(
   Rcpp::Named("model") = "bayesc",
   Rcpp::Named("backend") = "csr_bayesc",
   Rcpp::Named("data_level") = "summary",
   Rcpp::Named("prior_type") = "global",
   Rcpp::Named("m") = m,
   Rcpp::Named("nt") = nt,
   Rcpp::Named("n_trace") = n_trace,
   Rcpp::Named("nit") = nit,
   Rcpp::Named("nburn") = nburn,
   Rcpp::Named("nthin") = nthin,
   Rcpp::Named("nchains") = nchains,
   Rcpp::Named("keep_chains") = keep_chains,
   Rcpp::Named("n_components") = 2,
   Rcpp::Named("n_annotations") = 0,
   Rcpp::Named("n_groups") = 0
  ),
  Rcpp::Named("marker") = marker,
  Rcpp::Named("trace") = trace,
  Rcpp::Named("variance") = variance,
  Rcpp::Named("pi") = Rcpp::List::create(
   Rcpp::Named("final") = pi_final,
   Rcpp::Named("mean") = pi_mean,
   Rcpp::Named("names") = Rcpp::CharacterVector::create("pi0", "pi1")
  ),
  Rcpp::Named("diagnostics") = diagnostics,
  Rcpp::Named("chains") = chains,
  Rcpp::Named("prior") = Rcpp::List::create(),
  Rcpp::Named("group") = Rcpp::List::create(),
  Rcpp::Named("annotation") = Rcpp::List::create(),
  Rcpp::Named("component") = Rcpp::List::create(),
  Rcpp::Named("selection") = selection
 );
 raw.attr("class") = Rcpp::CharacterVector::create("stblr_raw_v1", "stblr_raw", "list");
 return raw;
 }
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> d_init,
  bool use_d_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  std::string ld_prefix,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<double> pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  double adjE,
  std::vector<int> n,
  int nit,
  int nburn,
  int nthin,
  double pi_prior_a,
  double pi_prior_b,
  int ncores,
  int seed,
  int nchains,
  bool keep_chains,
  std::vector<int> chain_seeds,
  bool updateLDswap = false,
  double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8,
  int ld_swap_max_friends = 50,
  int ld_swap_moves = 1,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale = R_NilValue,
  bool estimate_maf_effect_s = false,
  double maf_effect_s_init = 0.0,
  Rcpp::NumericVector maf_effect_s_prior = Rcpp::NumericVector::create(-3.0, 2.0),
  double maf_effect_s_proposal_sd = 0.35,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h = R_NilValue,
  Rcpp::IntegerVector convergence_markers = Rcpp::IntegerVector::create(),
  bool convergence_b = false,
  bool convergence_d = false,
  Rcpp::Nullable<Rcpp::List> execution_contract = R_NilValue
) {
 auto make_csr_operator = [&](int m,
                              const std::vector<double>& xx,
                              const arma::rowvec& xx_row,
                              arma::mat& wy_mat,
                              bool update_ld_swap,
                              double ld_swap_r2_value,
                              int ld_swap_max_friends_value) {
  (void)wy_mat;
  STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
  CsrOperator op(ld, xx_row);
  LDLDFriends ld_swap_friends;
  if (update_ld_swap) {
   ld_swap_friends = build_ld_swap_friends_st_csr(
    m,
    ld,
    xx,
    ld_swap_r2_value,
    ld_swap_max_friends_value
   );
  } else {
   ld_swap_friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  }
  return BayescOperatorContext<CsrOperator>(
   op,
   ld_swap_friends,
   Rcpp::List::create()
  );
 };

 const int expected_tasks = static_cast<int>(wy.size()) * nchains;
 const BlrPhase3ExecutionContract phase3 =
  parse_blr_phase3_execution_contract(execution_contract, expected_tasks, nit);
 return stblr_cpg_omp_csr_impl(
  wy, ww, yy, b_init, d_init, use_d_init, r_init, use_r_init,
  rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, pi,
  nub, nue, updateB, updateE, updatePi, adjE, n, nit, nburn, nthin,
  pi_prior_a, pi_prior_b, ncores, seed, nchains, keep_chains, chain_seeds,
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves,
  maf_effect_s_prior_scale, estimate_maf_effect_s, maf_effect_s_init,
  maf_effect_s_prior, maf_effect_s_proposal_sd, maf_effect_s_log_h,
  Rcpp::as<std::vector<int>>(convergence_markers), convergence_b,
  convergence_d,
  0, nullptr, phase3, make_csr_operator
 );
}

Rcpp::List stblr_cpg_omp_csr_block_eigen_with_policy(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> d_init,
  bool use_d_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  std::string ld_prefix,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<double> pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  double adjE,
  std::vector<int> n,
  int nit,
  int nburn,
  int nthin,
  double pi_prior_a,
  double pi_prior_b,
  int ncores,
  int seed,
  int nchains,
  bool keep_chains,
  std::vector<int> chain_seeds,
  bool updateLDswap,
  double ld_swap_prob,
  double ld_swap_r2,
  int ld_swap_max_friends,
  int ld_swap_moves,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale,
  bool estimate_maf_effect_s,
  double maf_effect_s_init,
  Rcpp::NumericVector maf_effect_s_prior,
  double maf_effect_s_proposal_sd,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h,
  Rcpp::IntegerVector convergence_markers,
  bool convergence_b,
  bool convergence_d,
  Rcpp::CharacterVector bed_files,
  int n_bed,
  Rcpp::List cls,
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::NumericVector af,
  Rcpp::IntegerVector block_start,
  std::string eigen_filter,
  double eigen_tau,
  double eigen_eta,
  std::string representation,
  double eigen_prop,
  int low_rank_residual_rebuild_every,
  CsrBayesCPolicyFactory* policy_factory,
  const BlrPhase3ExecutionContract& execution_contract
) {
 if (updateLDswap) {
  throw std::runtime_error(
   "LD-swap is not yet supported with the experimental block-eigen operator."
  );
 }

 if (bed_files.size() <= 0) {
  throw std::runtime_error("bed_files must contain at least one BED file.");
 }
 if (n_bed <= 0) {
  throw std::runtime_error("n_bed must be positive.");
 }
 if (cls.size() != bed_files.size()) {
  throw std::runtime_error("cls must have one element per BED file.");
 }

 const std::vector<std::string> bed_files_cpp =
  stblr_copy_character_vector(bed_files, "bed_files");
 const std::vector<std::vector<int>> cls_cpp =
  stblr_copy_int_list(cls, "cls");
 const std::vector<int> rows0 = stblr_copy_rows0_or_empty(rows, n_bed);
 const std::vector<double> af_cpp = stblr_copy_numeric_vector(af, "af");
 const std::vector<int> block_start_cpp =
  stblr_copy_integer_vector(block_start, "block_start");
 const EigenFilterMode mode = parse_block_eigen_filter_mode(eigen_filter);
 if (representation != "low_rank" && representation != "dense_reconstructed")
  throw std::runtime_error("unknown block-eigen representation.");
 if (representation == "low_rank" && use_r_init)
  throw std::runtime_error(
   "r_init is not supported for representation = \"retained_low_rank\"; "
   "a reduced-residual restart contract has not been implemented."
  );

 BlockEigenExecutionInput block_input;
 block_input.bed_files = bed_files_cpp;
 block_input.n_bed = n_bed;
 block_input.cls = cls_cpp;
 block_input.rows0 = rows0;
 block_input.af = af_cpp;
 block_input.block_start = block_start_cpp;
 block_input.filter_mode = mode;
 block_input.eigen_tau = eigen_tau;
 block_input.eigen_eta = eigen_eta;
 block_input.low_rank = representation == "low_rank";
 block_input.eigen_prop = eigen_prop;
 block_input.ncores = ncores;

 auto make_block_eigen_operator = [&](int m,
                                      const std::vector<double>& xx,
                                      const arma::rowvec& xx_row,
                                      arma::mat& wy_mat,
                                      bool update_ld_swap,
                                      double ld_swap_r2_value,
                                      int ld_swap_max_friends_value) {
  (void)xx;
  (void)xx_row;
  (void)ld_swap_r2_value;
  (void)ld_swap_max_friends_value;
  if (update_ld_swap) {
   throw std::runtime_error(
    "LD-swap is not yet supported with the experimental block-eigen operator."
   );
  }
  if (static_cast<int>(af_cpp.size()) != m) {
   throw std::runtime_error("af length must equal m for block-eigen BayesC.");
  }

  PreparedBlockEigenOperator prepared = prepare_block_eigen_operator(
   block_input, m, wy_mat
  );
  Rcpp::List diagnostics;
  if (prepared.op.low_rank) {
   diagnostics = Rcpp::List::create(
    Rcpp::Named("blocks") = block_low_rank_diagnostics_to_data_frame(
     prepared.low_rank_diagnostics),
    Rcpp::Named("operator_contract") = "block_low_rank_v1",
    Rcpp::Named("operator_representation") = "low_rank",
    Rcpp::Named("operator_scale_contract") = "general_cross_product",
    Rcpp::Named("eigen_policy") = "cumulative_positive_mass",
    Rcpp::Named("eigen_prop") = eigen_prop,
    Rcpp::Named("build") = block_low_rank_build_metadata(prepared.op.retained)
   );
  } else {
   diagnostics = Rcpp::List::create(
    Rcpp::Named("blocks") = block_eigen_diagnostics_to_data_frame(
     prepared.dense_diagnostics),
    Rcpp::Named("operator_contract") = "block_dense_reconstructed_v1",
    Rcpp::Named("operator_representation") = "dense_reconstructed"
   );
  }
  LDLDFriends ld_swap_friends;
  ld_swap_friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  return BayescOperatorContext<BlockEigenDispatchOperator>(
   std::move(prepared.op),
   std::move(ld_swap_friends),
   diagnostics
  );
 };

 return stblr_cpg_omp_csr_impl(
  wy, ww, yy, b_init, d_init, use_d_init, r_init, use_r_init,
  rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, pi,
  nub, nue, updateB, updateE, updatePi, adjE, n, nit, nburn, nthin,
  pi_prior_a, pi_prior_b, ncores, seed, nchains, keep_chains, chain_seeds,
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves,
  maf_effect_s_prior_scale, estimate_maf_effect_s, maf_effect_s_init,
  maf_effect_s_prior, maf_effect_s_proposal_sd, maf_effect_s_log_h,
  Rcpp::as<std::vector<int>>(convergence_markers), convergence_b,
  convergence_d,
  low_rank_residual_rebuild_every, policy_factory, execution_contract,
  make_block_eigen_operator
 );
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_block_eigen(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> d_init,
  bool use_d_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  std::string ld_prefix,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<double> pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  double adjE,
  std::vector<int> n,
  int nit,
  int nburn,
  int nthin,
  double pi_prior_a,
  double pi_prior_b,
  int ncores,
  int seed,
  int nchains,
  bool keep_chains,
  std::vector<int> chain_seeds,
  bool updateLDswap = false,
  double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8,
  int ld_swap_max_friends = 50,
  int ld_swap_moves = 1,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale = R_NilValue,
  bool estimate_maf_effect_s = false,
  double maf_effect_s_init = 0.0,
  Rcpp::NumericVector maf_effect_s_prior = Rcpp::NumericVector::create(-3.0, 2.0),
  double maf_effect_s_proposal_sd = 0.35,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h = R_NilValue,
  Rcpp::IntegerVector convergence_markers = Rcpp::IntegerVector::create(),
  bool convergence_b = false,
  bool convergence_d = false,
  Rcpp::CharacterVector bed_files = Rcpp::CharacterVector::create(),
  int n_bed = 0,
  Rcpp::List cls = R_NilValue,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::NumericVector af = Rcpp::NumericVector::create(),
  Rcpp::IntegerVector block_start = Rcpp::IntegerVector::create(),
  std::string eigen_filter = "hard_truncate",
  double eigen_tau = 0.01,
  double eigen_eta = 0.0,
  std::string representation = "dense_reconstructed",
  double eigen_prop = 0.995,
  int low_rank_residual_rebuild_every = 100,
  Rcpp::Nullable<Rcpp::List> execution_contract = R_NilValue
) {
 const int expected_tasks = static_cast<int>(wy.size()) * nchains;
 const BlrPhase3ExecutionContract phase3 =
  parse_blr_phase3_execution_contract(execution_contract, expected_tasks, nit);
 return stblr_cpg_omp_csr_block_eigen_with_policy(
  std::move(wy), std::move(ww), std::move(yy), std::move(b_init),
  std::move(d_init), use_d_init, std::move(r_init), use_r_init,
  rebuild_r_before_updateE, std::move(ld_prefix), std::move(B), std::move(E),
  std::move(ssb_prior), std::move(sse_prior), std::move(pi), nub, nue,
  updateB, updateE, updatePi, adjE, std::move(n), nit, nburn, nthin,
  pi_prior_a, pi_prior_b, ncores, seed, nchains, keep_chains,
  std::move(chain_seeds), updateLDswap, ld_swap_prob, ld_swap_r2,
  ld_swap_max_friends, ld_swap_moves, maf_effect_s_prior_scale,
  estimate_maf_effect_s, maf_effect_s_init, maf_effect_s_prior,
  maf_effect_s_proposal_sd, maf_effect_s_log_h, convergence_markers,
  convergence_b, convergence_d, bed_files, n_bed, cls, rows, af,
  block_start, std::move(eigen_filter), eigen_tau, eigen_eta,
  std::move(representation), eigen_prop, low_rank_residual_rebuild_every,
  nullptr, phase3
 );
}

// Development-only benchmark of the real retained-low-rank BayesC hot path.
// Factor construction is intentionally outside every measured region.
// [[Rcpp::export]]
Rcpp::List stblr_low_rank_bayesc_hot_path_benchmark_internal(
  int markers = 1000, int rank = 250, int repetitions = 7,
  int warmup = 2, int seed = 7301
) {
 BlockLowRankOperator op =
  sblr::core::make_block_low_rank_hot_path_fixture(markers, rank);
 const arma::rowvec& diagonal = op.diag();
 arma::rowvec score(static_cast<arma::uword>(markers), arma::fill::zeros);
 const sblr::core::BlockLowRankBlock& block = op.blocks[0];
 for (int marker = 0; marker < markers; ++marker) {
  const float* factor_ptr = block.factor.data() +
   static_cast<std::size_t>(marker) * static_cast<std::size_t>(rank);
  for (int k = 0; k < rank; ++k) {
   score(static_cast<arma::uword>(marker)) +=
    static_cast<double>(factor_ptr[k]) *
    block.transformed_score(0, static_cast<arma::uword>(k));
  }
 }
 arma::rowvec zero_effects(static_cast<arma::uword>(markers), arma::fill::zeros);
 arma::rowvec base_residual;
 op.rebuild(0, score, zero_effects, base_residual);
 volatile double sink = 0.0;

 auto measure = [&](const std::function<void()>& operation,
                    double divisor = 1.0) {
  if (repetitions <= 0 || warmup < 0)
   throw std::runtime_error("benchmark repetitions/warmup are invalid.");
  for (int i = 0; i < warmup; ++i) operation();
  std::vector<double> seconds(static_cast<std::size_t>(repetitions));
  for (int i = 0; i < repetitions; ++i) {
   const auto start = std::chrono::steady_clock::now();
   operation();
   const auto stop = std::chrono::steady_clock::now();
   seconds[static_cast<std::size_t>(i)] =
    std::chrono::duration<double>(stop - start).count() / divisor;
  }
  std::sort(seconds.begin(), seconds.end());
  return Rcpp::NumericVector::create(
   Rcpp::Named("median_seconds") =
    seconds[static_cast<std::size_t>(repetitions / 2)],
   Rcpp::Named("minimum_seconds") = seconds.front()
  );
 };

 const int coordinate_batches = 10000;
 Rcpp::List timings;
 timings["corrected_rhs_1000"] = measure([&]() {
  double value = 0.0;
  for (int call = 0; call < 1000; ++call)
   value += op.corrected_rhs(call % markers, 0.0, base_residual);
  sink = value;
 });
 timings["apply_difference_1000"] = measure([&]() {
  arma::rowvec residual = base_residual;
  for (int call = 0; call < 1000; ++call)
   op.apply_difference(call % markers, (call & 1) ? 1e-8 : -1e-8, residual);
  sink = residual(0);
 });
 timings["fitted_quadratic"] = measure([&]() {
  double value = 0.0;
  for (int call = 0; call < coordinate_batches; ++call)
   value += op.fitted_quadratic(0, zero_effects, score, base_residual);
  sink = value;
 }, static_cast<double>(coordinate_batches));
 arma::rowvec reference_effects(static_cast<arma::uword>(markers));
 for (int marker = 0; marker < markers; ++marker)
  reference_effects(static_cast<arma::uword>(marker)) =
   0.001 * std::sin(0.1 * static_cast<double>(marker + 1));
 timings["quadratic_form_reference"] = measure([&]() {
  sink = op.quadratic_form(reference_effects);
 });
 timings["residual_rebuild"] = measure([&]() {
  arma::rowvec residual;
  op.rebuild(0, score, reference_effects, residual);
  sink = residual(0);
 });

 auto bayesc_sweep = [&]() {
  arma::rowvec effects = zero_effects;
  arma::rowvec residual = base_residual;
  arma::Row<int> state(static_cast<arma::uword>(markers), arma::fill::zeros);
  std::vector<double> pi = {0.5, 0.5};
  std::mt19937 generator(static_cast<unsigned int>(seed));
  for (int marker = 0; marker < markers; ++marker) {
   sampleBetaC_ST_csr_unscaled(
    marker, pi, 0.2, 1.0, diagonal, residual, effects, state, op, generator
   );
  }
  sink = effects(0) + residual(0);
 };
 timings["bayesc_marker_sweep"] = measure(bayesc_sweep);
 timings["bayesc_gibbs_iteration"] = measure([&]() {
  arma::rowvec effects = zero_effects;
  arma::rowvec residual = base_residual;
  arma::Row<int> state(static_cast<arma::uword>(markers), arma::fill::zeros);
  std::vector<double> pi = {0.5, 0.5};
  double vb = 0.2;
  double ve = 1.0;
  std::mt19937 generator(static_cast<unsigned int>(seed));
  for (int marker = 0; marker < markers; ++marker) {
   sampleBetaC_ST_csr_unscaled(
    marker, pi, vb, ve, diagonal, residual, effects, state, op, generator
   );
  }
  sampleB_ST_csr_unscaled(markers, 4.0, vb, effects, state, 0.1, generator);
  sampleE_ST_operator(
   op, 0, 4.0, ve, effects, score, residual, 0.1,
   2.0 * op.transformed_score_norm_squared(0) + 1.0,
   markers, generator
  );
  samplePi_ST(state, pi, 2.0, 2.0, generator);
  sink = op.fitted_quadratic(0, effects, score, residual) /
   static_cast<double>(markers) + vb + ve + pi[1];
 });
 timings["bayesc_gibbs_iteration_direct_reference"] = measure([&]() {
  arma::rowvec effects = zero_effects;
  arma::rowvec residual = base_residual;
  arma::Row<int> state(static_cast<arma::uword>(markers), arma::fill::zeros);
  std::vector<double> pi = {0.5, 0.5};
  double vb = 0.2;
  double ve = 1.0;
  std::mt19937 generator(static_cast<unsigned int>(seed));
  for (int marker = 0; marker < markers; ++marker) {
   sampleBetaC_ST_csr_unscaled(
    marker, pi, vb, ve, diagonal, residual, effects, state, op, generator
   );
  }
  sampleB_ST_csr_unscaled(markers, 4.0, vb, effects, state, 0.1, generator);
  sampleE_ST_operator(
   op, 0, 4.0, ve, effects, score, residual, 0.1,
   2.0 * op.transformed_score_norm_squared(0) + 1.0,
   markers, generator
  );
  samplePi_ST(state, pi, 2.0, 2.0, generator);
  sink = op.quadratic_form(effects) / static_cast<double>(markers) +
   vb + ve + pi[1];
 });

 return Rcpp::List::create(
  Rcpp::Named("markers") = markers,
  Rcpp::Named("rank") = rank,
  Rcpp::Named("repetitions") = repetitions,
  Rcpp::Named("warmup") = warmup,
  Rcpp::Named("updates_enabled") =
   "marker effects/states, effect variance, residual variance, pi, genetic variance",
  Rcpp::Named("timings") = timings,
  Rcpp::Named("sink") = static_cast<double>(sink)
 );
}
