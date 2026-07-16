// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "st_chain_utils.h"
#include "st_csr_common.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

// =============================================================================
// STBLR summary-stat CSR with group-level annotation priors
// =============================================================================
//
// Exported function:
//
//   stblr_cpg_omp_csr_group_annot(...)
//
// Purpose:
//   Group-level annotation prior model for summary-statistic CSR BayesC.
//
// Model:
//   Each marker i belongs to exactly one group g_i in {0, ..., G-1}.
//
//   pi_g,t ~ Beta(pi_group_prior_a[g], pi_group_prior_b[g])
//   d_i,t | g_i ~ Bernoulli(pi_g_i,t)
//
//   b_i,t | d_i,t = 1, g_i ~ N(0, vb_t * lambda_g_i,t)
//
// where lambda_g is a group-specific variance multiplier.
//
// Options:
//   updatePi = TRUE:
//     updates group-specific pi_g using marker inclusion counts by group.
//
//   updateGroupVb = TRUE:
//     updates group-specific variance multipliers lambda_g by inverse-chi-square
//     style conditional updates, normalized to have weighted mean 1 across
//     groups. This keeps vb_t interpretable as the global scale.
//
//   updateB = TRUE:
//     updates global vb_t conditional on current lambda_g.
//
// This is intentionally a full-sweep CSR sampler. No sparse scheduling.
//
// Return layout, aligned with current CSR BLR conventions plus group outputs:
//
//   0  bm
//   1  dm
//   2  wy
//   3  r
//   4  b
//   5  d
//   6  marker index
//   7  vbs
//   8  vgs
//   9  ves
//   10 covb
//   11 covg
//   12 cove
//   13 vb
//   14 vg
//   15 ve
//   16 final global pi, c(pi0, pi1), where pi1 is marker-weighted mean pi_g
//   17 posterior mean global pi, c(pi0, pi1)
//   18 diagnostics: reserved/log-CPO placeholders + runtime, length 4
//   19 diagnostics: nsamples, n, length 2
//   20 vle trace
//   21 vld trace
//   22 posterior mean group pi, flattened ngroup x nt by trait slot
//   23 posterior mean group variance multiplier, flattened ngroup x nt by trait slot
//   24 posterior mean group inclusion counts, length ngroup by trait slot
//   25 group sizes, length ngroup by trait slot
//
// R indexing:
//   raw_fit[[23]] = group_pi
//   raw_fit[[24]] = group_vb_multiplier
//   raw_fit[[25]] = group_nincluded
//   raw_fit[[26]] = group_size
//
// Neutral model check:
//   If all markers have group = 1 in R, ngroup = 1, and group_vb_multiplier = 1,
//   this reduces to the ordinary CSR BayesC sampler, apart from harmless extra
//   bookkeeping.
// =============================================================================

inline double clamp_prob_group(double x) {
 if (!std::isfinite(x)) {
  throw std::runtime_error("clamp_prob_group: probability is NaN/Inf.");
 }
 return std::min(std::max(x, 1e-300), 1.0 - 1e-12);
}

inline double computeLE_ST_csr_group(
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

inline void sampleBetaC_ST_csr_group(
  int i,
  double pi1_i,
  double vb_t,
  double vb_mult_i,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const STLDCSR& ld,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double wi = ww(iu);

 if (!std::isfinite(wi) || wi <= 0.0) {
  throw std::runtime_error("sampleBetaC_ST_csr_group: invalid ww value.");
 }
 if (!std::isfinite(vb_t) || vb_t <= 0.0) {
  throw std::runtime_error("sampleBetaC_ST_csr_group: invalid global vb.");
 }
 if (!std::isfinite(vb_mult_i) || vb_mult_i <= 0.0) {
  throw std::runtime_error("sampleBetaC_ST_csr_group: invalid group vb multiplier.");
 }

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm01(0.0, 1.0);

 pi1_i = clamp_prob_group(pi1_i);
 const double pi0_i = std::max(1.0 - pi1_i, 1e-300);
 const double vbi = std::max(vb_t * vb_mult_i, 1e-12);
 const double vei_safe = std::max(vei_i, 1e-300);

 const double score = r(iu) + wi * b(iu);
 const double denom = std::max(vei_safe + wi * vbi, 1e-300);

 const double logBF =
  0.5 * std::log(vei_safe / denom) +
  0.5 * score * score * vbi / (vei_safe * denom);

 const double logp1 = std::log(pi1_i) + logBF;
 const double logp0 = std::log(pi0_i);
 const double delta_log = logp0 - logp1;

 double p1 = 0.0;
 if (delta_log > 35.0) p1 = 0.0;
 else if (delta_log < -35.0) p1 = 1.0;
 else p1 = 1.0 / (1.0 + std::exp(delta_log));

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
  r(iu) -= wi * diff;

  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   r(static_cast<arma::uword>(j)) -=
    static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
  }
 }

 b(iu) = b_new;
 d(iu) = di;
}

inline void sampleB_ST_csr_group(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  const arma::Row<int>& group,
  const arma::rowvec& group_vb_multiplier,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb_scaled = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (d(iu) > 0) {
   const int g = group(iu);
   const double mult = group_vb_multiplier(static_cast<arma::uword>(g));

   if (!std::isfinite(mult) || mult <= 0.0) {
    throw std::runtime_error("sampleB_ST_csr_group: invalid group vb multiplier.");
   }

   ssb_scaled += b(iu) * b(iu) / mult;
   dfb += 1.0;
  }
 }

 const double scale = ssb_scaled + nub * ssb_prior;
 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleB_ST_csr_group: invalid scale.");
 }

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 vb = std::max(scale / chi2, 1e-12);
}

inline void samplePiGroups_ST_csr_group(
  const arma::Row<int>& d,
  const arma::Row<int>& group,
  arma::rowvec& group_pi,
  const arma::rowvec& prior_a,
  const arma::rowvec& prior_b,
  int ngroup,
  std::mt19937& gen
) {
 arma::rowvec c1 = prior_a;
 arma::rowvec c0 = prior_b;

 if (static_cast<int>(prior_a.n_elem) != ngroup ||
     static_cast<int>(prior_b.n_elem) != ngroup ||
     static_cast<int>(group_pi.n_elem) != ngroup) {
  throw std::runtime_error("samplePiGroups_ST_csr_group: invalid group prior dimensions.");
 }

 for (arma::uword i = 0; i < d.n_elem; ++i) {
  const int g = group(i);
  if (g < 0 || g >= ngroup) {
   throw std::runtime_error("samplePiGroups_ST_csr_group: group index out of range.");
  }

  if (d(i) > 0) c1(static_cast<arma::uword>(g)) += 1.0;
  else          c0(static_cast<arma::uword>(g)) += 1.0;
 }

 for (int g = 0; g < ngroup; ++g) {
  const arma::uword gu = static_cast<arma::uword>(g);
  if (!std::isfinite(c1(gu)) || c1(gu) <= 0.0 ||
      !std::isfinite(c0(gu)) || c0(gu) <= 0.0) {
      throw std::runtime_error("samplePiGroups_ST_csr_group: invalid Beta posterior shape.");
  }

  std::gamma_distribution<double> rg0(c0(gu), 1.0);
  std::gamma_distribution<double> rg1(c1(gu), 1.0);

  const double g0 = std::max(rg0(gen), 1e-300);
  const double g1 = std::max(rg1(gen), 1e-300);
  group_pi(gu) = g1 / (g0 + g1);
 }
}

inline void sampleGroupVbMultipliers_ST_csr_group(
  int m,
  int ngroup,
  double nub_group,
  double ssb_group_prior,
  double vb_global,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  const arma::Row<int>& group,
  const arma::rowvec& group_size,
  arma::rowvec& group_vb_multiplier,
  bool normalize_group_vb,
  std::mt19937& gen
) {
 if (!std::isfinite(vb_global) || vb_global <= 0.0) {
  throw std::runtime_error("sampleGroupVbMultipliers_ST_csr_group: invalid global vb.");
 }
 if (!std::isfinite(nub_group) || nub_group <= 0.0) {
  throw std::runtime_error("sampleGroupVbMultipliers_ST_csr_group: nub_group must be positive.");
 }
 if (!std::isfinite(ssb_group_prior) || ssb_group_prior <= 0.0) {
  throw std::runtime_error("sampleGroupVbMultipliers_ST_csr_group: ssb_group_prior must be positive.");
 }

 arma::rowvec ssb(ngroup, arma::fill::zeros);
 arma::rowvec df(ngroup, arma::fill::zeros);

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (d(iu) > 0) {
   const int g = group(iu);
   if (g < 0 || g >= ngroup) {
    throw std::runtime_error("sampleGroupVbMultipliers_ST_csr_group: group index out of range.");
   }
   const arma::uword gu = static_cast<arma::uword>(g);
   ssb(gu) += b(iu) * b(iu) / vb_global;
   df(gu) += 1.0;
  }
 }

 for (int g = 0; g < ngroup; ++g) {
  const arma::uword gu = static_cast<arma::uword>(g);
  const double scale = ssb(gu) + nub_group * ssb_group_prior;
  std::chi_squared_distribution<double> rchisq(df(gu) + nub_group);
  const double chi2 = std::max(rchisq(gen), 1e-300);
  group_vb_multiplier(gu) = std::max(scale / chi2, 1e-12);
 }

 if (normalize_group_vb) {
  double weighted_mean = 0.0;
  double total_size = 0.0;

  for (int g = 0; g < ngroup; ++g) {
   const arma::uword gu = static_cast<arma::uword>(g);
   weighted_mean += group_size(gu) * group_vb_multiplier(gu);
   total_size += group_size(gu);
  }

  if (total_size > 0.0) weighted_mean /= total_size;

  if (std::isfinite(weighted_mean) && weighted_mean > 0.0) {
   group_vb_multiplier /= weighted_mean;
  }
 }
}

inline double marker_weighted_mean_group_value(
  const arma::rowvec& x,
  const arma::rowvec& group_size
) {
 if (x.n_elem != group_size.n_elem) {
  throw std::runtime_error("marker_weighted_mean_group_value: dimension mismatch.");
 }

 double s = 0.0;
 double n = 0.0;

 for (arma::uword g = 0; g < x.n_elem; ++g) {
  s += x(g) * group_size(g);
  n += group_size(g);
 }

 if (n <= 0.0) return NA_REAL;
 return s / n;
}

inline arma::rowvec count_group_inclusions(
  const arma::Row<int>& d,
  const arma::Row<int>& group,
  int ngroup
) {
 arma::rowvec out(ngroup, arma::fill::zeros);

 for (arma::uword i = 0; i < d.n_elem; ++i) {
  if (d(i) > 0) {
   const int g = group(i);
   if (g < 0 || g >= ngroup) {
    throw std::runtime_error("count_group_inclusions: group index out of range.");
   }
   out(static_cast<arma::uword>(g)) += 1.0;
  }
 }

 return out;
}

struct LDLDFriendsGroup {
 std::vector<uint64_t> ptr;
 std::vector<int> idx;
 std::vector<double> r2;
};

inline LDLDFriendsGroup build_ld_swap_friends_st_csr_group(
  int m,
  const STLDCSR& ld,
  const std::vector<double>& xx,
  double min_r2,
  int max_friends
) {
 std::vector<std::vector<std::pair<int, double>>> rows(static_cast<std::size_t>(m));

 for (int i = 0; i < m; ++i) {
  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   if (j <= i) continue;

   const double denom = xx[static_cast<std::size_t>(i)] * xx[static_cast<std::size_t>(j)];
   if (!std::isfinite(denom) || denom <= 0.0) continue;

   const double xij = static_cast<double>(ld.xij[static_cast<std::size_t>(p)]);
   const double r2 = (xij * xij) / denom;
   if (!std::isfinite(r2) || r2 < min_r2) continue;

   rows[static_cast<std::size_t>(i)].push_back(std::make_pair(j, r2));
   rows[static_cast<std::size_t>(j)].push_back(std::make_pair(i, r2));
  }
 }

 LDLDFriendsGroup friends;
 friends.ptr.resize(static_cast<std::size_t>(m) + 1);
 friends.ptr[0] = 0;

 for (int i = 0; i < m; ++i) {
  std::vector<std::pair<int, double>>& row = rows[static_cast<std::size_t>(i)];
  std::sort(row.begin(), row.end(),
            [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
             if (a.second == b.second) return a.first < b.first;
             return a.second > b.second;
            });

  std::vector<std::pair<int, double>> unique_row;
  unique_row.reserve(row.size());
  for (std::size_t k = 0; k < row.size(); ++k) {
   if (!unique_row.empty() && unique_row.back().first == row[k].first) continue;
   unique_row.push_back(row[k]);
  }
  if (static_cast<int>(unique_row.size()) > max_friends) {
   unique_row.resize(static_cast<std::size_t>(max_friends));
  }
  row.swap(unique_row);
  friends.ptr[static_cast<std::size_t>(i + 1)] =
   friends.ptr[static_cast<std::size_t>(i)] + row.size();
 }

 const uint64_t nfriend = friends.ptr[static_cast<std::size_t>(m)];
 friends.idx.resize(static_cast<std::size_t>(nfriend));
 friends.r2.resize(static_cast<std::size_t>(nfriend));
 for (int i = 0; i < m; ++i) {
  const uint64_t offset = friends.ptr[static_cast<std::size_t>(i)];
  const std::vector<std::pair<int, double>>& row = rows[static_cast<std::size_t>(i)];
  for (std::size_t k = 0; k < row.size(); ++k) {
   friends.idx[static_cast<std::size_t>(offset + k)] = row[k].first;
   friends.r2[static_cast<std::size_t>(offset + k)] = row[k].second;
  }
 }

 return friends;
}

inline void set_marker_effect_st_csr_group(
  int i,
  double b_new,
  int d_new,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const STLDCSR& ld
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  r(iu) -= ww(iu) * diff;
  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];
  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   r(static_cast<arma::uword>(j)) -=
    static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
  }
 }

 b(iu) = b_new;
 d(iu) = d_new;
}

inline double residual_sse_st_csr_group(
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

inline int count_excluded_ld_friends_group(
  int i,
  const arma::Row<int>& d,
  const LDLDFriendsGroup& friends
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

inline int collect_ld_swap_candidates_group(
  int m,
  const arma::Row<int>& d,
  const LDLDFriendsGroup& friends,
  std::vector<int>& candidates,
  std::vector<int>& n_excluded
) {
 candidates.clear();
 n_excluded.clear();
 for (int i = 0; i < m; ++i) {
  if (d(static_cast<arma::uword>(i)) <= 0) continue;
  const int nf = count_excluded_ld_friends_group(i, d, friends);
  if (nf > 0) {
   candidates.push_back(i);
   n_excluded.push_back(nf);
  }
 }
 return static_cast<int>(candidates.size());
}

inline double log_group_prior_ratio_group(
  double b_move,
  int j,
  int k,
  double vb_t,
  const arma::Row<int>& group,
  const arma::rowvec& group_pi,
  const arma::rowvec& group_vb_multiplier
) {
 const int gj = group(static_cast<arma::uword>(j));
 const int gk = group(static_cast<arma::uword>(k));
 const arma::uword gju = static_cast<arma::uword>(gj);
 const arma::uword gku = static_cast<arma::uword>(gk);

 const double pi_j = clamp_prob_group(group_pi(gju));
 const double pi_k = clamp_prob_group(group_pi(gku));
 double log_ratio =
  std::log(pi_k) + std::log(std::max(1.0 - pi_j, 1e-300)) -
  std::log(pi_j) - std::log(std::max(1.0 - pi_k, 1e-300));

 const double mult_j = group_vb_multiplier(gju);
 const double mult_k = group_vb_multiplier(gku);
 if (!std::isfinite(vb_t) || vb_t <= 0.0 ||
     !std::isfinite(mult_j) || mult_j <= 0.0 ||
     !std::isfinite(mult_k) || mult_k <= 0.0) {
  return -std::numeric_limits<double>::infinity();
 }
 const double sigma2_j = std::max(vb_t * mult_j, 1e-300);
 const double sigma2_k = std::max(vb_t * mult_k, 1e-300);
 log_ratio +=
  0.5 * (std::log(sigma2_j) - std::log(sigma2_k)) -
  0.5 * b_move * b_move * (1.0 / sigma2_k - 1.0 / sigma2_j);

 return log_ratio;
}

inline bool attempt_ld_swap_st_csr_group(
  int m,
  double vei,
  double yy,
  double vb_t,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  const arma::Row<int>& group,
  const arma::rowvec& group_pi,
  const arma::rowvec& group_vb_multiplier,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const STLDCSR& ld,
  const LDLDFriendsGroup& friends,
  std::mt19937& gen
) {
 if (!std::isfinite(vei) || vei <= 0.0) return false;

 std::vector<int> candidates;
 std::vector<int> n_excluded;
 const int n_candidates = collect_ld_swap_candidates_group(m, d, friends, candidates, n_excluded);
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

 const double sse_old = residual_sse_st_csr_group(m, b, wy, r, yy);
 if (!std::isfinite(sse_old)) return false;

 const arma::rowvec r_old = r;
 set_marker_effect_st_csr_group(j, 0.0, 0, ww, r, b, d, ld);
 set_marker_effect_st_csr_group(k, b_j_old, 1, ww, r, b, d, ld);

 const double sse_new = residual_sse_st_csr_group(m, b, wy, r, yy);
 bool accept = false;

 if (std::isfinite(sse_new)) {
  std::vector<int> reverse_candidates;
  std::vector<int> reverse_n_excluded;
  const int n_reverse_candidates =
   collect_ld_swap_candidates_group(m, d, friends, reverse_candidates, reverse_n_excluded);

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
   const double log_prior_ratio = log_group_prior_ratio_group(
    b_j_old, j, k, vb_t, group, group_pi, group_vb_multiplier
   );
   const double log_alpha =
    -0.5 * (sse_new - sse_old) / vei +
    log_prior_ratio +
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

std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_group_annot_single(
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
  std::vector<int> group_index,
  int ngroup,
  std::vector<std::vector<double>> group_pi_init,
  std::vector<double> pi_group_prior_a,
  std::vector<double> pi_group_prior_b,
  std::vector<std::vector<double>> group_vb_multiplier_init,
  bool updateGroupVb,
  double nub_group,
  double ssb_group_prior,
  bool normalize_group_vb,
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
  int ncores,
  int seed,
  bool updateLDswap = false,
  double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8,
  int ld_swap_max_friends = 50,
  int ld_swap_moves = 1
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) throw std::runtime_error("stblr_cpg_omp_csr_group_annot: nt must be positive.");
 const int m = static_cast<int>(wy[0].size());
 if (m <= 0) throw std::runtime_error("stblr_cpg_omp_csr_group_annot: m must be positive.");

 if (nit <= 0 || nburn < 0 || nthin <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: invalid nit/nburn/nthin.");
 }

 if (ngroup <= 0) throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ngroup must be positive.");
 if (static_cast<int>(group_index.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: group_index must have length m.");
 }
 if (!std::isfinite(ld_swap_prob) || ld_swap_prob < 0.0 || ld_swap_prob > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ld_swap_prob must be in [0, 1].");
 }
 if (!std::isfinite(ld_swap_r2) || ld_swap_r2 < 0.0 || ld_swap_r2 > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ld_swap_r2 must be in [0, 1].");
 }
 if (ld_swap_max_friends <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ld_swap_max_friends must be positive.");
 }
 if (ld_swap_moves < 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ld_swap_moves must be non-negative.");
 }

 if (pi.size() != 2) throw std::runtime_error("stblr_cpg_omp_csr_group_annot: pi must have length 2.");
 if ((int)ww.size() != nt || (int)b_init.size() != nt || (int)yy.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: inconsistent trait dimensions.");
 }
 if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: priors must be nt x nt.");
 }
 if ((int)B.n_rows != nt || (int)B.n_cols != nt) throw std::runtime_error("B must be nt x nt.");
 if ((int)E.n_rows != nt || (int)E.n_cols != nt) throw std::runtime_error("E must be nt x nt.");

 if (static_cast<int>(group_pi_init.size()) != nt ||
     static_cast<int>(group_vb_multiplier_init.size()) != nt) {
  throw std::runtime_error("group_pi_init and group_vb_multiplier_init must have length nt.");
 }
 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(group_pi_init[static_cast<std::size_t>(t)].size()) != ngroup ||
      static_cast<int>(group_vb_multiplier_init[static_cast<std::size_t>(t)].size()) != ngroup) {
   throw std::runtime_error("group initial vectors must have length ngroup for each trait.");
  }
 }
 if (static_cast<int>(pi_group_prior_a.size()) != ngroup ||
     static_cast<int>(pi_group_prior_b.size()) != ngroup) {
  throw std::runtime_error("pi_group_prior_a and pi_group_prior_b must have length ngroup.");
 }

 arma::Row<int> group(m, arma::fill::zeros);
 arma::rowvec group_size(ngroup, arma::fill::zeros);

 for (int i = 0; i < m; ++i) {
  const int g0 = group_index[static_cast<std::size_t>(i)];
  if (g0 < 0 || g0 >= ngroup) {
   throw std::runtime_error("group_index must be 0-based and in [0, ngroup-1].");
  }
  group(static_cast<arma::uword>(i)) = g0;
  group_size(static_cast<arma::uword>(g0)) += 1.0;
 }
 for (int g = 0; g < ngroup; ++g) {
  if (group_size(static_cast<arma::uword>(g)) <= 0.0) {
   throw std::runtime_error("Every group must contain at least one marker.");
  }
 }

 arma::rowvec prior_a(ngroup, arma::fill::zeros);
 arma::rowvec prior_b(ngroup, arma::fill::zeros);
 for (int g = 0; g < ngroup; ++g) {
  const double a = pi_group_prior_a[static_cast<std::size_t>(g)];
  const double b = pi_group_prior_b[static_cast<std::size_t>(g)];
  if (!std::isfinite(a) || a <= 0.0 || !std::isfinite(b) || b <= 0.0) {
   throw std::runtime_error("group Beta prior shapes must be finite and positive.");
  }
  prior_a(static_cast<arma::uword>(g)) = a;
  prior_b(static_cast<arma::uword>(g)) = b;
 }

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m || (int)ww[t].size() != m || (int)b_init[t].size() != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_group_annot: inconsistent marker dimensions.");
  }
 }
 if (use_r_init) {
  if (static_cast<int>(r_init.size()) != nt) throw std::runtime_error("r_init must have length nt.");
  for (int t = 0; t < nt; ++t) if (static_cast<int>(r_init[t].size()) != m) throw std::runtime_error("r_init[t] must have length m.");
 }
 if (use_d_init) {
  if (static_cast<int>(d_init.size()) != nt) throw std::runtime_error("d_init must have length nt.");
  for (int t = 0; t < nt; ++t) if (static_cast<int>(d_init[t].size()) != m) throw std::runtime_error("d_init[t] must have length m.");
 }

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);
 arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
 arma::vec yy_vec(nt, arma::fill::zeros);
 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  yy_vec(static_cast<arma::uword>(t)) = yy[t];
  for (int i = 0; i < m; ++i) {
   wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
   ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
   b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = b_init[t][i];
  }
  if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("priors must be nt x nt.");
  }
  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
   sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
  }
 }

 for (int t = 1; t < nt; ++t) {
  if (n[t] != n[0]) {
   throw std::runtime_error("stblr_cpg_omp_csr_group_annot: shared LD scaling assumes equal n across traits.");
  }
 }
 for (int t = 1; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   const double w0 = ww_mat(0, static_cast<arma::uword>(i));
   const double wt = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   const double tol = 1e-8 * std::max(1.0, std::abs(w0));
   if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
    throw std::runtime_error("ww contains invalid value before LD pre-scaling.");
   }
   if (std::abs(w0 - wt) > tol) {
    throw std::runtime_error("ww differs across traits; pre-scaled shared ST LD is invalid.");
   }
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
 for (int i = 0; i < m; ++i) {
  const double wi = ww_mat(0, static_cast<arma::uword>(i));
  if (!std::isfinite(wi) || wi <= 0.0) throw std::runtime_error("ww contains invalid value in trait 0.");
  xx[static_cast<std::size_t>(i)] = wi;
 }
 STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);

 LDLDFriendsGroup ld_swap_friends;
 if (updateLDswap) {
  ld_swap_friends = build_ld_swap_friends_st_csr_group(
   m,
   ld,
   xx,
   ld_swap_r2,
   ld_swap_max_friends
  );
 } else {
  ld_swap_friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
 }

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
 std::sort(order.begin(), order.end(), [&](int a, int b) {
  return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
 });

#define SBLR_CSR_GROUP_BAYESC_CORE_IMPL_TRANSLATION_UNIT 1
#include "blr_csr_group_bayesc_core_impl.h"
#undef SBLR_CSR_GROUP_BAYESC_CORE_IMPL_TRANSLATION_UNIT
}

static Rcpp::NumericMatrix cpg_raw_marker_matrix(
 const std::vector<std::vector<std::vector<double>>>& raw,
 std::size_t slot,
 int m,
 int nt
) {
 Rcpp::NumericMatrix out(m, nt);
 if (slot >= raw.size()) return out;
 for (int t = 0; t < nt; ++t) {
  const std::vector<double>& x = raw[slot][static_cast<std::size_t>(t)];
  for (int i = 0; i < m && static_cast<std::size_t>(i) < x.size(); ++i) out(i, t) = x[i];
 }
 return out;
}

static Rcpp::NumericMatrix cpg_raw_trace_matrix(
 const std::vector<std::vector<std::vector<double>>>& raw,
 std::size_t slot,
 int n_trace,
 int nt
) {
 Rcpp::NumericMatrix out(n_trace, nt);
 if (slot >= raw.size()) return out;
 for (int t = 0; t < nt; ++t) {
  const std::vector<double>& x = raw[slot][static_cast<std::size_t>(t)];
  for (int it = 0; it < n_trace && static_cast<std::size_t>(it) < x.size(); ++it) out(it, t) = x[it];
 }
 return out;
}

static Rcpp::NumericMatrix cpg_raw_trait_matrix(
 const std::vector<std::vector<std::vector<double>>>& raw,
 std::size_t slot,
 int nt,
 int ncol
) {
 Rcpp::NumericMatrix out(nt, ncol);
 if (slot >= raw.size()) return out;
 for (int t = 0; t < nt; ++t) {
  const std::vector<double>& x = raw[slot][static_cast<std::size_t>(t)];
  for (int j = 0; j < ncol && static_cast<std::size_t>(j) < x.size(); ++j) out(t, j) = x[j];
 }
 return out;
}

static Rcpp::NumericMatrix cpg_raw_group_matrix(
 const std::vector<std::vector<std::vector<double>>>& raw,
 std::size_t slot,
 int ngroup,
 int nt
) {
 Rcpp::NumericMatrix out(ngroup, nt);
 if (slot >= raw.size()) return out;
 for (int t = 0; t < nt; ++t) {
  const std::vector<double>& x = raw[slot][static_cast<std::size_t>(t)];
  for (int g = 0; g < ngroup && static_cast<std::size_t>(g) < x.size(); ++g) out(g, t) = x[g];
 }
 return out;
}

static Rcpp::List cpg_disabled_selection() {
 return Rcpp::List::create(
  Rcpp::_["enabled"] = false,
  Rcpp::_["fixed"] = false,
  Rcpp::_["scale"] = "standardized_genotype_effect",
  Rcpp::_["trace"] = R_NilValue,
  Rcpp::_["mean"] = R_NilValue,
  Rcpp::_["sd"] = R_NilValue,
  Rcpp::_["min"] = R_NilValue,
  Rcpp::_["max"] = R_NilValue,
  Rcpp::_["acceptance"] = R_NilValue
 );
}

static Rcpp::List cpg_group_chains_raw_v1(
 const std::vector<std::vector<std::vector<double>>>& raw,
 int m,
 int nt,
 int ngroup,
 int nchains
) {
 if (raw.size() <= 39) return Rcpp::List::create();
 Rcpp::List traits(nt);
 for (int t = 0; t < nt; ++t) {
  Rcpp::List chains(nchains);
  for (int cc = 0; cc < nchains; ++cc) {
   Rcpp::NumericVector dm(m), bm(m);
   Rcpp::NumericVector group_pi(ngroup), group_vb(ngroup), group_nincluded(ngroup);
   for (int i = 0; i < m; ++i) {
    const std::size_t idx = static_cast<std::size_t>(cc * m + i);
    if (idx < raw[34][t].size()) dm[i] = raw[34][t][idx];
    if (idx < raw[35][t].size()) bm[i] = raw[35][t][idx];
   }
   for (int g = 0; g < ngroup; ++g) {
    const std::size_t idx = static_cast<std::size_t>(cc * ngroup + g);
    if (idx < raw[37][t].size()) group_pi[g] = raw[37][t][idx];
    if (idx < raw[38][t].size()) group_vb[g] = raw[38][t][idx];
    if (idx < raw[39][t].size()) group_nincluded[g] = raw[39][t][idx];
   }
   Rcpp::NumericMatrix ld(1, 3);
   const std::size_t didx = static_cast<std::size_t>(cc * 4);
   if (didx + 2 < raw[36][t].size()) {
    ld(0, 0) = raw[36][t][didx];
    ld(0, 1) = raw[36][t][didx + 1];
    ld(0, 2) = raw[36][t][didx + 2];
   }
   chains[cc] = Rcpp::List::create(
    Rcpp::_["marker"] = Rcpp::List::create(Rcpp::_["bm"] = bm, Rcpp::_["dm"] = dm),
    Rcpp::_["group"] = Rcpp::List::create(
     Rcpp::_["pi_mean"] = group_pi,
     Rcpp::_["vb_multiplier_mean"] = group_vb,
     Rcpp::_["n_included_mean"] = group_nincluded
    ),
    Rcpp::_["diagnostics"] = Rcpp::List::create(Rcpp::_["ld_swap"] = ld)
   );
  }
  traits[t] = chains;
 }
 return traits;
}

static Rcpp::List cpg_group_raw_v1(
 const std::vector<std::vector<std::vector<double>>>& raw,
 const std::vector<int>& group_index,
 bool updateLDswap,
 int m,
 int nt,
 int ngroup,
 int nit,
 int nburn,
 int nthin,
 int nchains,
 bool keep_chains
) {
 const int n_trace = nit + nburn;
 Rcpp::NumericVector nsamples(nt), n_used(nt), log_cpo(nt), mean_log_cpo(nt), seconds_mean(nt), seconds_max(nt);
 for (int t = 0; t < nt; ++t) {
  if (raw.size() > 19 && raw[19][t].size() > 0) nsamples[t] = raw[19][t][0];
  if (raw.size() > 19 && raw[19][t].size() > 1) n_used[t] = raw[19][t][1];
  if (raw.size() > 18 && raw[18][t].size() > 0) log_cpo[t] = raw[18][t][0];
  if (raw.size() > 18 && raw[18][t].size() > 1) mean_log_cpo[t] = raw[18][t][1];
  if (raw.size() > 18 && raw[18][t].size() > 2) seconds_mean[t] = raw[18][t][2];
  if (raw.size() > 18 && raw[18][t].size() > 3) seconds_max[t] = raw[18][t][3];
 }
 Rcpp::RObject ld_swap = R_NilValue;
 if (updateLDswap && raw.size() > 26) {
  Rcpp::NumericMatrix ld(nt, 3);
  for (int t = 0; t < nt; ++t) {
   for (int j = 0; j < 3 && static_cast<std::size_t>(j) < raw[26][t].size(); ++j) ld(t, j) = raw[26][t][j];
  }
  ld_swap = ld;
 }
 Rcpp::IntegerVector group_index_r(m);
 Rcpp::IntegerVector group_size(ngroup);
 for (int i = 0; i < m; ++i) {
  group_index_r[i] = group_index[static_cast<std::size_t>(i)];
  if (group_index_r[i] >= 0 && group_index_r[i] < ngroup) ++group_size[group_index_r[i]];
 }
 Rcpp::RObject chains = keep_chains ? Rcpp::RObject(cpg_group_chains_raw_v1(raw, m, nt, ngroup, nchains)) : Rcpp::RObject(R_NilValue);
 Rcpp::List raw_out = Rcpp::List::create(
  Rcpp::_["schema"] = Rcpp::List::create(Rcpp::_["class"] = "stblr_raw", Rcpp::_["version"] = 1),
  Rcpp::_["meta"] = Rcpp::List::create(
   Rcpp::_["model"] = "bayesc",
   Rcpp::_["backend"] = "csr_group_bayesc",
   Rcpp::_["data_level"] = "summary",
   Rcpp::_["prior_type"] = "group",
   Rcpp::_["m"] = m,
   Rcpp::_["nt"] = nt,
   Rcpp::_["n_trace"] = n_trace,
   Rcpp::_["nit"] = nit,
   Rcpp::_["nburn"] = nburn,
   Rcpp::_["nthin"] = nthin,
   Rcpp::_["nchains"] = nchains,
   Rcpp::_["keep_chains"] = keep_chains,
   Rcpp::_["n_components"] = 2,
   Rcpp::_["n_annotations"] = 0,
   Rcpp::_["n_groups"] = ngroup
  ),
  Rcpp::_["marker"] = Rcpp::List::create(
   Rcpp::_["bm"] = cpg_raw_marker_matrix(raw, 0, m, nt),
   Rcpp::_["dm"] = cpg_raw_marker_matrix(raw, 1, m, nt),
   Rcpp::_["wy"] = cpg_raw_marker_matrix(raw, 2, m, nt),
   Rcpp::_["r"] = cpg_raw_marker_matrix(raw, 3, m, nt),
   Rcpp::_["b"] = cpg_raw_marker_matrix(raw, 4, m, nt),
   Rcpp::_["state"] = cpg_raw_marker_matrix(raw, 5, m, nt),
   Rcpp::_["bm_sd"] = raw.size() > 28 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 28, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["bm_min"] = raw.size() > 29 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 29, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["bm_max"] = raw.size() > 30 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 30, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["dm_sd"] = raw.size() > 31 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 31, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["dm_min"] = raw.size() > 32 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 32, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["dm_max"] = raw.size() > 33 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 33, m, nt)) : Rcpp::RObject(R_NilValue)
  ),
  Rcpp::_["trace"] = Rcpp::List::create(
   Rcpp::_["vbs"] = cpg_raw_trace_matrix(raw, 7, n_trace, nt),
   Rcpp::_["vgs"] = cpg_raw_trace_matrix(raw, 8, n_trace, nt),
   Rcpp::_["ves"] = cpg_raw_trace_matrix(raw, 9, n_trace, nt),
   Rcpp::_["vle"] = cpg_raw_trace_matrix(raw, 20, n_trace, nt),
   Rcpp::_["vld"] = cpg_raw_trace_matrix(raw, 21, n_trace, nt),
   Rcpp::_["pis"] = cpg_raw_trace_matrix(raw, 27, n_trace, nt)
  ),
  Rcpp::_["variance"] = Rcpp::List::create(
   Rcpp::_["covb"] = cpg_raw_trait_matrix(raw, 10, nt, nt),
   Rcpp::_["covg"] = cpg_raw_trait_matrix(raw, 11, nt, nt),
   Rcpp::_["cove"] = cpg_raw_trait_matrix(raw, 12, nt, nt),
   Rcpp::_["vb"] = cpg_raw_trait_matrix(raw, 13, nt, nt),
   Rcpp::_["vg"] = cpg_raw_trait_matrix(raw, 14, nt, nt),
   Rcpp::_["ve"] = cpg_raw_trait_matrix(raw, 15, nt, nt)
  ),
  Rcpp::_["pi"] = Rcpp::List::create(
   Rcpp::_["final"] = cpg_raw_trait_matrix(raw, 16, nt, 2),
   Rcpp::_["mean"] = cpg_raw_trait_matrix(raw, 17, nt, 2),
   Rcpp::_["names"] = Rcpp::CharacterVector::create("pi0", "pi1")
  ),
  Rcpp::_["diagnostics"] = Rcpp::List::create(
   Rcpp::_["nsamples"] = nsamples,
   Rcpp::_["n_used"] = n_used,
   Rcpp::_["log_cpo"] = log_cpo,
   Rcpp::_["mean_log_cpo"] = mean_log_cpo,
   Rcpp::_["seconds_mean"] = seconds_mean,
   Rcpp::_["seconds_max"] = seconds_max,
   Rcpp::_["ld_swap"] = ld_swap
  ),
  Rcpp::_["chains"] = chains,
  Rcpp::_["prior"] = Rcpp::List::create(),
  Rcpp::_["group"] = Rcpp::List::create(
   Rcpp::_["group_index"] = group_index_r,
   Rcpp::_["pi_mean"] = cpg_raw_group_matrix(raw, 22, ngroup, nt),
   Rcpp::_["pi_final"] = R_NilValue,
   Rcpp::_["vb_multiplier_mean"] = cpg_raw_group_matrix(raw, 23, ngroup, nt),
   Rcpp::_["vb_multiplier_final"] = R_NilValue,
   Rcpp::_["n_included_mean"] = cpg_raw_group_matrix(raw, 24, ngroup, nt),
   Rcpp::_["size"] = group_size
  ),
  Rcpp::_["annotation"] = Rcpp::List::create(),
  Rcpp::_["component"] = Rcpp::List::create(),
  Rcpp::_["selection"] = cpg_disabled_selection()
 );
 raw_out.attr("class") = Rcpp::CharacterVector::create("stblr_raw_v1", "stblr_raw", "list");
 return raw_out;
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_group_annot(
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
  std::vector<int> group_index,
  int ngroup,
  std::vector<std::vector<double>> group_pi_init,
  std::vector<double> pi_group_prior_a,
  std::vector<double> pi_group_prior_b,
  std::vector<std::vector<double>> group_vb_multiplier_init,
  bool updateGroupVb,
  double nub_group,
  double ssb_group_prior,
  bool normalize_group_vb,
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
  int ncores,
  int seed,
  int nchains = 1,
  bool keep_chains = false,
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds = R_NilValue,
  bool updateLDswap = false,
  double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8,
  int ld_swap_max_friends = 50,
  int ld_swap_moves = 1
) {
 if (nchains <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: nchains must be positive.");
 }
 std::vector<int> chain_seeds_vec;
 if (chain_seeds.isNotNull()) {
  Rcpp::IntegerVector chain_seeds_r(chain_seeds);
  chain_seeds_vec = Rcpp::as<std::vector<int>>(chain_seeds_r);
 }
 if (!chain_seeds_vec.empty() && static_cast<int>(chain_seeds_vec.size()) != nchains) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: chain_seeds must have length nchains.");
 }
 if (!std::isfinite(ld_swap_prob) || ld_swap_prob < 0.0 || ld_swap_prob > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ld_swap_prob must be in [0, 1].");
 }
 if (!std::isfinite(ld_swap_r2) || ld_swap_r2 < 0.0 || ld_swap_r2 > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ld_swap_r2 must be in [0, 1].");
 }
 if (ld_swap_max_friends <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ld_swap_max_friends must be positive.");
 }
 if (ld_swap_moves < 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ld_swap_moves must be non-negative.");
 }

 std::vector<std::vector<std::vector<double>>> out;
 std::vector<std::vector<std::vector<double>>> sumsq;
 std::vector<std::vector<std::vector<double>>> minv;
 std::vector<std::vector<std::vector<double>>> maxv;
 std::vector<std::vector<double>> chain_dm_flat;
 std::vector<std::vector<double>> chain_bm_flat;
 std::vector<std::vector<double>> chain_ld_swap_flat;
 std::vector<std::vector<double>> chain_group_pi_flat;
 std::vector<std::vector<double>> chain_group_vb_flat;
 std::vector<std::vector<double>> chain_group_nincluded_flat;
 std::vector<std::vector<double>> ld_swap_attempted_sum;
 std::vector<std::vector<double>> ld_swap_accepted_sum;

 for (int chain = 0; chain < nchains; ++chain) {
  int chain_seed = seed;
  if (!chain_seeds_vec.empty()) {
   chain_seed = chain_seeds_vec[static_cast<std::size_t>(chain)];
  } else if (nchains > 1) {
   chain_seed = seed + 9176 * (chain + 1);
  }

  std::vector<std::vector<std::vector<double>>> raw =
   stblr_cpg_omp_csr_group_annot_single(
    wy, ww, yy, b_init, d_init, use_d_init, r_init, use_r_init,
    rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, pi,
    group_index, ngroup, group_pi_init, pi_group_prior_a, pi_group_prior_b,
    group_vb_multiplier_init, updateGroupVb, nub_group, ssb_group_prior,
    normalize_group_vb, nub, nue, updateB, updateE, updatePi, adjE, n,
    nit, nburn, nthin, ncores, chain_seed, updateLDswap, ld_swap_prob,
    ld_swap_r2, ld_swap_max_friends, ld_swap_moves
   );

  if (chain == 0) {
   out = raw;
   sumsq = raw;
   minv = raw;
   maxv = raw;
   for (std::size_t k = 0; k < out.size(); ++k) {
    for (std::size_t t = 0; t < out[k].size(); ++t) {
     for (std::size_t i = 0; i < out[k][t].size(); ++i) {
      sumsq[k][t][i] = raw[k][t][i] * raw[k][t][i];
     }
    }
   }
   chain_dm_flat.resize(raw[1].size());
   chain_bm_flat.resize(raw[0].size());
   chain_ld_swap_flat.resize(raw[26].size());
   chain_group_pi_flat.resize(raw[22].size());
   chain_group_vb_flat.resize(raw[23].size());
   chain_group_nincluded_flat.resize(raw[24].size());
   ld_swap_attempted_sum.assign(raw[0].size(), std::vector<double>(1, 0.0));
   ld_swap_accepted_sum.assign(raw[0].size(), std::vector<double>(1, 0.0));
  } else {
   for (std::size_t k = 0; k < out.size(); ++k) {
    for (std::size_t t = 0; t < out[k].size(); ++t) {
     for (std::size_t i = 0; i < out[k][t].size(); ++i) {
      out[k][t][i] += raw[k][t][i];
      sumsq[k][t][i] += raw[k][t][i] * raw[k][t][i];
      minv[k][t][i] = std::min(minv[k][t][i], raw[k][t][i]);
      maxv[k][t][i] = std::max(maxv[k][t][i], raw[k][t][i]);
     }
    }
   }
  }

  for (std::size_t t = 0; t < raw[0].size(); ++t) {
   if (raw.size() > 26 && raw[26][t].size() >= 2) {
    ld_swap_attempted_sum[t][0] += raw[26][t][0];
    ld_swap_accepted_sum[t][0] += raw[26][t][1];
   }
  }

  if (keep_chains) {
   for (std::size_t t = 0; t < raw[0].size(); ++t) {
    chain_dm_flat[t].insert(chain_dm_flat[t].end(), raw[1][t].begin(), raw[1][t].end());
    chain_bm_flat[t].insert(chain_bm_flat[t].end(), raw[0][t].begin(), raw[0][t].end());
    if (raw.size() > 26 && raw[26][t].size() >= 4) {
     chain_ld_swap_flat[t].insert(
      chain_ld_swap_flat[t].end(),
      raw[26][t].begin(),
      raw[26][t].begin() + 4
     );
    } else {
     chain_ld_swap_flat[t].insert(chain_ld_swap_flat[t].end(), 4, 0.0);
    }
    chain_group_pi_flat[t].insert(chain_group_pi_flat[t].end(), raw[22][t].begin(), raw[22][t].end());
    chain_group_vb_flat[t].insert(chain_group_vb_flat[t].end(), raw[23][t].begin(), raw[23][t].end());
    chain_group_nincluded_flat[t].insert(
     chain_group_nincluded_flat[t].end(),
     raw[24][t].begin(),
     raw[24][t].end()
    );
   }
  }
 }

 const double inv_chains = 1.0 / static_cast<double>(nchains);
 for (std::size_t k = 0; k < out.size(); ++k) {
  for (std::size_t t = 0; t < out[k].size(); ++t) {
   for (std::size_t i = 0; i < out[k][t].size(); ++i) {
    out[k][t][i] *= inv_chains;
   }
  }
 }

 if (out.size() > 26) {
  for (std::size_t t = 0; t < out[26].size(); ++t) {
   if (out[26][t].size() < 4) out[26][t].resize(4);
   const double attempted = ld_swap_attempted_sum[t][0];
   const double accepted = ld_swap_accepted_sum[t][0];
   out[26][t][0] = attempted;
   out[26][t][1] = accepted;
   out[26][t][2] = attempted > 0.0 ? accepted / attempted : 0.0;
   out[26][t][3] = 1.0;
  }
 }

 const bool return_chain_summaries = (nchains > 1) || keep_chains;
 if (return_chain_summaries) {
  std::vector<std::vector<std::vector<double>>> extended(keep_chains ? 40 : 34);
  for (std::size_t k = 0; k < out.size(); ++k) extended[k] = out[k];
  for (int slot = 28; slot <= 33; ++slot) extended[static_cast<std::size_t>(slot)].resize(out[0].size());
  for (std::size_t t = 0; t < out[0].size(); ++t) {
   const std::size_t m_t = out[0][t].size();
   for (int slot = 28; slot <= 33; ++slot) extended[static_cast<std::size_t>(slot)][t].resize(m_t);
   for (std::size_t i = 0; i < m_t; ++i) {
    const double bm_mean = out[0][t][i];
    const double dm_mean = out[1][t][i];
    double bm_sd = 0.0;
    double dm_sd = 0.0;
    if (nchains > 1) {
     bm_sd = std::sqrt(std::max(0.0, (sumsq[0][t][i] - nchains * bm_mean * bm_mean) / static_cast<double>(nchains - 1)));
     dm_sd = std::sqrt(std::max(0.0, (sumsq[1][t][i] - nchains * dm_mean * dm_mean) / static_cast<double>(nchains - 1)));
    }
    extended[28][t][i] = bm_sd;
    extended[29][t][i] = minv[0][t][i];
    extended[30][t][i] = maxv[0][t][i];
    extended[31][t][i] = dm_sd;
    extended[32][t][i] = minv[1][t][i];
    extended[33][t][i] = maxv[1][t][i];
   }
  }
  if (keep_chains) {
   extended[34] = chain_dm_flat;
   extended[35] = chain_bm_flat;
   extended[36] = chain_ld_swap_flat;
   extended[37] = chain_group_pi_flat;
   extended[38] = chain_group_vb_flat;
   extended[39] = chain_group_nincluded_flat;
  }
  return cpg_group_raw_v1(
   extended, group_index, updateLDswap, static_cast<int>(wy[0].size()),
   static_cast<int>(wy.size()), ngroup, nit, nburn, nthin, nchains, keep_chains
  );
 }

 return cpg_group_raw_v1(
  out, group_index, updateLDswap, static_cast<int>(wy[0].size()),
  static_cast<int>(wy.size()), ngroup, nit, nburn, nthin, nchains, keep_chains
 );
}
