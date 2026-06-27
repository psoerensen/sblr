// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "st_csr_common.h"


#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

// -----------------------------------------------------------------------------
// ST-specific LD structure: flat symmetric CSR
// Stores pre-scaled X_i'X_j, not raw LD correlation.
// Disk input is expected to be upper-triangular or otherwise non-symmetric CSR
// with raw correlations r_ij. This builder symmetrizes it.
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Single-trait BayesC marker update
// -----------------------------------------------------------------------------

inline void sampleBetaC_ST_csr(
  int i,
  const std::vector<double>& pi,
  double vb,
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

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm01(0.0, 1.0);

 const double pi0 = std::max(pi[0], 1e-300);
 const double pi1 = std::max(pi[1], 1e-300);

 const double vei_safe = std::max(vei_i, 1e-300);

 // score = x_i' residual_without_i
 const double score = r(iu) + wi * b(iu);

 // Same BayesC scalar marginal likelihood as old sbayes(),
 // but using scalar adjusted residual variance vei_i.
 const double denom = std::max(vei_safe + wi * vb, 1e-300);

 const double logBF =
  0.5 * std::log(vei_safe / denom)
  + 0.5 * score * score * vb / (vei_safe * denom);

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
  const double lhs = wi + vei_safe / vb;
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

// -----------------------------------------------------------------------------
// Single-trait variance updates
// -----------------------------------------------------------------------------

inline void sampleB_ST_csr(
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

struct LDLDFriends {
 std::vector<uint64_t> ptr;
 std::vector<int> idx;
 std::vector<double> r2;
};

inline LDLDFriends build_ld_swap_friends_st_csr(
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

 LDLDFriends friends;
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

inline void set_marker_effect_st_csr(
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

inline bool attempt_ld_swap_st_csr(
  int m,
  double vei,
  double yy,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const STLDCSR& ld,
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

 set_marker_effect_st_csr(j, 0.0, 0, ww, r, b, d, ld);
 set_marker_effect_st_csr(k, b_j_old, 1, ww, r, b, d, ld);

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
    -0.5 * (sse_new - sse_old) / vei + log_q_reverse - log_q_forward;

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

// -----------------------------------------------------------------------------
// Main exported function: parallel single-trait BayesC over traits
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr(
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
  bool updateLDswap = false,
  double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8,
  int ld_swap_max_friends = 50,
  int ld_swap_moves = 1
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
 arma::Mat<int> d_mat(nt, m, arma::fill::zeros);

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

 STLDCSR ld = read_and_build_st_ld_csr(
  ld_prefix,
  m,
  xx
 );

 LDLDFriends ld_swap_friends;
 if (updateLDswap) {
  ld_swap_friends = build_ld_swap_friends_st_csr(
   m,
   ld,
   xx,
   ld_swap_r2,
   ld_swap_max_friends
  );
 } else {
  ld_swap_friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
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

 // --------------------------------------------------------------------------
 // Output storage
 // --------------------------------------------------------------------------

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);

 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vles_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vlds_mat(nt, nit + nburn, arma::fill::zeros);

 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi(nt, arma::fill::zeros);
 arma::vec final_vle(nt, arma::fill::zeros);
 arma::vec final_vld(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);
 arma::vec ld_swap_attempted_vec(nt, arma::fill::zeros);
 arma::vec ld_swap_accepted_vec(nt, arma::fill::zeros);

 // --------------------------------------------------------------------------
 // Parallel over traits
 // --------------------------------------------------------------------------

 std::vector<int> failed(static_cast<std::size_t>(nt), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(nt));
 std::vector<int> thread_used(static_cast<std::size_t>(nt), 0);
 std::vector<double> trait_seconds(static_cast<std::size_t>(nt), 0.0);

 int nthreads = 1;

#ifdef _OPENMP
 omp_set_dynamic(0);

 nthreads = std::max(1, std::min(ncores, nt));

 // Explicitly set the OpenMP thread count for this runtime.
 omp_set_num_threads(nthreads);

 Rcpp::Rcout
 << "STBLR OpenMP requested threads = "
 << nthreads
 << ", omp_get_max_threads = "
 << omp_get_max_threads()
 << ", num procs = "
 << omp_get_num_procs()
 << ", pi_prior_a = "
 << pi_prior_a
 << ", pi_prior_b = "
 << pi_prior_b
 << "\n";
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
 for (int t = 0; t < nt; ++t) {

#ifdef _OPENMP
  const double wall_start = omp_get_wtime();
  thread_used[static_cast<std::size_t>(t)] = omp_get_thread_num();
#else
  const double wall_start = 0.0;
  thread_used[static_cast<std::size_t>(t)] = 0;
#endif

  try {
   std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));

   arma::rowvec b_t(m, arma::fill::zeros);
   for (int i = 0; i < m; ++i) {
    b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   }

   arma::rowvec r_t(m, arma::fill::zeros);
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
    rebuild_residual_st_csr(
     m,
     wy_t,
     ww_t,
     b_t,
     r_t,
     ld
    );
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

   double nsamples_t = 0.0;
   double ld_swap_attempted_t = 0.0;
   double ld_swap_accepted_t = 0.0;

   for (int it = 0; it < nit + nburn; ++it) {

    // -------------------------------------------------------
    // Marker updates
    // -------------------------------------------------------
    for (int isort = 0; isort < m; ++isort) {
     const int i = order[static_cast<std::size_t>(isort)];

     sampleBetaC_ST_csr(
      i,
      pi_t,
      vb_t,
      vei_t,
      ww_t,
      r_t,
      b_t,
      d_t,
      ld,
      gen_t
     );
    }

    if (updateLDswap && ld_swap_moves > 0 && ld_swap_prob > 0.0) {
     std::uniform_real_distribution<double> runif(0.0, 1.0);

     if (runif(gen_t) < ld_swap_prob) {
      for (int move = 0; move < ld_swap_moves; ++move) {
       ld_swap_attempted_t += 1.0;
       if (attempt_ld_swap_st_csr(
            m,
            vei_t,
            yy_vec(static_cast<arma::uword>(t)),
            ww_t,
            wy_t,
            r_t,
            b_t,
            d_t,
            ld,
            ld_swap_friends,
            gen_t
           )) {
        ld_swap_accepted_t += 1.0;
       }
      }
     }
    }

    // -------------------------------------------------------
    // Variance updates
    // -------------------------------------------------------
    if (updateB) {
     sampleB_ST_csr(
      m,
      nub,
      vb_t,
      b_t,
      d_t,
      ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      gen_t
     );

     if (!std::isfinite(vb_t) || vb_t <= 0.0) {
      throw std::runtime_error(
        "vb became invalid after sampleB. iter=" +
         std::to_string(it) +
         ", vb=" + std::to_string(vb_t)
      );
     }
    }

    if (updateE) {
     if (rebuild_r_before_updateE) {
      rebuild_residual_st_csr(
       m,
       wy_t,
       ww_t,
       b_t,
       r_t,
       ld
      );
     }

     sampleE_ST_csr(
      m,
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

    vg_t = computeG_ST_csr(
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

    // -------------------------------------------------------
    // Store posterior summaries
    // -------------------------------------------------------
    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;

     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += static_cast<double>(d_t(iu));
     }
    }
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

   bm_mat.row(static_cast<arma::uword>(t)) = bm_t;
   dm_mat.row(static_cast<arma::uword>(t)) = dm_t;
   b_mat.row(static_cast<arma::uword>(t))  = b_t;
   r_mat.row(static_cast<arma::uword>(t))  = r_t;
   d_mat.row(static_cast<arma::uword>(t))  = d_t;

   vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
   vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
   ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
   pis_mat.row(static_cast<arma::uword>(t)) = pis_t;
   vles_mat.row(static_cast<arma::uword>(t)) = vles_t;
   vlds_mat.row(static_cast<arma::uword>(t)) = vlds_t;

   final_vb(static_cast<arma::uword>(t)) = vb_t;
   final_ve(static_cast<arma::uword>(t)) = ve_t;
   final_vg(static_cast<arma::uword>(t)) = vg_t;
   final_vle(static_cast<arma::uword>(t)) = vle_t;
   final_vld(static_cast<arma::uword>(t)) = vld_t;
   final_pi(static_cast<arma::uword>(t)) = pi_t[1];
   nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;
   ld_swap_attempted_vec(static_cast<arma::uword>(t)) = ld_swap_attempted_t;
   ld_swap_accepted_vec(static_cast<arma::uword>(t)) = ld_swap_accepted_t;

#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif

  } catch (const std::exception& e) {
   failed[static_cast<std::size_t>(t)] = 1;
   errors[static_cast<std::size_t>(t)] = e.what();
#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif
  } catch (...) {
   failed[static_cast<std::size_t>(t)] = 1;
   errors[static_cast<std::size_t>(t)] = "unknown error";
#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif
  }
 }

#ifdef _OPENMP
 for (int t = 0; t < nt; ++t) {
  Rcpp::Rcout
  << "trait " << t
  << " used thread " << thread_used[static_cast<std::size_t>(t)]
  << ", seconds = " << trait_seconds[static_cast<std::size_t>(t)]
  << "\n";
 }
#endif

 for (int t = 0; t < nt; ++t) {
  if (failed[static_cast<std::size_t>(t)]) {
   throw std::runtime_error(
     "stblr_cpg_omp_csr failed for trait " +
      std::to_string(t) +
      ": " +
      errors[static_cast<std::size_t>(t)]
   );
  }
 }

 // --------------------------------------------------------------------------
 // Build result with same style as MT output
 // --------------------------------------------------------------------------

 std::vector<std::vector<std::vector<double>>> result(23);

 for (int k = 0; k < 23; ++k) {
  result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
 }

 for (int t = 0; t < nt; ++t) {
  result[0][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // bm
  result[1][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // dm
  result[2][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // wy
  result[3][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // r
  result[4][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // b
  result[5][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // d
  result[6][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // marker index

  result[7][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vbs
  result[8][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vgs
  result[9][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // ves

  result[10][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // covb
  result[11][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // covg
  result[12][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // cove
  result[13][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final B
  result[14][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final G
  result[15][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final E

  result[16][static_cast<std::size_t>(t)].resize(2);                                     // final pi
  result[17][static_cast<std::size_t>(t)].resize(2);                                     // posterior mean pi approx

  result[18][static_cast<std::size_t>(t)].resize(4);
  result[19][static_cast<std::size_t>(t)].resize(2);
  result[20][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vle
  result[21][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vld = vg - vle
  result[22][static_cast<std::size_t>(t)].resize(4);                                       // LD-swap diagnostics
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int i = 0; i < m; ++i) {
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword iu = static_cast<arma::uword>(i);
   const std::size_t is = static_cast<std::size_t>(i);

   result[0][ts][is] = bm_mat(tu, iu);
   result[1][ts][is] = dm_mat(tu, iu);
   result[2][ts][is] = wy_mat(tu, iu);
   result[3][ts][is] = r_mat(tu, iu);
   result[4][ts][is] = b_mat(tu, iu);
   result[5][ts][is] = static_cast<double>(d_mat(tu, iu));
   result[6][ts][is] = static_cast<double>(i);
  }
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int it = 0; it < nit + nburn; ++it) {
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword itu = static_cast<arma::uword>(it);
   const std::size_t its = static_cast<std::size_t>(it);

   result[7][ts][its] = vbs_mat(tu, itu);
   result[8][ts][its] = vgs_mat(tu, itu);
   result[9][ts][its] = ves_mat(tu, itu);
   result[20][ts][its] = vles_mat(tu, itu);
   result[21][ts][its] = vlds_mat(tu, itu);
  }
 }

 for (int t1 = 0; t1 < nt; ++t1) {
  const std::size_t t1s = static_cast<std::size_t>(t1);

  for (int t2 = 0; t2 < nt; ++t2) {
   const std::size_t t2s = static_cast<std::size_t>(t2);

   result[10][t1s][t2s] = 0.0;
   result[11][t1s][t2s] = 0.0;
   result[12][t1s][t2s] = 0.0;

   result[13][t1s][t2s] = 0.0;
   result[14][t1s][t2s] = 0.0;
   result[15][t1s][t2s] = 0.0;
  }

  result[10][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
  result[11][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
  result[12][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));

  result[13][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
  result[14][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
  result[15][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  result[16][ts][0] = 1.0 - final_pi(static_cast<arma::uword>(t));
  result[16][ts][1] = final_pi(static_cast<arma::uword>(t));

  // Approximate posterior mean inclusion pi from saved trace.
  double mean_pi = 0.0;
  int npi = 0;

  for (int it = nburn; it < nit + nburn; ++it) {
   mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
   ++npi;
  }

  if (npi > 0) mean_pi /= static_cast<double>(npi);
  else mean_pi = final_pi(static_cast<arma::uword>(t));

  result[17][ts][0] = 1.0 - mean_pi;
  result[17][ts][1] = mean_pi;

  for (int i = 0; i < 4; ++i) {
   result[18][ts][static_cast<std::size_t>(i)] = 0.0;
  }

  for (int i = 0; i < 2; ++i) {
   result[19][ts][static_cast<std::size_t>(i)] = 0.0;
  }

  result[22][ts][0] = ld_swap_attempted_vec(static_cast<arma::uword>(t));
  result[22][ts][1] = ld_swap_accepted_vec(static_cast<arma::uword>(t));
  result[22][ts][2] =
   (ld_swap_attempted_vec(static_cast<arma::uword>(t)) > 0.0)
   ? ld_swap_accepted_vec(static_cast<arma::uword>(t)) /
     ld_swap_attempted_vec(static_cast<arma::uword>(t))
   : 0.0;
  result[22][ts][3] = 1.0;
 }

 return result;
}


// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
//
// #include "cpg_samplers.h"
// #include "distributions.h"
// #include "st_csr_common.h"
//
//
// #include <algorithm>
// #include <cmath>
// #include <cstdio>
// #include <cstdint>
// #include <fstream>
// #include <limits>
// #include <numeric>
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
// #include <unordered_map>
// #include <unordered_set>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // -----------------------------------------------------------------------------
// // ST-specific LD structure: flat symmetric CSR
// // Stores pre-scaled X_i'X_j, not raw LD correlation.
// // Disk input is expected to be upper-triangular or otherwise non-symmetric CSR
// // with raw correlations r_ij. This builder symmetrizes it.
// // -----------------------------------------------------------------------------
//
//
// // -----------------------------------------------------------------------------
// // Single-trait BayesC marker update
// // -----------------------------------------------------------------------------
//
// inline void sampleBetaC_ST_csr(
//   int i,
//   const std::vector<double>& pi,
//   double vb,
//   double vei_i,
//   const arma::rowvec& ww,
//   arma::rowvec& r,
//   arma::rowvec& b,
//   arma::Row<int>& d,
//   const STLDCSR& ld,
//   std::mt19937& gen
// ) {
//  const arma::uword iu = static_cast<arma::uword>(i);
//
//  const double wi = ww(iu);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm01(0.0, 1.0);
//
//  const double pi0 = std::max(pi[0], 1e-300);
//  const double pi1 = std::max(pi[1], 1e-300);
//
//  const double vei_safe = std::max(vei_i, 1e-300);
//
//  // score = x_i' residual_without_i
//  const double score = r(iu) + wi * b(iu);
//
//  // Same BayesC scalar marginal likelihood as old sbayes(),
//  // but using scalar adjusted residual variance vei_i.
//  const double denom = std::max(vei_safe + wi * vb, 1e-300);
//
//  const double logBF =
//   0.5 * std::log(vei_safe / denom)
//   + 0.5 * score * score * vb / (vei_safe * denom);
//
//  const double logp1 = std::log(pi1) + logBF;
//  const double logp0 = std::log(pi0);
//
//  const double delta_log = logp0 - logp1;
//
//  double p1 = 0.0;
//
//  if (delta_log > 35.0) {
//   p1 = 0.0;
//  } else if (delta_log < -35.0) {
//   p1 = 1.0;
//  } else {
//   p1 = 1.0 / (1.0 + std::exp(delta_log));
//  }
//
//  const int di = (runif(gen) < p1) ? 1 : 0;
//
//  double b_new = 0.0;
//
//  if (di == 1) {
//   const double lhs = wi + vei_safe / vb;
//   const double mean = score / lhs;
//   const double sd = std::sqrt(vei_safe / lhs);
//
//   b_new = mean + sd * norm01(gen);
//  }
//
//  const double diff = b_new - b(iu);
//
//  if (diff != 0.0) {
//   r(iu) -= wi * diff;
//
//   const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
//   const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];
//
//   for (uint64_t p = start; p < end; ++p) {
//    const int j = ld.idx[static_cast<std::size_t>(p)];
//    r(static_cast<arma::uword>(j)) -=
//     static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
//   }
//  }
//
//  b(iu) = b_new;
//  d(iu) = di;
// }
//
// // -----------------------------------------------------------------------------
// // Single-trait variance updates
// // -----------------------------------------------------------------------------
//
// inline void sampleB_ST_csr(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   double ssb_prior,
//   std::mt19937& gen
// ) {
//  double ssb = 0.0;
//  double dfb = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//   if (d(iu) > 0) {
//    ssb += b(iu) * b(iu);
//    dfb += 1.0;
//   }
//  }
//
//  const double scale = ssb + nub * ssb_prior;
//
//  std::chi_squared_distribution<double> rchisq(dfb + nub);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  vb = std::max(scale / chi2, 1e-12);
// }
//
//
// inline void samplePi_ST(
//   const arma::Row<int>& d,
//   std::vector<double>& pi,
//   std::mt19937& gen
// ) {
//  double c0 = 1.0;
//  double c1 = 1.0;
//
//  for (arma::uword i = 0; i < d.n_elem; ++i) {
//   if (d(i) > 0) c1 += 1.0;
//   else c0 += 1.0;
//  }
//
//  std::gamma_distribution<double> rg0(c0, 1.0);
//  std::gamma_distribution<double> rg1(c1, 1.0);
//
//  const double g0 = std::max(rg0(gen), 1e-300);
//  const double g1 = std::max(rg1(gen), 1e-300);
//  const double s = g0 + g1;
//
//  pi[0] = g0 / s;
//  pi[1] = g1 / s;
// }
//
// // -----------------------------------------------------------------------------
// // Linkage-equilibrium and linkage-disequilibrium variance components
// // -----------------------------------------------------------------------------
//
// inline double computeLE_ST_csr(
//   int m,
//   const arma::rowvec& b,
//   const arma::rowvec& ww,
//   int n
// ) {
//  double vle = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//   const double bi = b(iu);
//
//   if (bi != 0.0) {
//    vle += ww(iu) * bi * bi;
//   }
//  }
//
//  return vle / static_cast<double>(n);
// }
//
// // -----------------------------------------------------------------------------
// // Main exported function: parallel single-trait BayesC over traits
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr(
//   std::vector<std::vector<double>> wy,
//   std::vector<std::vector<double>> ww,
//   std::vector<double> yy,
//   std::vector<std::vector<double>> b_init,
//   std::vector<std::vector<double>> d_init,
//   bool use_d_init,
//   std::vector<std::vector<double>> r_init,
//   bool use_r_init,
//   bool rebuild_r_before_updateE,
//   std::string ld_prefix,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   std::vector<int> n,
//   int nit,
//   int nburn,
//   int nthin,
//   int ncores,
//   int seed
// ) {
//  const int nt = static_cast<int>(wy.size());
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nt must be positive.");
//  }
//
//  const int m = static_cast<int>(wy[0].size());
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: m must be positive.");
//  }
//
//  if (nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nit must be positive.");
//  }
//
//  if (nburn < 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nburn must be non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nthin must be positive.");
//  }
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt ||
//      (int)yy.size() != nt || (int)n.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: inconsistent trait dimensions.");
//  }
//
//  if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: priors must be nt x nt.");
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_csr: pi must have length 2, c(pi0, pi1).");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m ||
//       (int)ww[t].size() != m ||
//       (int)b_init[t].size() != m) {
//    throw std::runtime_error("stblr_cpg_omp_csr: inconsistent marker dimensions.");
//   }
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: E must be nt x nt.");
//  }
//
//  if (use_r_init) {
//   if (static_cast<int>(r_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_state: r_init must have length nt when use_r_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(r_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_state: each r_init[t] must have length m.");
//    }
//   }
//  }
//
//  if (use_d_init) {
//   if (static_cast<int>(d_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_state: d_init must have length nt when use_d_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(d_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_state: each d_init[t] must have length m.");
//    }
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Convert inputs to Armadillo
//  // --------------------------------------------------------------------------
//
//  arma::mat wy_mat(nt, m, arma::fill::zeros);
//  arma::mat ww_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::mat r_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
//
//  arma::vec yy_vec(nt, arma::fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   yy_vec(static_cast<arma::uword>(t)) = yy[t];
//
//   for (int i = 0; i < m; ++i) {
//    wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
//    ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
//    b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i))  = b_init[t][i];
//   }
//
//   if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
//    sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Validate shared scaling and build shared flat LD object
//  // --------------------------------------------------------------------------
//
//  for (int t = 1; t < nt; ++t) {
//   if (n[t] != n[0]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_csr: current shared-LD scaling assumes equal n across traits."
//    );
//   }
//  }
//
//  for (int t = 1; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    const double w0 = ww_mat(0, static_cast<arma::uword>(i));
//    const double wt = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    const double tol = 1e-8 * std::max(1.0, std::abs(w0));
//
//    if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr: ww contains invalid value before LD pre-scaling."
//     );
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr: ww differs across traits; pre-scaled shared ST LD is invalid."
//     );
//    }
//   }
//  }
//
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//
//  for (int i = 0; i < m; ++i) {
//   const double wi = ww_mat(0, static_cast<arma::uword>(i));
//   if (!std::isfinite(wi) || wi <= 0.0) {
//    throw std::runtime_error("stblr_cpg_omp_csr: ww contains invalid value in trait 0.");
//   }
//   xx[static_cast<std::size_t>(i)] = wi;
//  }
//
//  STLDCSR ld = read_and_build_st_ld_csr(
//   ld_prefix,
//   m,
//   xx
//  );
//
//  // --------------------------------------------------------------------------
//  // Marker update order based on max single-trait marginal effect
//  // --------------------------------------------------------------------------
//
//  std::vector<double> x2(static_cast<std::size_t>(m), 0.0);
//  std::vector<int> order(static_cast<std::size_t>(m));
//
//  for (int i = 0; i < m; ++i) {
//   double best = 0.0;
//
//   for (int t = 0; t < nt; ++t) {
//    const double wi = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    if (wi > 0.0) {
//     const double bhat = wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) / wi;
//     best = std::max(best, bhat * bhat);
//    }
//   }
//
//   x2[static_cast<std::size_t>(i)] = best;
//   order[static_cast<std::size_t>(i)] = i;
//  }
//
//  std::sort(order.begin(), order.end(),
//            [&](int a, int b) { return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)]; });
//
//  // --------------------------------------------------------------------------
//  // Output storage
//  // --------------------------------------------------------------------------
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//
//  arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vles_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vlds_mat(nt, nit + nburn, arma::fill::zeros);
//
//  arma::vec final_vb(nt, arma::fill::zeros);
//  arma::vec final_vg(nt, arma::fill::zeros);
//  arma::vec final_ve(nt, arma::fill::zeros);
//  arma::vec final_pi(nt, arma::fill::zeros);
//  arma::vec final_vle(nt, arma::fill::zeros);
//  arma::vec final_vld(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  // --------------------------------------------------------------------------
//  // Parallel over traits
//  // --------------------------------------------------------------------------
//
//  std::vector<int> failed(static_cast<std::size_t>(nt), 0);
//  std::vector<std::string> errors(static_cast<std::size_t>(nt));
//  std::vector<int> thread_used(static_cast<std::size_t>(nt), 0);
//  std::vector<double> trait_seconds(static_cast<std::size_t>(nt), 0.0);
//
//  int nthreads = 1;
//
// #ifdef _OPENMP
//  omp_set_dynamic(0);
//
//  nthreads = std::max(1, std::min(ncores, nt));
//
//  // Explicitly set the OpenMP thread count for this runtime.
//  omp_set_num_threads(nthreads);
//
//  Rcpp::Rcout
//  << "STBLR OpenMP requested threads = "
//  << nthreads
//  << ", omp_get_max_threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
//  // #ifdef _OPENMP
//  //  omp_set_dynamic(0);
//  //  nthreads = std::max(1, std::min(ncores, nt));
//  //  nthreads = std::min(nthreads, omp_get_max_threads());
//  //
//  //  Rcpp::Rcout
//  //  << "STBLR OpenMP threads = "
//  //  << nthreads
//  //  << ", max threads = "
//  //  << omp_get_max_threads()
//  //  << ", num procs = "
//  //  << omp_get_num_procs()
//  //  << "\n";
//  // #endif
//
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(nthreads) schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
//
// #ifdef _OPENMP
//   const double wall_start = omp_get_wtime();
//   thread_used[static_cast<std::size_t>(t)] = omp_get_thread_num();
// #else
//   const double wall_start = 0.0;
//   thread_used[static_cast<std::size_t>(t)] = 0;
// #endif
//
//   try {
//    std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//    arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
//    arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));
//
//    arma::rowvec b_t(m, arma::fill::zeros);
//    for (int i = 0; i < m; ++i) {
//     b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    }
//
//    arma::rowvec r_t(m, arma::fill::zeros);
//    arma::Row<int> d_t(m, arma::fill::zeros);
//
//    if (use_d_init) {
//     for (int i = 0; i < m; ++i) {
//      d_t(static_cast<arma::uword>(i)) = d_init[t][i] > 0 ? 1 : 0;
//     }
//    } else {
//     for (int i = 0; i < m; ++i) {
//      d_t(static_cast<arma::uword>(i)) =
//       (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
//     }
//    }
//
//    if (use_r_init) {
//     for (int i = 0; i < m; ++i) {
//      r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
//     }
//
//     if (!r_t.is_finite()) {
//      throw std::runtime_error("stblr_cpg_omp_csr_state: r_init contains NaN/Inf.");
//     }
//    } else {
//     rebuild_residual_st_csr(
//      m,
//      wy_t,
//      ww_t,
//      b_t,
//      r_t,
//      ld
//     );
//    }
//
//    double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double vg_t = 0.0;
//    double vle_t = computeLE_ST_csr(m, b_t, ww_t, n[t]);
//    double vld_t = vg_t - vle_t;
//    double vei_t = ve_t + adjE * vg_t;
//
//    std::vector<double> pi_t = pi;
//
//    if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//        pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//     throw std::runtime_error(
//       "invalid initial pi: pi0=" + std::to_string(pi_t[0]) +
//        ", pi1=" + std::to_string(pi_t[1])
//     );
//    }
//
//    {
//     const double psum = pi_t[0] + pi_t[1];
//
//     if (!std::isfinite(psum) || psum <= 0.0) {
//      throw std::runtime_error("invalid initial pi sum.");
//     }
//
//     pi_t[0] /= psum;
//     pi_t[1] /= psum;
//    }
//
//    arma::rowvec bm_t(m, arma::fill::zeros);
//    arma::rowvec dm_t(m, arma::fill::zeros);
//    arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//
//     // -------------------------------------------------------
//     // Marker updates
//     // -------------------------------------------------------
//     for (int isort = 0; isort < m; ++isort) {
//      const int i = order[static_cast<std::size_t>(isort)];
//
//      sampleBetaC_ST_csr(
//       i,
//       pi_t,
//       vb_t,
//       vei_t,
//       ww_t,
//       r_t,
//       b_t,
//       d_t,
//       ld,
//       gen_t
//      );
//     }
//
//     // -------------------------------------------------------
//     // Variance updates
//     // -------------------------------------------------------
//     if (updateB) {
//      sampleB_ST_csr(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
//       gen_t
//      );
//
//      if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//       throw std::runtime_error(
//         "vb became invalid after sampleB. iter=" +
//          std::to_string(it) +
//          ", vb=" + std::to_string(vb_t)
//       );
//      }
//     }
//
//     if (updateE) {
//      if (rebuild_r_before_updateE) {
//       rebuild_residual_st_csr(
//        m,
//        wy_t,
//        ww_t,
//        b_t,
//        r_t,
//        ld
//       );
//      }
//
//      sampleE_ST_csr(
//       m,
//       nue,
//       ve_t,
//       b_t,
//       wy_t,
//       r_t,
//       sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
//       yy_vec(static_cast<arma::uword>(t)),
//       n[t],
//        gen_t
//      );
//     }
//
//     if (updatePi) {
//      samplePi_ST(d_t, pi_t, gen_t);
//
//      if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//          pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//       throw std::runtime_error(
//         "pi became invalid after samplePi. iter=" +
//          std::to_string(it) +
//          ", pi0=" + std::to_string(pi_t[0]) +
//          ", pi1=" + std::to_string(pi_t[1])
//       );
//      }
//     }
//
//     vg_t = computeG_ST_csr(
//      b_t,
//      wy_t,
//      r_t,
//      n[t]
//     );
//
//     vle_t = computeLE_ST_csr(
//      m,
//      b_t,
//      ww_t,
//      n[t]
//     );
//
//     vld_t = vg_t - vle_t;
//
//     if (!std::isfinite(vg_t)) {
//      throw std::runtime_error(
//        "vg became NaN/Inf after computeG. iter=" +
//         std::to_string(it)
//      );
//     }
//
//     if (!std::isfinite(vle_t)) {
//      throw std::runtime_error(
//        "vle became NaN/Inf after computeLE. iter=" +
//         std::to_string(it)
//      );
//     }
//
//     if (!std::isfinite(vld_t)) {
//      throw std::runtime_error(
//        "vld became NaN/Inf after computeLE. iter=" +
//         std::to_string(it)
//      );
//     }
//
//     vei_t = ve_t + adjE * vg_t;
//
//     if (!std::isfinite(vei_t) || vei_t <= 0.0) {
//      throw std::runtime_error(
//        "adjusted residual variance vei became invalid. iter=" +
//         std::to_string(it) +
//         ", vei=" + std::to_string(vei_t)
//      );
//     }
//
//     vbs_t(static_cast<arma::uword>(it)) = vb_t;
//     ves_t(static_cast<arma::uword>(it)) = ve_t;
//     vgs_t(static_cast<arma::uword>(it)) = vg_t;
//     pis_t(static_cast<arma::uword>(it)) = pi_t[1];
//     vles_t(static_cast<arma::uword>(it)) = vle_t;
//     vlds_t(static_cast<arma::uword>(it)) = vld_t;
//
//     // -------------------------------------------------------
//     // Store posterior summaries
//     // -------------------------------------------------------
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       const arma::uword iu = static_cast<arma::uword>(i);
//       bm_t(iu) += b_t(iu);
//       dm_t(iu) += static_cast<double>(d_t(iu));
//      }
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    bm_t /= nsamples_t;
//    dm_t /= nsamples_t;
//
//    if (!bm_t.is_finite()) {
//     throw std::runtime_error("posterior mean bm contains NaN/Inf.");
//    }
//
//    if (!dm_t.is_finite()) {
//     throw std::runtime_error("posterior mean dm contains NaN/Inf.");
//    }
//
//    bm_mat.row(static_cast<arma::uword>(t)) = bm_t;
//    dm_mat.row(static_cast<arma::uword>(t)) = dm_t;
//    b_mat.row(static_cast<arma::uword>(t))  = b_t;
//    r_mat.row(static_cast<arma::uword>(t))  = r_t;
//    d_mat.row(static_cast<arma::uword>(t))  = d_t;
//
//    vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
//    vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
//    ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
//    pis_mat.row(static_cast<arma::uword>(t)) = pis_t;
//    vles_mat.row(static_cast<arma::uword>(t)) = vles_t;
//    vlds_mat.row(static_cast<arma::uword>(t)) = vlds_t;
//
//    final_vb(static_cast<arma::uword>(t)) = vb_t;
//    final_ve(static_cast<arma::uword>(t)) = ve_t;
//    final_vg(static_cast<arma::uword>(t)) = vg_t;
//    final_vle(static_cast<arma::uword>(t)) = vle_t;
//    final_vld(static_cast<arma::uword>(t)) = vld_t;
//    final_pi(static_cast<arma::uword>(t)) = pi_t[1];
//    nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;
//
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//
//   } catch (const std::exception& e) {
//    failed[static_cast<std::size_t>(t)] = 1;
//    errors[static_cast<std::size_t>(t)] = e.what();
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//   } catch (...) {
//    failed[static_cast<std::size_t>(t)] = 1;
//    errors[static_cast<std::size_t>(t)] = "unknown error";
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//   }
//  }
//
// #ifdef _OPENMP
//  for (int t = 0; t < nt; ++t) {
//   Rcpp::Rcout
//   << "trait " << t
//   << " used thread " << thread_used[static_cast<std::size_t>(t)]
//   << ", seconds = " << trait_seconds[static_cast<std::size_t>(t)]
//   << "\n";
//  }
// #endif
//
//  for (int t = 0; t < nt; ++t) {
//   if (failed[static_cast<std::size_t>(t)]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_csr failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[static_cast<std::size_t>(t)]
//    );
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Build result with same style as MT output
//  // --------------------------------------------------------------------------
//
//  std::vector<std::vector<std::vector<double>>> result(22);
//
//  for (int k = 0; k < 22; ++k) {
//   result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // bm
//   result[1][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // dm
//   result[2][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // wy
//   result[3][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // r
//   result[4][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // b
//   result[5][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // d
//   result[6][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // marker index
//
//   result[7][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vbs
//   result[8][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vgs
//   result[9][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // ves
//
//   result[10][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // covb
//   result[11][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // covg
//   result[12][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // cove
//   result[13][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final B
//   result[14][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final G
//   result[15][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final E
//
//   result[16][static_cast<std::size_t>(t)].resize(2);                                     // final pi
//   result[17][static_cast<std::size_t>(t)].resize(2);                                     // posterior mean pi approx
//
//   result[18][static_cast<std::size_t>(t)].resize(4);
//   result[19][static_cast<std::size_t>(t)].resize(2);
//   result[20][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vle
//   result[21][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vld = vg - vle
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int i = 0; i < m; ++i) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword iu = static_cast<arma::uword>(i);
//    const std::size_t is = static_cast<std::size_t>(i);
//
//    result[0][ts][is] = bm_mat(tu, iu);
//    result[1][ts][is] = dm_mat(tu, iu);
//    result[2][ts][is] = wy_mat(tu, iu);
//    result[3][ts][is] = r_mat(tu, iu);
//    result[4][ts][is] = b_mat(tu, iu);
//    result[5][ts][is] = static_cast<double>(d_mat(tu, iu));
//    result[6][ts][is] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int it = 0; it < nit + nburn; ++it) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword itu = static_cast<arma::uword>(it);
//    const std::size_t its = static_cast<std::size_t>(it);
//
//    result[7][ts][its] = vbs_mat(tu, itu);
//    result[8][ts][its] = vgs_mat(tu, itu);
//    result[9][ts][its] = ves_mat(tu, itu);
//    result[20][ts][its] = vles_mat(tu, itu);
//    result[21][ts][its] = vlds_mat(tu, itu);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   const std::size_t t1s = static_cast<std::size_t>(t1);
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    const std::size_t t2s = static_cast<std::size_t>(t2);
//
//    result[10][t1s][t2s] = 0.0;
//    result[11][t1s][t2s] = 0.0;
//    result[12][t1s][t2s] = 0.0;
//
//    result[13][t1s][t2s] = 0.0;
//    result[14][t1s][t2s] = 0.0;
//    result[15][t1s][t2s] = 0.0;
//   }
//
//   result[10][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
//   result[11][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
//   result[12][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
//
//   result[13][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
//   result[14][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
//   result[15][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   result[16][ts][0] = 1.0 - final_pi(static_cast<arma::uword>(t));
//   result[16][ts][1] = final_pi(static_cast<arma::uword>(t));
//
//   // Approximate posterior mean inclusion pi from saved trace.
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi(static_cast<arma::uword>(t));
//
//   result[17][ts][0] = 1.0 - mean_pi;
//   result[17][ts][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) {
//    result[18][ts][static_cast<std::size_t>(i)] = 0.0;
//   }
//
//   for (int i = 0; i < 2; ++i) {
//    result[19][ts][static_cast<std::size_t>(i)] = 0.0;
//   }
//  }
//
//  return result;
// }
//




//
// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
//
// #include "cpg_samplers.h"
// #include "distributions.h"
// #include "st_csr_common.h"
//
//
// #include <algorithm>
// #include <cmath>
// #include <cstdio>
// #include <cstdint>
// #include <fstream>
// #include <limits>
// #include <numeric>
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
// #include <unordered_map>
// #include <unordered_set>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // -----------------------------------------------------------------------------
// // ST-specific LD structure: flat symmetric CSR
// // Stores pre-scaled X_i'X_j, not raw LD correlation.
// // Disk input is expected to be upper-triangular or otherwise non-symmetric CSR
// // with raw correlations r_ij. This builder symmetrizes it.
// // -----------------------------------------------------------------------------
//
//
// // -----------------------------------------------------------------------------
// // Single-trait BayesC marker update
// // -----------------------------------------------------------------------------
//
// inline void sampleBetaC_ST_csr(
//   int i,
//   const std::vector<double>& pi,
//   double vb,
//   double vei_i,
//   const arma::rowvec& ww,
//   arma::rowvec& r,
//   arma::rowvec& b,
//   arma::Row<int>& d,
//   const STLDCSR& ld,
//   std::mt19937& gen
// ) {
//  const arma::uword iu = static_cast<arma::uword>(i);
//
//  const double wi = ww(iu);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm01(0.0, 1.0);
//
//  const double pi0 = std::max(pi[0], 1e-300);
//  const double pi1 = std::max(pi[1], 1e-300);
//
//  const double vei_safe = std::max(vei_i, 1e-300);
//
//  // score = x_i' residual_without_i
//  const double score = r(iu) + wi * b(iu);
//
//  // Same BayesC scalar marginal likelihood as old sbayes(),
//  // but using scalar adjusted residual variance vei_i.
//  const double denom = std::max(vei_safe + wi * vb, 1e-300);
//
//  const double logBF =
//   0.5 * std::log(vei_safe / denom)
//   + 0.5 * score * score * vb / (vei_safe * denom);
//
//  const double logp1 = std::log(pi1) + logBF;
//  const double logp0 = std::log(pi0);
//
//  const double delta_log = logp0 - logp1;
//
//  double p1 = 0.0;
//
//  if (delta_log > 35.0) {
//   p1 = 0.0;
//  } else if (delta_log < -35.0) {
//   p1 = 1.0;
//  } else {
//   p1 = 1.0 / (1.0 + std::exp(delta_log));
//  }
//
//  const int di = (runif(gen) < p1) ? 1 : 0;
//
//  double b_new = 0.0;
//
//  if (di == 1) {
//   const double lhs = wi + vei_safe / vb;
//   const double mean = score / lhs;
//   const double sd = std::sqrt(vei_safe / lhs);
//
//   b_new = mean + sd * norm01(gen);
//  }
//
//  const double diff = b_new - b(iu);
//
//  if (diff != 0.0) {
//   r(iu) -= wi * diff;
//
//   const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
//   const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];
//
//   for (uint64_t p = start; p < end; ++p) {
//    const int j = ld.idx[static_cast<std::size_t>(p)];
//    r(static_cast<arma::uword>(j)) -=
//     static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
//   }
//  }
//
//  b(iu) = b_new;
//  d(iu) = di;
// }
//
// // -----------------------------------------------------------------------------
// // Single-trait variance updates
// // -----------------------------------------------------------------------------
//
// inline void sampleB_ST_csr(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   double ssb_prior,
//   std::mt19937& gen
// ) {
//  double ssb = 0.0;
//  double dfb = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//   if (d(iu) > 0) {
//    ssb += b(iu) * b(iu);
//    dfb += 1.0;
//   }
//  }
//
//  const double scale = ssb + nub * ssb_prior;
//
//  std::chi_squared_distribution<double> rchisq(dfb + nub);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  vb = std::max(scale / chi2, 1e-12);
// }
//
//
// inline void samplePi_ST(
//   const arma::Row<int>& d,
//   std::vector<double>& pi,
//   std::mt19937& gen
// ) {
//  double c0 = 1.0;
//  double c1 = 1.0;
//
//  for (arma::uword i = 0; i < d.n_elem; ++i) {
//   if (d(i) > 0) c1 += 1.0;
//   else c0 += 1.0;
//  }
//
//  std::gamma_distribution<double> rg0(c0, 1.0);
//  std::gamma_distribution<double> rg1(c1, 1.0);
//
//  const double g0 = std::max(rg0(gen), 1e-300);
//  const double g1 = std::max(rg1(gen), 1e-300);
//  const double s = g0 + g1;
//
//  pi[0] = g0 / s;
//  pi[1] = g1 / s;
// }
//
// // -----------------------------------------------------------------------------
// // Main exported function: parallel single-trait BayesC over traits
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr(
//   std::vector<std::vector<double>> wy,
//   std::vector<std::vector<double>> ww,
//   std::vector<double> yy,
//   std::vector<std::vector<double>> b_init,
//   std::vector<std::vector<double>> d_init,
//   bool use_d_init,
//   std::vector<std::vector<double>> r_init,
//   bool use_r_init,
//   bool rebuild_r_before_updateE,
//   std::string ld_prefix,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   std::vector<int> n,
//   int nit,
//   int nburn,
//   int nthin,
//   int ncores,
//   int seed
// ) {
//  const int nt = static_cast<int>(wy.size());
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nt must be positive.");
//  }
//
//  const int m = static_cast<int>(wy[0].size());
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: m must be positive.");
//  }
//
//  if (nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nit must be positive.");
//  }
//
//  if (nburn < 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nburn must be non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nthin must be positive.");
//  }
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt ||
//      (int)yy.size() != nt || (int)n.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: inconsistent trait dimensions.");
//  }
//
//  if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: priors must be nt x nt.");
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_csr: pi must have length 2, c(pi0, pi1).");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m ||
//       (int)ww[t].size() != m ||
//       (int)b_init[t].size() != m) {
//    throw std::runtime_error("stblr_cpg_omp_csr: inconsistent marker dimensions.");
//   }
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: E must be nt x nt.");
//  }
//
//  if (use_r_init) {
//   if (static_cast<int>(r_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_state: r_init must have length nt when use_r_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(r_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_state: each r_init[t] must have length m.");
//    }
//   }
//  }
//
//  if (use_d_init) {
//   if (static_cast<int>(d_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_state: d_init must have length nt when use_d_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(d_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_state: each d_init[t] must have length m.");
//    }
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Convert inputs to Armadillo
//  // --------------------------------------------------------------------------
//
//  arma::mat wy_mat(nt, m, arma::fill::zeros);
//  arma::mat ww_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::mat r_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
//
//  arma::vec yy_vec(nt, arma::fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   yy_vec(static_cast<arma::uword>(t)) = yy[t];
//
//   for (int i = 0; i < m; ++i) {
//    wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
//    ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
//    b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i))  = b_init[t][i];
//   }
//
//   if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
//    sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Validate shared scaling and build shared flat LD object
//  // --------------------------------------------------------------------------
//
//  for (int t = 1; t < nt; ++t) {
//   if (n[t] != n[0]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_csr: current shared-LD scaling assumes equal n across traits."
//    );
//   }
//  }
//
//  for (int t = 1; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    const double w0 = ww_mat(0, static_cast<arma::uword>(i));
//    const double wt = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    const double tol = 1e-8 * std::max(1.0, std::abs(w0));
//
//    if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr: ww contains invalid value before LD pre-scaling."
//     );
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr: ww differs across traits; pre-scaled shared ST LD is invalid."
//     );
//    }
//   }
//  }
//
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//
//  for (int i = 0; i < m; ++i) {
//   const double wi = ww_mat(0, static_cast<arma::uword>(i));
//   if (!std::isfinite(wi) || wi <= 0.0) {
//    throw std::runtime_error("stblr_cpg_omp_csr: ww contains invalid value in trait 0.");
//   }
//   xx[static_cast<std::size_t>(i)] = wi;
//  }
//
//  STLDCSR ld = read_and_build_st_ld_csr(
//   ld_prefix,
//   m,
//   xx
//  );
//
//  // --------------------------------------------------------------------------
//  // Marker update order based on max single-trait marginal effect
//  // --------------------------------------------------------------------------
//
//  std::vector<double> x2(static_cast<std::size_t>(m), 0.0);
//  std::vector<int> order(static_cast<std::size_t>(m));
//
//  for (int i = 0; i < m; ++i) {
//   double best = 0.0;
//
//   for (int t = 0; t < nt; ++t) {
//    const double wi = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    if (wi > 0.0) {
//     const double bhat = wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) / wi;
//     best = std::max(best, bhat * bhat);
//    }
//   }
//
//   x2[static_cast<std::size_t>(i)] = best;
//   order[static_cast<std::size_t>(i)] = i;
//  }
//
//  std::sort(order.begin(), order.end(),
//            [&](int a, int b) { return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)]; });
//
//  // --------------------------------------------------------------------------
//  // Output storage
//  // --------------------------------------------------------------------------
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//
//  arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
//
//  arma::vec final_vb(nt, arma::fill::zeros);
//  arma::vec final_vg(nt, arma::fill::zeros);
//  arma::vec final_ve(nt, arma::fill::zeros);
//  arma::vec final_pi(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  // --------------------------------------------------------------------------
//  // Parallel over traits
//  // --------------------------------------------------------------------------
//
//  std::vector<int> failed(static_cast<std::size_t>(nt), 0);
//  std::vector<std::string> errors(static_cast<std::size_t>(nt));
//  std::vector<int> thread_used(static_cast<std::size_t>(nt), 0);
//  std::vector<double> trait_seconds(static_cast<std::size_t>(nt), 0.0);
//
//  int nthreads = 1;
//
// #ifdef _OPENMP
//  omp_set_dynamic(0);
//
//  nthreads = std::max(1, std::min(ncores, nt));
//
//  // Explicitly set the OpenMP thread count for this runtime.
//  omp_set_num_threads(nthreads);
//
//  Rcpp::Rcout
//  << "STBLR OpenMP requested threads = "
//  << nthreads
//  << ", omp_get_max_threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
// // #ifdef _OPENMP
// //  omp_set_dynamic(0);
// //  nthreads = std::max(1, std::min(ncores, nt));
// //  nthreads = std::min(nthreads, omp_get_max_threads());
// //
// //  Rcpp::Rcout
// //  << "STBLR OpenMP threads = "
// //  << nthreads
// //  << ", max threads = "
// //  << omp_get_max_threads()
// //  << ", num procs = "
// //  << omp_get_num_procs()
// //  << "\n";
// // #endif
//
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(nthreads) schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
//
// #ifdef _OPENMP
//   const double wall_start = omp_get_wtime();
//   thread_used[static_cast<std::size_t>(t)] = omp_get_thread_num();
// #else
//   const double wall_start = 0.0;
//   thread_used[static_cast<std::size_t>(t)] = 0;
// #endif
//
//   try {
//    std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//    arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
//    arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));
//
//    arma::rowvec b_t(m, arma::fill::zeros);
//    for (int i = 0; i < m; ++i) {
//     b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    }
//
//    arma::rowvec r_t(m, arma::fill::zeros);
//    arma::Row<int> d_t(m, arma::fill::zeros);
//
//    if (use_d_init) {
//     for (int i = 0; i < m; ++i) {
//      d_t(static_cast<arma::uword>(i)) = d_init[t][i] > 0 ? 1 : 0;
//     }
//    } else {
//     for (int i = 0; i < m; ++i) {
//      d_t(static_cast<arma::uword>(i)) =
//       (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
//     }
//    }
//
//    if (use_r_init) {
//     for (int i = 0; i < m; ++i) {
//      r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
//     }
//
//     if (!r_t.is_finite()) {
//      throw std::runtime_error("stblr_cpg_omp_csr_state: r_init contains NaN/Inf.");
//     }
//    } else {
//     rebuild_residual_st_csr(
//      m,
//      wy_t,
//      ww_t,
//      b_t,
//      r_t,
//      ld
//     );
//    }
//
//    double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double vg_t = 0.0;
//    double vei_t = ve_t + adjE * vg_t;
//
//    std::vector<double> pi_t = pi;
//
//    if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//        pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//     throw std::runtime_error(
//       "invalid initial pi: pi0=" + std::to_string(pi_t[0]) +
//        ", pi1=" + std::to_string(pi_t[1])
//     );
//    }
//
//    {
//     const double psum = pi_t[0] + pi_t[1];
//
//     if (!std::isfinite(psum) || psum <= 0.0) {
//      throw std::runtime_error("invalid initial pi sum.");
//     }
//
//     pi_t[0] /= psum;
//     pi_t[1] /= psum;
//    }
//
//    arma::rowvec bm_t(m, arma::fill::zeros);
//    arma::rowvec dm_t(m, arma::fill::zeros);
//    arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//
//     // -------------------------------------------------------
//     // Marker updates
//     // -------------------------------------------------------
//     for (int isort = 0; isort < m; ++isort) {
//      const int i = order[static_cast<std::size_t>(isort)];
//
//      sampleBetaC_ST_csr(
//       i,
//       pi_t,
//       vb_t,
//       vei_t,
//       ww_t,
//       r_t,
//       b_t,
//       d_t,
//       ld,
//       gen_t
//      );
//     }
//
//     // -------------------------------------------------------
//     // Variance updates
//     // -------------------------------------------------------
//     if (updateB) {
//      sampleB_ST_csr(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
//       gen_t
//      );
//
//      if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//       throw std::runtime_error(
//         "vb became invalid after sampleB. iter=" +
//          std::to_string(it) +
//          ", vb=" + std::to_string(vb_t)
//       );
//      }
//     }
//
//     if (updateE) {
//      if (rebuild_r_before_updateE) {
//       rebuild_residual_st_csr(
//        m,
//        wy_t,
//        ww_t,
//        b_t,
//        r_t,
//        ld
//       );
//      }
//
//      sampleE_ST_csr(
//       m,
//       nue,
//       ve_t,
//       b_t,
//       wy_t,
//       r_t,
//       sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
//       yy_vec(static_cast<arma::uword>(t)),
//       n[t],
//        gen_t
//      );
//     }
//
//     if (updatePi) {
//      samplePi_ST(d_t, pi_t, gen_t);
//
//      if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//          pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//       throw std::runtime_error(
//         "pi became invalid after samplePi. iter=" +
//          std::to_string(it) +
//          ", pi0=" + std::to_string(pi_t[0]) +
//          ", pi1=" + std::to_string(pi_t[1])
//       );
//      }
//     }
//
//     vg_t = computeG_ST_csr(
//      b_t,
//      wy_t,
//      r_t,
//      n[t]
//     );
//
//     if (!std::isfinite(vg_t)) {
//      throw std::runtime_error(
//        "vg became NaN/Inf after computeG. iter=" +
//         std::to_string(it)
//      );
//     }
//
//     vei_t = ve_t + adjE * vg_t;
//
//     if (!std::isfinite(vei_t) || vei_t <= 0.0) {
//      throw std::runtime_error(
//        "adjusted residual variance vei became invalid. iter=" +
//         std::to_string(it) +
//         ", vei=" + std::to_string(vei_t)
//      );
//     }
//
//     vbs_t(static_cast<arma::uword>(it)) = vb_t;
//     ves_t(static_cast<arma::uword>(it)) = ve_t;
//     vgs_t(static_cast<arma::uword>(it)) = vg_t;
//     pis_t(static_cast<arma::uword>(it)) = pi_t[1];
//
//     // -------------------------------------------------------
//     // Store posterior summaries
//     // -------------------------------------------------------
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       const arma::uword iu = static_cast<arma::uword>(i);
//       bm_t(iu) += b_t(iu);
//       dm_t(iu) += static_cast<double>(d_t(iu));
//      }
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    bm_t /= nsamples_t;
//    dm_t /= nsamples_t;
//
//    if (!bm_t.is_finite()) {
//     throw std::runtime_error("posterior mean bm contains NaN/Inf.");
//    }
//
//    if (!dm_t.is_finite()) {
//     throw std::runtime_error("posterior mean dm contains NaN/Inf.");
//    }
//
//    bm_mat.row(static_cast<arma::uword>(t)) = bm_t;
//    dm_mat.row(static_cast<arma::uword>(t)) = dm_t;
//    b_mat.row(static_cast<arma::uword>(t))  = b_t;
//    r_mat.row(static_cast<arma::uword>(t))  = r_t;
//    d_mat.row(static_cast<arma::uword>(t))  = d_t;
//
//    vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
//    vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
//    ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
//    pis_mat.row(static_cast<arma::uword>(t)) = pis_t;
//
//    final_vb(static_cast<arma::uword>(t)) = vb_t;
//    final_ve(static_cast<arma::uword>(t)) = ve_t;
//    final_vg(static_cast<arma::uword>(t)) = vg_t;
//    final_pi(static_cast<arma::uword>(t)) = pi_t[1];
//    nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;
//
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//
//   } catch (const std::exception& e) {
//    failed[static_cast<std::size_t>(t)] = 1;
//    errors[static_cast<std::size_t>(t)] = e.what();
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//   } catch (...) {
//    failed[static_cast<std::size_t>(t)] = 1;
//    errors[static_cast<std::size_t>(t)] = "unknown error";
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//   }
//  }
//
// #ifdef _OPENMP
//  for (int t = 0; t < nt; ++t) {
//   Rcpp::Rcout
//   << "trait " << t
//   << " used thread " << thread_used[static_cast<std::size_t>(t)]
//   << ", seconds = " << trait_seconds[static_cast<std::size_t>(t)]
//   << "\n";
//  }
// #endif
//
//  for (int t = 0; t < nt; ++t) {
//   if (failed[static_cast<std::size_t>(t)]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_csr failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[static_cast<std::size_t>(t)]
//    );
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Build result with same style as MT output
//  // --------------------------------------------------------------------------
//
//  std::vector<std::vector<std::vector<double>>> result(20);
//
//  for (int k = 0; k < 20; ++k) {
//   result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // bm
//   result[1][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // dm
//   result[2][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // wy
//   result[3][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // r
//   result[4][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // b
//   result[5][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // d
//   result[6][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // marker index
//
//   result[7][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vbs
//   result[8][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vgs
//   result[9][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // ves
//
//   result[10][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // covb
//   result[11][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // covg
//   result[12][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // cove
//   result[13][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final B
//   result[14][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final G
//   result[15][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final E
//
//   result[16][static_cast<std::size_t>(t)].resize(2);                                     // final pi
//   result[17][static_cast<std::size_t>(t)].resize(2);                                     // posterior mean pi approx
//
//   result[18][static_cast<std::size_t>(t)].resize(4);
//   result[19][static_cast<std::size_t>(t)].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int i = 0; i < m; ++i) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword iu = static_cast<arma::uword>(i);
//    const std::size_t is = static_cast<std::size_t>(i);
//
//    result[0][ts][is] = bm_mat(tu, iu);
//    result[1][ts][is] = dm_mat(tu, iu);
//    result[2][ts][is] = wy_mat(tu, iu);
//    result[3][ts][is] = r_mat(tu, iu);
//    result[4][ts][is] = b_mat(tu, iu);
//    result[5][ts][is] = static_cast<double>(d_mat(tu, iu));
//    result[6][ts][is] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int it = 0; it < nit + nburn; ++it) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword itu = static_cast<arma::uword>(it);
//    const std::size_t its = static_cast<std::size_t>(it);
//
//    result[7][ts][its] = vbs_mat(tu, itu);
//    result[8][ts][its] = vgs_mat(tu, itu);
//    result[9][ts][its] = ves_mat(tu, itu);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   const std::size_t t1s = static_cast<std::size_t>(t1);
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    const std::size_t t2s = static_cast<std::size_t>(t2);
//
//    result[10][t1s][t2s] = 0.0;
//    result[11][t1s][t2s] = 0.0;
//    result[12][t1s][t2s] = 0.0;
//
//    result[13][t1s][t2s] = 0.0;
//    result[14][t1s][t2s] = 0.0;
//    result[15][t1s][t2s] = 0.0;
//   }
//
//   result[10][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
//   result[11][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
//   result[12][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
//
//   result[13][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
//   result[14][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
//   result[15][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   result[16][ts][0] = 1.0 - final_pi(static_cast<arma::uword>(t));
//   result[16][ts][1] = final_pi(static_cast<arma::uword>(t));
//
//   // Approximate posterior mean inclusion pi from saved trace.
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi(static_cast<arma::uword>(t));
//
//   result[17][ts][0] = 1.0 - mean_pi;
//   result[17][ts][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) {
//    result[18][ts][static_cast<std::size_t>(i)] = 0.0;
//   }
//
//   for (int i = 0; i < 2; ++i) {
//    result[19][ts][static_cast<std::size_t>(i)] = 0.0;
//   }
//  }
//
//  return result;
// }
//
//
//
// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
//
// #include "cpg_samplers.h"
// #include "distributions.h"
// #include "st_csr_common.h"
//
// #include <algorithm>
// #include <cmath>
// #include <cstdio>
// #include <cstdint>
// #include <fstream>
// #include <limits>
// #include <numeric>
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
// #include <unordered_map>
// #include <unordered_set>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // -----------------------------------------------------------------------------
// // ST-specific LD structure: flat symmetric CSR
// // Stores pre-scaled X_i'X_j, not raw LD correlation.
// // Disk input is expected to be upper-triangular or otherwise non-symmetric CSR
// // with raw correlations r_ij. This builder symmetrizes it.
// // -----------------------------------------------------------------------------
//
// // -----------------------------------------------------------------------------
// // Static LD-risk correction
// // -----------------------------------------------------------------------------
// //
// // The sparse-LD approximation error for marker i can accumulate when many
// // omitted or retained LD links have coherent sign. This is not necessarily
// // captured by local LD score alone. Therefore, we allow the caller to pass a
// // precomputed marker-level LD-risk vector q_i, for example
// //
// //   q_i = |sum_j r_ij| / sum_j |r_ij| * log(1 + sum_j |r_ij|)
// //
// // scaled to median 1.
// //
// // The marker-specific effective residual variance is then
// //
// //   vei_i = ve + adjE * vg + ld_error_scale * vg * q_i,
// //
// // optionally capped relative to the base variance.
// //
// // If ld_error_scale = 0 or ld_risk is empty, this reduces exactly to the old
// // behavior.
// // -----------------------------------------------------------------------------
//
// static inline double marker_specific_vei_ST_csr(
//   int i,
//   double base_vei,
//   double vg,
//   const std::vector<double>& ld_risk,
//   double ld_error_scale,
//   double max_ld_inflation
// ) {
//  const double vei_safe = std::max(base_vei, 1e-300);
//
//  if (ld_error_scale <= 0.0 || ld_risk.empty() || vg <= 0.0) {
//   return vei_safe;
//  }
//
//  const double q = std::max(0.0, ld_risk[static_cast<std::size_t>(i)]);
//
//  if (!std::isfinite(q)) {
//   throw std::runtime_error("ld_risk contains NaN/Inf.");
//  }
//
//  double extra = ld_error_scale * vg * q;
//
//  if (!std::isfinite(extra) || extra < 0.0) {
//   extra = 0.0;
//  }
//
//  if (max_ld_inflation >= 0.0) {
//   const double cap = max_ld_inflation * vei_safe;
//   extra = std::min(extra, cap);
//  }
//
//  return vei_safe + extra;
// }
//
// static std::vector<double> copy_ld_risk_or_empty(
//   Rcpp::Nullable<Rcpp::NumericVector> ld_risk,
//   int m
// ) {
//  std::vector<double> out;
//
//  if (ld_risk.isNotNull()) {
//   Rcpp::NumericVector x(ld_risk);
//
//   if (x.size() == 0) {
//    return out;
//   }
//
//   if (x.size() != m) {
//    throw std::runtime_error("ld_risk must be NULL, empty, or have length m.");
//   }
//
//   out.resize(static_cast<std::size_t>(m));
//
//   for (int i = 0; i < m; ++i) {
//    const double qi = x[i];
//
//    if (!std::isfinite(qi) || qi < 0.0) {
//     throw std::runtime_error("ld_risk must contain finite non-negative values.");
//    }
//
//    out[static_cast<std::size_t>(i)] = qi;
//   }
//  }
//
//  return out;
// }
//
// // -----------------------------------------------------------------------------
// // Single-trait BayesC marker update
// // -----------------------------------------------------------------------------
//
// inline void sampleBetaC_ST_csr(
//   int i,
//   const std::vector<double>& pi,
//   double vb,
//   double base_vei,
//   double vg,
//   const std::vector<double>& ld_risk,
//   double ld_error_scale,
//   double max_ld_inflation,
//   const arma::rowvec& ww,
//   arma::rowvec& r,
//   arma::rowvec& b,
//   arma::Row<int>& d,
//   const STLDCSR& ld,
//   std::mt19937& gen
// ) {
//  const arma::uword iu = static_cast<arma::uword>(i);
//  const double wi = ww(iu);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm01(0.0, 1.0);
//
//  const double pi0 = std::max(pi[0], 1e-300);
//  const double pi1 = std::max(pi[1], 1e-300);
//
//  // score = x_i' residual_without_i
//  const double score = r(iu) + wi * b(iu);
//
//  // Old behavior: vei_i = base_vei.
//  // Optional static sparse-LD error correction: add ld_error_scale * vg * q_i.
//  const double vei_safe = marker_specific_vei_ST_csr(
//   i,
//   base_vei,
//   vg,
//   ld_risk,
//   ld_error_scale,
//   max_ld_inflation
//  );
//
//  const double denom = std::max(vei_safe + wi * vb, 1e-300);
//
//  const double logBF =
//   0.5 * std::log(vei_safe / denom)
//   + 0.5 * score * score * vb / (vei_safe * denom);
//
//  const double logp1 = std::log(pi1) + logBF;
//  const double logp0 = std::log(pi0);
//
//  const double delta_log = logp0 - logp1;
//
//  double p1 = 0.0;
//
//  if (delta_log > 35.0) {
//   p1 = 0.0;
//  } else if (delta_log < -35.0) {
//   p1 = 1.0;
//  } else {
//   p1 = 1.0 / (1.0 + std::exp(delta_log));
//  }
//
//  const int di = (runif(gen) < p1) ? 1 : 0;
//
//  double b_new = 0.0;
//
//  if (di == 1) {
//   const double lhs = wi + vei_safe / vb;
//   const double mean = score / lhs;
//   const double sd = std::sqrt(vei_safe / lhs);
//
//   b_new = mean + sd * norm01(gen);
//  }
//
//  const double diff = b_new - b(iu);
//
//  if (diff != 0.0) {
//   r(iu) -= wi * diff;
//
//   const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
//   const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];
//
//   for (uint64_t p = start; p < end; ++p) {
//    const int j = ld.idx[static_cast<std::size_t>(p)];
//    r(static_cast<arma::uword>(j)) -=
//     static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
//   }
//  }
//
//  b(iu) = b_new;
//  d(iu) = di;
// }
//
// // -----------------------------------------------------------------------------
// // Single-trait variance updates
// // -----------------------------------------------------------------------------
//
// inline void sampleB_ST_csr(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   double ssb_prior,
//   std::mt19937& gen
// ) {
//  double ssb = 0.0;
//  double dfb = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//   if (d(iu) > 0) {
//    ssb += b(iu) * b(iu);
//    dfb += 1.0;
//   }
//  }
//
//  const double scale = ssb + nub * ssb_prior;
//
//  std::chi_squared_distribution<double> rchisq(dfb + nub);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  vb = std::max(scale / chi2, 1e-12);
// }
//
// inline void samplePi_ST(
//   const arma::Row<int>& d,
//   std::vector<double>& pi,
//   std::mt19937& gen
// ) {
//  double c0 = 1.0;
//  double c1 = 1.0;
//
//  for (arma::uword i = 0; i < d.n_elem; ++i) {
//   if (d(i) > 0) c1 += 1.0;
//   else c0 += 1.0;
//  }
//
//  std::gamma_distribution<double> rg0(c0, 1.0);
//  std::gamma_distribution<double> rg1(c1, 1.0);
//
//  const double g0 = std::max(rg0(gen), 1e-300);
//  const double g1 = std::max(rg1(gen), 1e-300);
//  const double s = g0 + g1;
//
//  pi[0] = g0 / s;
//  pi[1] = g1 / s;
// }
//
// // -----------------------------------------------------------------------------
// // Main exported function: parallel single-trait BayesC over traits
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr(
//   std::vector<std::vector<double>> wy,
//   std::vector<std::vector<double>> ww,
//   std::vector<double> yy,
//   std::vector<std::vector<double>> b_init,
//   std::vector<std::vector<double>> d_init,
//   bool use_d_init,
//   std::vector<std::vector<double>> r_init,
//   bool use_r_init,
//   bool rebuild_r_before_updateE,
//   std::string ld_prefix,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   std::vector<int> n,
//   int nit,
//   int nburn,
//   int nthin,
//   int ncores,
//   int seed,
//   Rcpp::Nullable<Rcpp::NumericVector> ld_risk = R_NilValue,
//   double ld_error_scale = 0.0,
//   double max_ld_inflation = 5.0,
//   std::string update_order = "marginal"
// ) {
//  const int nt = static_cast<int>(wy.size());
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nt must be positive.");
//  }
//
//  const int m = static_cast<int>(wy[0].size());
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: m must be positive.");
//  }
//
//  if (nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nit must be positive.");
//  }
//
//  if (nburn < 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nburn must be non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr: nthin must be positive.");
//  }
//
//  if (!std::isfinite(ld_error_scale) || ld_error_scale < 0.0) {
//   throw std::runtime_error("ld_error_scale must be finite and non-negative.");
//  }
//
//  if (!std::isfinite(max_ld_inflation) || max_ld_inflation < 0.0) {
//   throw std::runtime_error("max_ld_inflation must be finite and non-negative.");
//  }
//
//  if (update_order != "marginal" && update_order != "natural" && update_order != "random") {
//   throw std::runtime_error("update_order must be 'marginal', 'natural', or 'random'.");
//  }
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt ||
//      (int)yy.size() != nt || (int)n.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: inconsistent trait dimensions.");
//  }
//
//  if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: priors must be nt x nt.");
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_csr: pi must have length 2, c(pi0, pi1).");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m ||
//       (int)ww[t].size() != m ||
//       (int)b_init[t].size() != m) {
//    throw std::runtime_error("stblr_cpg_omp_csr: inconsistent marker dimensions.");
//   }
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr: E must be nt x nt.");
//  }
//
//  if (use_r_init) {
//   if (static_cast<int>(r_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_state: r_init must have length nt when use_r_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(r_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_state: each r_init[t] must have length m.");
//    }
//   }
//  }
//
//  if (use_d_init) {
//   if (static_cast<int>(d_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_state: d_init must have length nt when use_d_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(d_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_state: each d_init[t] must have length m.");
//    }
//   }
//  }
//
//  const std::vector<double> ld_risk_vec = copy_ld_risk_or_empty(ld_risk, m);
//
//  // --------------------------------------------------------------------------
//  // Convert inputs to Armadillo
//  // --------------------------------------------------------------------------
//
//  arma::mat wy_mat(nt, m, arma::fill::zeros);
//  arma::mat ww_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::mat r_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
//
//  arma::vec yy_vec(nt, arma::fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   yy_vec(static_cast<arma::uword>(t)) = yy[t];
//
//   for (int i = 0; i < m; ++i) {
//    wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
//    ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
//    b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i))  = b_init[t][i];
//   }
//
//   if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
//    sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Validate shared scaling and build shared flat LD object
//  // --------------------------------------------------------------------------
//
//  for (int t = 1; t < nt; ++t) {
//   if (n[t] != n[0]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_csr: current shared-LD scaling assumes equal n across traits."
//    );
//   }
//  }
//
//  for (int t = 1; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    const double w0 = ww_mat(0, static_cast<arma::uword>(i));
//    const double wt = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    const double tol = 1e-8 * std::max(1.0, std::abs(w0));
//
//    if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr: ww contains invalid value before LD pre-scaling."
//     );
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr: ww differs across traits; pre-scaled shared ST LD is invalid."
//     );
//    }
//   }
//  }
//
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//
//  for (int i = 0; i < m; ++i) {
//   const double wi = ww_mat(0, static_cast<arma::uword>(i));
//   if (!std::isfinite(wi) || wi <= 0.0) {
//    throw std::runtime_error("stblr_cpg_omp_csr: ww contains invalid value in trait 0.");
//   }
//   xx[static_cast<std::size_t>(i)] = wi;
//  }
//
//  STLDCSR ld = read_and_build_st_ld_csr(
//   ld_prefix,
//   m,
//   xx
//  );
//
//  // --------------------------------------------------------------------------
//  // Marker update order
//  // --------------------------------------------------------------------------
//
//  std::vector<double> x2(static_cast<std::size_t>(m), 0.0);
//  std::vector<int> order(static_cast<std::size_t>(m));
//
//  for (int i = 0; i < m; ++i) {
//   double best = 0.0;
//
//   for (int t = 0; t < nt; ++t) {
//    const double wi = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    if (wi > 0.0) {
//     const double bhat = wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) / wi;
//     best = std::max(best, bhat * bhat);
//    }
//   }
//
//   x2[static_cast<std::size_t>(i)] = best;
//   order[static_cast<std::size_t>(i)] = i;
//  }
//
//  if (update_order == "marginal") {
//   std::sort(order.begin(), order.end(),
//             [&](int a, int b) { return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)]; });
//  }
//
//  // --------------------------------------------------------------------------
//  // Output storage
//  // --------------------------------------------------------------------------
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//
//  arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
//
//  arma::vec final_vb(nt, arma::fill::zeros);
//  arma::vec final_vg(nt, arma::fill::zeros);
//  arma::vec final_ve(nt, arma::fill::zeros);
//  arma::vec final_pi(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  // --------------------------------------------------------------------------
//  // Parallel over traits
//  // --------------------------------------------------------------------------
//
//  std::vector<int> failed(static_cast<std::size_t>(nt), 0);
//  std::vector<std::string> errors(static_cast<std::size_t>(nt));
//  std::vector<int> thread_used(static_cast<std::size_t>(nt), 0);
//  std::vector<double> trait_seconds(static_cast<std::size_t>(nt), 0.0);
//
//  int nthreads = 1;
//
// #ifdef _OPENMP
//  omp_set_dynamic(0);
//
//  nthreads = std::max(1, std::min(ncores, nt));
//
//  // Explicitly set the OpenMP thread count for this runtime.
//  omp_set_num_threads(nthreads);
//
//  Rcpp::Rcout
//  << "STBLR OpenMP requested threads = "
//  << nthreads
//  << ", omp_get_max_threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
//  Rcpp::Rcout
//  << "STBLR CSR static LD-risk correction: ld_error_scale=" << ld_error_scale
//  << ", max_ld_inflation=" << max_ld_inflation
//  << ", ld_risk_length=" << ld_risk_vec.size()
//  << ", update_order=" << update_order
//  << "\n";
//
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(nthreads) schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
//
// #ifdef _OPENMP
//   const double wall_start = omp_get_wtime();
//   thread_used[static_cast<std::size_t>(t)] = omp_get_thread_num();
// #else
//   const double wall_start = 0.0;
//   thread_used[static_cast<std::size_t>(t)] = 0;
// #endif
//
//   try {
//    std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//    arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
//    arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));
//
//    arma::rowvec b_t(m, arma::fill::zeros);
//    for (int i = 0; i < m; ++i) {
//     b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//    }
//
//    arma::rowvec r_t(m, arma::fill::zeros);
//    arma::Row<int> d_t(m, arma::fill::zeros);
//
//    if (use_d_init) {
//     for (int i = 0; i < m; ++i) {
//      d_t(static_cast<arma::uword>(i)) = d_init[t][i] > 0 ? 1 : 0;
//     }
//    } else {
//     for (int i = 0; i < m; ++i) {
//      d_t(static_cast<arma::uword>(i)) =
//       (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
//     }
//    }
//
//    if (use_r_init) {
//     for (int i = 0; i < m; ++i) {
//      r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
//     }
//
//     if (!r_t.is_finite()) {
//      throw std::runtime_error("stblr_cpg_omp_csr_state: r_init contains NaN/Inf.");
//     }
//    } else {
//     rebuild_residual_st_csr(
//      m,
//      wy_t,
//      ww_t,
//      b_t,
//      r_t,
//      ld
//     );
//    }
//
//    double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double vg_t = computeG_ST_csr(
//     b_t,
//     wy_t,
//     r_t,
//     n[t]
//    );
//    double vei_t = ve_t + adjE * vg_t;
//
//    std::vector<double> pi_t = pi;
//
//    if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//        pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//     throw std::runtime_error(
//       "invalid initial pi: pi0=" + std::to_string(pi_t[0]) +
//        ", pi1=" + std::to_string(pi_t[1])
//     );
//    }
//
//    {
//     const double psum = pi_t[0] + pi_t[1];
//
//     if (!std::isfinite(psum) || psum <= 0.0) {
//      throw std::runtime_error("invalid initial pi sum.");
//     }
//
//     pi_t[0] /= psum;
//     pi_t[1] /= psum;
//    }
//
//    arma::rowvec bm_t(m, arma::fill::zeros);
//    arma::rowvec dm_t(m, arma::fill::zeros);
//    arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//
//     // -------------------------------------------------------
//     // Marker updates
//     // -------------------------------------------------------
//     if (update_order == "random") {
//      std::shuffle(order.begin(), order.end(), gen_t);
//     }
//
//     for (int isort = 0; isort < m; ++isort) {
//      const int i = order[static_cast<std::size_t>(isort)];
//
//      sampleBetaC_ST_csr(
//       i,
//       pi_t,
//       vb_t,
//       vei_t,
//       vg_t,
//       ld_risk_vec,
//       ld_error_scale,
//       max_ld_inflation,
//       ww_t,
//       r_t,
//       b_t,
//       d_t,
//       ld,
//       gen_t
//      );
//     }
//
//     // -------------------------------------------------------
//     // Variance updates
//     // -------------------------------------------------------
//     if (updateB) {
//      sampleB_ST_csr(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
//       gen_t
//      );
//
//      if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//       throw std::runtime_error(
//         "vb became invalid after sampleB. iter=" +
//          std::to_string(it) +
//          ", vb=" + std::to_string(vb_t)
//       );
//      }
//     }
//
//     if (updateE) {
//      if (rebuild_r_before_updateE) {
//       rebuild_residual_st_csr(
//        m,
//        wy_t,
//        ww_t,
//        b_t,
//        r_t,
//        ld
//       );
//      }
//
//      sampleE_ST_csr(
//       m,
//       nue,
//       ve_t,
//       b_t,
//       wy_t,
//       r_t,
//       sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
//       yy_vec(static_cast<arma::uword>(t)),
//       n[t],
//        gen_t
//      );
//     }
//
//     if (updatePi) {
//      samplePi_ST(d_t, pi_t, gen_t);
//
//      if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//          pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//       throw std::runtime_error(
//         "pi became invalid after samplePi. iter=" +
//          std::to_string(it) +
//          ", pi0=" + std::to_string(pi_t[0]) +
//          ", pi1=" + std::to_string(pi_t[1])
//       );
//      }
//     }
//
//     vg_t = computeG_ST_csr(
//      b_t,
//      wy_t,
//      r_t,
//      n[t]
//     );
//
//     if (!std::isfinite(vg_t)) {
//      throw std::runtime_error(
//        "vg became NaN/Inf after computeG. iter=" +
//         std::to_string(it)
//      );
//     }
//
//     vei_t = ve_t + adjE * vg_t;
//
//     if (!std::isfinite(vei_t) || vei_t <= 0.0) {
//      throw std::runtime_error(
//        "adjusted residual variance vei became invalid. iter=" +
//         std::to_string(it) +
//         ", vei=" + std::to_string(vei_t)
//      );
//     }
//
//     vbs_t(static_cast<arma::uword>(it)) = vb_t;
//     ves_t(static_cast<arma::uword>(it)) = ve_t;
//     vgs_t(static_cast<arma::uword>(it)) = vg_t;
//     pis_t(static_cast<arma::uword>(it)) = pi_t[1];
//
//     // -------------------------------------------------------
//     // Store posterior summaries
//     // -------------------------------------------------------
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       const arma::uword iu = static_cast<arma::uword>(i);
//       bm_t(iu) += b_t(iu);
//       dm_t(iu) += static_cast<double>(d_t(iu));
//      }
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    bm_t /= nsamples_t;
//    dm_t /= nsamples_t;
//
//    if (!bm_t.is_finite()) {
//     throw std::runtime_error("posterior mean bm contains NaN/Inf.");
//    }
//
//    if (!dm_t.is_finite()) {
//     throw std::runtime_error("posterior mean dm contains NaN/Inf.");
//    }
//
//    bm_mat.row(static_cast<arma::uword>(t)) = bm_t;
//    dm_mat.row(static_cast<arma::uword>(t)) = dm_t;
//    b_mat.row(static_cast<arma::uword>(t))  = b_t;
//    r_mat.row(static_cast<arma::uword>(t))  = r_t;
//    d_mat.row(static_cast<arma::uword>(t))  = d_t;
//
//    vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
//    vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
//    ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
//    pis_mat.row(static_cast<arma::uword>(t)) = pis_t;
//
//    final_vb(static_cast<arma::uword>(t)) = vb_t;
//    final_ve(static_cast<arma::uword>(t)) = ve_t;
//    final_vg(static_cast<arma::uword>(t)) = vg_t;
//    final_pi(static_cast<arma::uword>(t)) = pi_t[1];
//    nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;
//
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//
//   } catch (const std::exception& e) {
//    failed[static_cast<std::size_t>(t)] = 1;
//    errors[static_cast<std::size_t>(t)] = e.what();
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//   } catch (...) {
//    failed[static_cast<std::size_t>(t)] = 1;
//    errors[static_cast<std::size_t>(t)] = "unknown error";
// #ifdef _OPENMP
//    trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
// #endif
//   }
//  }
//
// #ifdef _OPENMP
//  for (int t = 0; t < nt; ++t) {
//   Rcpp::Rcout
//   << "trait " << t
//   << " used thread " << thread_used[static_cast<std::size_t>(t)]
//   << ", seconds = " << trait_seconds[static_cast<std::size_t>(t)]
//   << "\n";
//  }
// #endif
//
//  for (int t = 0; t < nt; ++t) {
//   if (failed[static_cast<std::size_t>(t)]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_csr failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[static_cast<std::size_t>(t)]
//    );
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Build result with same style as MT output
//  // --------------------------------------------------------------------------
//
//  std::vector<std::vector<std::vector<double>>> result(20);
//
//  for (int k = 0; k < 20; ++k) {
//   result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // bm
//   result[1][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // dm
//   result[2][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // wy
//   result[3][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // r
//   result[4][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // b
//   result[5][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // d
//   result[6][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));             // marker index
//
//   result[7][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vbs
//   result[8][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // vgs
//   result[9][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));   // ves
//
//   result[10][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // covb
//   result[11][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // covg
//   result[12][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // cove
//   result[13][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final B
//   result[14][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final G
//   result[15][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));           // final E
//
//   result[16][static_cast<std::size_t>(t)].resize(2);                                     // final pi
//   result[17][static_cast<std::size_t>(t)].resize(2);                                     // posterior mean pi approx
//
//   result[18][static_cast<std::size_t>(t)].resize(4);
//   result[19][static_cast<std::size_t>(t)].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int i = 0; i < m; ++i) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword iu = static_cast<arma::uword>(i);
//    const std::size_t is = static_cast<std::size_t>(i);
//
//    result[0][ts][is] = bm_mat(tu, iu);
//    result[1][ts][is] = dm_mat(tu, iu);
//    result[2][ts][is] = wy_mat(tu, iu);
//    result[3][ts][is] = r_mat(tu, iu);
//    result[4][ts][is] = b_mat(tu, iu);
//    result[5][ts][is] = static_cast<double>(d_mat(tu, iu));
//    result[6][ts][is] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int it = 0; it < nit + nburn; ++it) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword itu = static_cast<arma::uword>(it);
//    const std::size_t its = static_cast<std::size_t>(it);
//
//    result[7][ts][its] = vbs_mat(tu, itu);
//    result[8][ts][its] = vgs_mat(tu, itu);
//    result[9][ts][its] = ves_mat(tu, itu);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   const std::size_t t1s = static_cast<std::size_t>(t1);
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    const std::size_t t2s = static_cast<std::size_t>(t2);
//
//    result[10][t1s][t2s] = 0.0;
//    result[11][t1s][t2s] = 0.0;
//    result[12][t1s][t2s] = 0.0;
//
//    result[13][t1s][t2s] = 0.0;
//    result[14][t1s][t2s] = 0.0;
//    result[15][t1s][t2s] = 0.0;
//   }
//
//   result[10][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
//   result[11][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
//   result[12][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
//
//   result[13][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
//   result[14][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
//   result[15][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   result[16][ts][0] = 1.0 - final_pi(static_cast<arma::uword>(t));
//   result[16][ts][1] = final_pi(static_cast<arma::uword>(t));
//
//   // Approximate posterior mean inclusion pi from saved trace.
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi(static_cast<arma::uword>(t));
//
//   result[17][ts][0] = 1.0 - mean_pi;
//   result[17][ts][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) {
//    result[18][ts][static_cast<std::size_t>(i)] = 0.0;
//   }
//
//   for (int i = 0; i < 2; ++i) {
//    result[19][ts][static_cast<std::size_t>(i)] = 0.0;
//   }
//  }
//
//  return result;
// }

