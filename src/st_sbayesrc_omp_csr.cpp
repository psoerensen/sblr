// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "st_chain_utils.h"
#include "st_bayesrc_annotation_prior.h"
#include "st_block_eigen.h"
#include "st_block_eigen_rcpp.h"
#include "st_csr_common.h"
#include "st_ld_operator.h"
#include "blr_csr_sbayesrc_types.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
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

static std::vector<std::string> stblr_sbayesrc_copy_character_vector(
  Rcpp::CharacterVector x,
  const char* label
) {
 std::vector<std::string> out(static_cast<std::size_t>(x.size()));
 for (int i = 0; i < x.size(); ++i) {
  if (x[i] == NA_STRING) throw std::runtime_error(std::string(label) + " contains NA.");
  out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(x[i]);
  if (out[static_cast<std::size_t>(i)].empty()) {
   throw std::runtime_error(std::string(label) + " contains an empty string.");
  }
 }
 return out;
}

static std::vector<std::vector<int>> stblr_sbayesrc_copy_int_list(
  Rcpp::List x,
  const char* label
) {
 std::vector<std::vector<int>> out(static_cast<std::size_t>(x.size()));
 for (int k = 0; k < x.size(); ++k) {
  Rcpp::IntegerVector v = x[k];
  out[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(v.size()));
  for (int i = 0; i < v.size(); ++i) {
   if (v[i] == NA_INTEGER || v[i] <= 0) {
    throw std::runtime_error(
     std::string(label) + " must contain positive 1-based marker indices."
    );
   }
   out[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = v[i];
  }
 }
 return out;
}

static std::vector<int> stblr_sbayesrc_copy_rows0_or_empty(
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  int n_bed
) {
 std::vector<int> out;
 if (rows.isNull()) return out;
 Rcpp::IntegerVector r(rows);
 out.resize(static_cast<std::size_t>(r.size()));
 bool identity = r.size() == n_bed;
 for (int i = 0; i < r.size(); ++i) {
  if (r[i] == NA_INTEGER || r[i] < 1 || r[i] > n_bed) {
   throw std::runtime_error("rows contains an index outside [1, n_bed].");
  }
  if (r[i] != i + 1) identity = false;
  out[static_cast<std::size_t>(i)] = r[i] - 1;
 }
 if (identity) out.clear();
 return out;
}


// =============================================================================
// STBLR summary-stat CSR: SBayesRC-style sampler with overlapping annotations
// =============================================================================
//
// Exported function:
//   stblr_cpg_omp_csr_sbayesrc(...)
//
// This is intentionally a full-sweep CSR sampler. No sparse scheduling is used.
//
// Model:
//   gamma = c(0, 0.01, 0.1, 1)             fixed mixture variance multipliers
//   b_i | z_i = k ~ N(0, vb * gamma[k])    for k > 0
//   b_i | z_i = 0 = 0
//
// Annotation model:
//   p_ij = Phi(A_i alpha_j), j = 1,...,Kgamma-1
//
// where p_ij is the conditional probability of stepping from component j-1
// to a larger component.
//
// For Kgamma = 4 components:
//   pi_i0 = 1 - p_i1
//   pi_i1 = p_i1 * (1 - p_i2)
//   pi_i2 = p_i1 * p_i2 * (1 - p_i3)
//   pi_i3 = p_i1 * p_i2 * p_i3
//
// This handles overlapping annotations directly through dense A (m x nAnno).
// If an intercept is desired, include a first column of 1s in A. The first
// annotation coefficient is updated with a flat prior when intercept_flat=true.
//
// Return structure:
//   0  bm
//   1  dm = P(component > 0)
//   2  wy
//   3  r
//   4  b
//   5  component index, 0-based
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
//   16 final active pi, length 2: c(pi0, 1-pi0)
//   17 posterior mean active pi, length 2
//   18 posterior mean alpha, flattened as nAnno x (Kgamma-1), column-major
//   19 posterior mean sigmaSqAlpha, length Kgamma-1
//   20 vle = linkage-equilibrium variance trace
//   21 vld = linkage-disequilibrium contribution trace, vg - vle
//   22 comp_prob = posterior component probabilities, flattened m x Kgamma
//   23 ncomp = posterior mean number of markers per component, length Kgamma
//   24 LD-swap diagnostics, length 4: attempted, accepted, rate, flag
//   25 bm_sd, optional when nchains > 1 or keep_chains = true
//   26 bm_min
//   27 bm_max
//   28 dm_sd
//   29 dm_min
//   30 dm_max
//   31 chain dm, optional when keep_chains = true, chain-major
//   32 chain bm, optional when keep_chains = true, chain-major
//   33 chain LD-swap diagnostics, optional when keep_chains = true, chain-major
//   34 chain comp_prob, optional when keep_chains = true, chain-major
//   35 chain alpha, optional when keep_chains = true, chain-major
//   36 chain sigmaSqAlpha, optional when keep_chains = true, chain-major
//   37 selection_s trace, optional when estimate_selection_s = true
//   38 selection_s acceptance, optional when estimate_selection_s = true
//   39 chain selection_s trace, optional when estimate_selection_s and keep_chains
//   40 chain selection_s acceptance, optional when estimate_selection_s and keep_chains
//
// Recommended R extraction for slot 18:
//   alpha <- matrix(fit[[18]][[t]], nrow=ncol(A), ncol=length(gamma)-1)
//
// =============================================================================

inline double logsumexp_vec(const std::vector<double>& x) {
 double mx = -std::numeric_limits<double>::infinity();
 for (double v : x) mx = std::max(mx, v);
 if (!std::isfinite(mx)) return mx;

 double s = 0.0;
 for (double v : x) s += std::exp(v - mx);
 return mx + std::log(s);
}

inline int sample_categorical_logprob(
  const std::vector<double>& logp,
  std::mt19937& gen
) {
 const double lse = logsumexp_vec(logp);
 if (!std::isfinite(lse)) {
  throw std::runtime_error("sample_categorical_logprob: invalid log probability vector.");
 }

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 const double u = runif(gen);
 double cum = 0.0;

 for (std::size_t k = 0; k < logp.size(); ++k) {
  cum += std::exp(logp[k] - lse);
  if (u <= cum) return static_cast<int>(k);
 }

 return static_cast<int>(logp.size() - 1);
}

template <class OpT>
inline void sampleBeta_SBayesRC_ST_csr(
  int i,
  const arma::rowvec& pi_i,
  const arma::vec& gamma,
  double vb_t,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const OpT& op,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double wi = ww(iu);

 if (!std::isfinite(wi) || wi <= 0.0) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: invalid ww value.");
 }

 if (!std::isfinite(vb_t) || vb_t <= 0.0) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: invalid vb.");
 }

 const int Kgamma = static_cast<int>(gamma.n_elem);
 if (static_cast<int>(pi_i.n_elem) != Kgamma) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: pi_i/gamma length mismatch.");
 }

 const double vei_safe = std::max(vei_i, 1e-300);
 const double score = r(iu) + wi * b(iu);

 std::vector<double> logp(static_cast<std::size_t>(Kgamma));

 for (int k = 0; k < Kgamma; ++k) {
  const double pik = std::max(static_cast<double>(pi_i(static_cast<arma::uword>(k))), 1e-300);
  const double gk = gamma(static_cast<arma::uword>(k));

  if (!std::isfinite(gk) || gk < 0.0) {
   throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: gamma must be non-negative.");
  }

  if (gk <= 0.0) {
   logp[static_cast<std::size_t>(k)] = std::log(pik);
  } else {
   const double vbk = std::max(vb_t * gk, 1e-300);
   const double denom = std::max(vei_safe + wi * vbk, 1e-300);

   const double logBF =
    0.5 * std::log(vei_safe / denom)
    + 0.5 * score * score * vbk / (vei_safe * denom);

   logp[static_cast<std::size_t>(k)] = std::log(pik) + logBF;
  }
 }

 const int k_new = sample_categorical_logprob(logp, gen);

 double b_new = 0.0;
 const double gamma_new = gamma(static_cast<arma::uword>(k_new));

 if (gamma_new > 0.0) {
  std::normal_distribution<double> norm01(0.0, 1.0);
  const double vbk = std::max(vb_t * gamma_new, 1e-300);
  const double lhs = wi + vei_safe / vbk;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  r(iu) -= wi * diff;

  op.apply_offdiag(i, diff, r);
 }

 b(iu) = b_new;
 comp(iu) = k_new;
}

template <class OpT>
inline void sampleBeta_SBayesRC_ST_csr_scaled(
  int i,
  const arma::rowvec& pi_i,
  const arma::vec& gamma,
  double vb_t,
  const arma::rowvec& prior_scale,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const OpT& op,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double wi = ww(iu);

 if (!std::isfinite(wi) || wi <= 0.0) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr_scaled: invalid ww value.");
 }

 if (!std::isfinite(vb_t) || vb_t <= 0.0) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr_scaled: invalid vb.");
 }

 const double scale_i = prior_scale(iu);
 if (!std::isfinite(scale_i) || scale_i <= 0.0) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr_scaled: invalid prior_scale value.");
 }

 const int Kgamma = static_cast<int>(gamma.n_elem);
 if (static_cast<int>(pi_i.n_elem) != Kgamma) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr_scaled: pi_i/gamma length mismatch.");
 }

 const double vei_safe = std::max(vei_i, 1e-300);
 const double score = r(iu) + wi * b(iu);

 std::vector<double> logp(static_cast<std::size_t>(Kgamma));

 for (int k = 0; k < Kgamma; ++k) {
  const double pik = std::max(static_cast<double>(pi_i(static_cast<arma::uword>(k))), 1e-300);
  const double gk = gamma(static_cast<arma::uword>(k));

  if (!std::isfinite(gk) || gk < 0.0) {
   throw std::runtime_error("sampleBeta_SBayesRC_ST_csr_scaled: gamma must be non-negative.");
  }

  if (gk <= 0.0) {
   logp[static_cast<std::size_t>(k)] = std::log(pik);
  } else {
   const double vbk = std::max(vb_t * gk * scale_i, 1e-300);
   const double denom = std::max(vei_safe + wi * vbk, 1e-300);

   const double logBF =
    0.5 * std::log(vei_safe / denom)
    + 0.5 * score * score * vbk / (vei_safe * denom);

   logp[static_cast<std::size_t>(k)] = std::log(pik) + logBF;
  }
 }

 const int k_new = sample_categorical_logprob(logp, gen);

 double b_new = 0.0;
 const double gamma_new = gamma(static_cast<arma::uword>(k_new));

 if (gamma_new > 0.0) {
  std::normal_distribution<double> norm01(0.0, 1.0);
  const double vbk = std::max(vb_t * gamma_new * scale_i, 1e-300);
  const double lhs = wi + vei_safe / vbk;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  r(iu) -= wi * diff;

  op.apply_offdiag(i, diff, r);
 }

 b(iu) = b_new;
 comp(iu) = k_new;
}

struct SBayesRCLDLDFriends {
 std::vector<uint64_t> ptr;
 std::vector<int> idx;
 std::vector<double> r2;
};

template <class OpT>
struct SBayesRCOperatorContext {
 OpT op;
 SBayesRCLDLDFriends ld_swap_friends;
 Rcpp::List diagnostics;

 SBayesRCOperatorContext(
   OpT&& op_,
   SBayesRCLDLDFriends&& friends_,
   const Rcpp::List& diagnostics_
 ) :
  op(std::move(op_)),
  ld_swap_friends(std::move(friends_)),
  diagnostics(diagnostics_) {}
};

inline SBayesRCLDLDFriends build_ld_swap_friends_sbayesrc_ST_csr(
  int m,
  const STLDCSR& ld,
  const std::vector<double>& xx,
  double min_r2,
  int max_friends
) {
 std::vector<std::vector<std::pair<int, double>>> rows(static_cast<std::size_t>(m));

 for (int i = 0; i < m; ++i) {
  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end = ld.ptr[static_cast<std::size_t>(i + 1)];

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

 SBayesRCLDLDFriends friends;
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

template <class OpT>
inline void set_marker_state_sbayesrc_ST_csr(
  int i,
  double b_new,
  int comp_new,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const OpT& op
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  r(iu) -= ww(iu) * diff;

  op.apply_offdiag(i, diff, r);
 }

 b(iu) = b_new;
 comp(iu) = comp_new;
}

inline int count_null_ld_friends_sbayesrc_ST_csr(
  int i,
  const arma::Row<int>& comp,
  const arma::rowvec& b,
  const SBayesRCLDLDFriends& friends
) {
 int n = 0;
 const uint64_t start = friends.ptr[static_cast<std::size_t>(i)];
 const uint64_t end = friends.ptr[static_cast<std::size_t>(i + 1)];

 for (uint64_t p = start; p < end; ++p) {
  const int j = friends.idx[static_cast<std::size_t>(p)];
  const arma::uword ju = static_cast<arma::uword>(j);
  if (comp(ju) == 0 && b(ju) == 0.0) ++n;
 }

 return n;
}

inline int collect_ld_swap_candidates_sbayesrc_ST_csr(
  int m,
  const arma::Row<int>& comp,
  const arma::rowvec& b,
  const SBayesRCLDLDFriends& friends,
  std::vector<int>& candidates,
  std::vector<int>& n_null
) {
 candidates.clear();
 n_null.clear();

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (comp(iu) <= 0 || b(iu) == 0.0) continue;

  const int nf = count_null_ld_friends_sbayesrc_ST_csr(i, comp, b, friends);
  if (nf > 0) {
   candidates.push_back(i);
   n_null.push_back(nf);
  }
 }

 return static_cast<int>(candidates.size());
}

inline double residual_sse_sbayesrc_ST_csr(
  int m,
  const arma::rowvec& b,
  const arma::rowvec& wy,
  const arma::rowvec& r,
  double yy
) {
 double bwy = 0.0;
 double br = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  bwy += b(iu) * wy(iu);
  br += b(iu) * r(iu);
 }

 return yy - bwy - br;
}

inline double clamp_sbayesrc_prob(double p) {
 if (!std::isfinite(p)) return 1e-300;
 return std::min(std::max(p, 1e-300), 1.0 - 1e-12);
}

inline double log_component_prior_ratio_sbayesrc_ST_csr(
  int source,
  int target,
  int component_moved,
  const arma::mat& snpPi
) {
 const arma::uword su = static_cast<arma::uword>(source);
 const arma::uword tu = static_cast<arma::uword>(target);
 const arma::uword cu = static_cast<arma::uword>(component_moved);

 return
  std::log(clamp_sbayesrc_prob(snpPi(tu, cu))) +
  std::log(clamp_sbayesrc_prob(snpPi(su, 0))) -
  std::log(clamp_sbayesrc_prob(snpPi(su, cu))) -
  std::log(clamp_sbayesrc_prob(snpPi(tu, 0)));
}

inline void throw_ld_swap_error_sbayesrc_ST_csr(
  int trait,
  int chain,
  int iter,
  int source,
  int target,
  int component_moved,
  double sse_old,
  double sse_new,
  double vei,
  double log_prior_ratio,
  double log_q_forward,
  double log_q_reverse,
  const arma::rowvec& b,
  const arma::rowvec& r
) {
 throw std::runtime_error(
  "SBayesRC CSR LD-swap invalid proposal. trait=" +
  std::to_string(trait) +
  ", chain=" + std::to_string(chain) +
  ", iter=" + std::to_string(iter) +
  ", source=" + std::to_string(source) +
  ", target=" + std::to_string(target) +
  ", component_moved=" + std::to_string(component_moved) +
  ", sse_old=" + std::to_string(sse_old) +
  ", sse_new=" + std::to_string(sse_new) +
  ", vei=" + std::to_string(vei) +
  ", log_prior_ratio=" + std::to_string(log_prior_ratio) +
  ", log_q_forward=" + std::to_string(log_q_forward) +
  ", log_q_reverse=" + std::to_string(log_q_reverse) +
  ", b_finite=" + std::to_string(b.is_finite() ? 1 : 0) +
  ", r_finite=" + std::to_string(r.is_finite() ? 1 : 0)
 );
}

template <class OpT>
inline bool attempt_ld_swap_sbayesrc_ST_csr(
  int m,
  int trait,
  int chain,
  int iter,
  double vei,
  double yy,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  const arma::mat& snpPi,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const OpT& op,
  const SBayesRCLDLDFriends& friends,
  std::mt19937& gen,
  bool& attempted
) {
 attempted = false;
 if (!std::isfinite(vei) || vei <= 0.0) return false;

 std::vector<int> candidates;
 std::vector<int> n_null;
 const int n_candidates =
  collect_ld_swap_candidates_sbayesrc_ST_csr(m, comp, b, friends, candidates, n_null);
 if (n_candidates <= 0) return false;

 std::uniform_int_distribution<int> pick_candidate(0, n_candidates - 1);
 const int cand_pos = pick_candidate(gen);
 const int j = candidates[static_cast<std::size_t>(cand_pos)];
 const int n_forward_friends = n_null[static_cast<std::size_t>(cand_pos)];
 if (n_forward_friends <= 0) return false;

 std::uniform_int_distribution<int> pick_friend(0, n_forward_friends - 1);
 const int friend_pos = pick_friend(gen);

 int k = -1;
 int seen = 0;
 const uint64_t start = friends.ptr[static_cast<std::size_t>(j)];
 const uint64_t end = friends.ptr[static_cast<std::size_t>(j + 1)];

 for (uint64_t p = start; p < end; ++p) {
  const int jj = friends.idx[static_cast<std::size_t>(p)];
  const arma::uword jju = static_cast<arma::uword>(jj);
  if (comp(jju) != 0 || b(jju) != 0.0) continue;

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
 const int comp_j_old = comp(ju);
 const int comp_k_old = comp(ku);
 if (comp_j_old <= 0 || b_j_old == 0.0 || comp_k_old != 0 || b_k_old != 0.0) return false;

 attempted = true;
 const double sse_old = residual_sse_sbayesrc_ST_csr(m, b, wy, r, yy);
 if (!std::isfinite(sse_old)) {
  throw_ld_swap_error_sbayesrc_ST_csr(
   trait, chain, iter, j, k, comp_j_old, sse_old,
   std::numeric_limits<double>::quiet_NaN(), vei,
   std::numeric_limits<double>::quiet_NaN(),
   std::numeric_limits<double>::quiet_NaN(),
   std::numeric_limits<double>::quiet_NaN(), b, r
  );
 }

 const arma::rowvec r_old = r;
 set_marker_state_sbayesrc_ST_csr(j, 0.0, 0, ww, r, b, comp, op);
 set_marker_state_sbayesrc_ST_csr(k, b_j_old, comp_j_old, ww, r, b, comp, op);

 const double sse_new = residual_sse_sbayesrc_ST_csr(m, b, wy, r, yy);

 std::vector<int> reverse_candidates;
 std::vector<int> reverse_n_null;
 const int n_reverse_candidates =
  collect_ld_swap_candidates_sbayesrc_ST_csr(
   m, comp, b, friends, reverse_candidates, reverse_n_null
  );

 int n_reverse_friends = 0;
 for (std::size_t pos = 0; pos < reverse_candidates.size(); ++pos) {
  if (reverse_candidates[pos] == k) {
   n_reverse_friends = reverse_n_null[pos];
   break;
  }
 }

 double log_q_forward = std::numeric_limits<double>::quiet_NaN();
 double log_q_reverse = std::numeric_limits<double>::quiet_NaN();
 if (n_candidates > 0 && n_forward_friends > 0) {
  log_q_forward =
   -std::log(static_cast<double>(n_candidates)) -
   std::log(static_cast<double>(n_forward_friends));
 }
 if (n_reverse_candidates > 0 && n_reverse_friends > 0) {
  log_q_reverse =
   -std::log(static_cast<double>(n_reverse_candidates)) -
   std::log(static_cast<double>(n_reverse_friends));
 }

 const double log_prior_ratio =
  log_component_prior_ratio_sbayesrc_ST_csr(j, k, comp_j_old, snpPi);

 if (!std::isfinite(sse_new) ||
     !std::isfinite(log_prior_ratio) ||
     !std::isfinite(log_q_forward) ||
     !std::isfinite(log_q_reverse)) {
  r = r_old;
  b(ju) = b_j_old;
  b(ku) = b_k_old;
  comp(ju) = comp_j_old;
  comp(ku) = comp_k_old;
  throw_ld_swap_error_sbayesrc_ST_csr(
   trait, chain, iter, j, k, comp_j_old, sse_old, sse_new,
   vei, log_prior_ratio, log_q_forward, log_q_reverse, b, r
  );
 }

 const double log_alpha =
  -0.5 * (sse_new - sse_old) / vei +
  log_prior_ratio +
  log_q_reverse - log_q_forward;

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 const bool accept = std::log(std::max(runif(gen), 1e-300)) < log_alpha;

 if (!accept) {
  r = r_old;
  b(ju) = b_j_old;
  b(ku) = b_k_old;
  comp(ju) = comp_j_old;
  comp(ku) = comp_k_old;
 }

 return accept;
}

template <class OpT>
inline bool attempt_ld_swap_sbayesrc_ST_csr_scaled(
  int m,
  int trait,
  int chain,
  int iter,
  double vei,
  double vb,
  double yy,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  const arma::vec& gamma,
  const arma::mat& snpPi,
  const arma::rowvec& prior_scale,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const OpT& op,
  const SBayesRCLDLDFriends& friends,
  std::mt19937& gen,
  bool& attempted
) {
 attempted = false;
 if (!std::isfinite(vei) || vei <= 0.0) return false;

 std::vector<int> candidates;
 std::vector<int> n_null;
 const int n_candidates =
  collect_ld_swap_candidates_sbayesrc_ST_csr(m, comp, b, friends, candidates, n_null);
 if (n_candidates <= 0) return false;

 std::uniform_int_distribution<int> pick_candidate(0, n_candidates - 1);
 const int cand_pos = pick_candidate(gen);
 const int j = candidates[static_cast<std::size_t>(cand_pos)];
 const int n_forward_friends = n_null[static_cast<std::size_t>(cand_pos)];
 if (n_forward_friends <= 0) return false;

 std::uniform_int_distribution<int> pick_friend(0, n_forward_friends - 1);
 const int friend_pos = pick_friend(gen);

 int k = -1;
 int seen = 0;
 const uint64_t start = friends.ptr[static_cast<std::size_t>(j)];
 const uint64_t end = friends.ptr[static_cast<std::size_t>(j + 1)];

 for (uint64_t p = start; p < end; ++p) {
  const int jj = friends.idx[static_cast<std::size_t>(p)];
  const arma::uword jju = static_cast<arma::uword>(jj);
  if (comp(jju) != 0 || b(jju) != 0.0) continue;

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
 const int comp_j_old = comp(ju);
 const int comp_k_old = comp(ku);
 if (comp_j_old <= 0 || b_j_old == 0.0 || comp_k_old != 0 || b_k_old != 0.0) return false;
 if (comp_j_old >= static_cast<int>(gamma.n_elem)) {
  throw std::runtime_error("attempt_ld_swap_sbayesrc_ST_csr_scaled: component index out of range.");
 }
 const double gk = gamma(static_cast<arma::uword>(comp_j_old));
 if (!std::isfinite(gk) || gk <= 0.0) {
  throw std::runtime_error("attempt_ld_swap_sbayesrc_ST_csr_scaled: active component has invalid gamma.");
 }
 const double vb_j = std::max(vb * gk * prior_scale(ju), 1e-300);
 const double vb_k = std::max(vb * gk * prior_scale(ku), 1e-300);

 attempted = true;
 const double sse_old = residual_sse_sbayesrc_ST_csr(m, b, wy, r, yy);
 if (!std::isfinite(sse_old)) {
  throw_ld_swap_error_sbayesrc_ST_csr(
   trait, chain, iter, j, k, comp_j_old, sse_old,
   std::numeric_limits<double>::quiet_NaN(), vei,
   std::numeric_limits<double>::quiet_NaN(),
   std::numeric_limits<double>::quiet_NaN(),
   std::numeric_limits<double>::quiet_NaN(), b, r
  );
 }

 const arma::rowvec r_old = r;
 set_marker_state_sbayesrc_ST_csr(j, 0.0, 0, ww, r, b, comp, op);
 set_marker_state_sbayesrc_ST_csr(k, b_j_old, comp_j_old, ww, r, b, comp, op);

 const double sse_new = residual_sse_sbayesrc_ST_csr(m, b, wy, r, yy);

 std::vector<int> reverse_candidates;
 std::vector<int> reverse_n_null;
 const int n_reverse_candidates =
  collect_ld_swap_candidates_sbayesrc_ST_csr(
   m, comp, b, friends, reverse_candidates, reverse_n_null
  );

 int n_reverse_friends = 0;
 for (std::size_t pos = 0; pos < reverse_candidates.size(); ++pos) {
  if (reverse_candidates[pos] == k) {
   n_reverse_friends = reverse_n_null[pos];
   break;
  }
 }

 double log_q_forward = std::numeric_limits<double>::quiet_NaN();
 double log_q_reverse = std::numeric_limits<double>::quiet_NaN();
 if (n_candidates > 0 && n_forward_friends > 0) {
  log_q_forward =
   -std::log(static_cast<double>(n_candidates)) -
   std::log(static_cast<double>(n_forward_friends));
 }
 if (n_reverse_candidates > 0 && n_reverse_friends > 0) {
  log_q_reverse =
   -std::log(static_cast<double>(n_reverse_candidates)) -
   std::log(static_cast<double>(n_reverse_friends));
 }

 const double log_component_prior_ratio =
  log_component_prior_ratio_sbayesrc_ST_csr(j, k, comp_j_old, snpPi);
 // Prior terms cancel only when marker-specific prior scales are equal.
 const double log_effect_prior_ratio =
  -0.5 * (std::log(vb_k / vb_j) +
          b_j_old * b_j_old * (1.0 / vb_k - 1.0 / vb_j));
 const double log_prior_ratio =
  log_component_prior_ratio + log_effect_prior_ratio;

 if (!std::isfinite(sse_new) ||
     !std::isfinite(log_prior_ratio) ||
     !std::isfinite(log_q_forward) ||
     !std::isfinite(log_q_reverse)) {
  r = r_old;
  b(ju) = b_j_old;
  b(ku) = b_k_old;
  comp(ju) = comp_j_old;
  comp(ku) = comp_k_old;
  throw_ld_swap_error_sbayesrc_ST_csr(
   trait, chain, iter, j, k, comp_j_old, sse_old, sse_new,
   vei, log_prior_ratio, log_q_forward, log_q_reverse, b, r
  );
 }

 const double log_alpha =
  -0.5 * (sse_new - sse_old) / vei +
  log_prior_ratio +
  log_q_reverse - log_q_forward;

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 const bool accept = std::log(std::max(runif(gen), 1e-300)) < log_alpha;

 if (!accept) {
  r = r_old;
  b(ju) = b_j_old;
  b(ku) = b_k_old;
  comp(ju) = comp_j_old;
  comp(ku) = comp_k_old;
 }

 return accept;
}

inline void sampleB_SBayesRC_ST_csr(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  const arma::vec& gamma,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb_scaled = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const int k = comp(iu);

  if (k > 0) {
   if (k >= static_cast<int>(gamma.n_elem)) {
    throw std::runtime_error("sampleB_SBayesRC_ST_csr: component index out of range.");
   }

   const double gk = gamma(static_cast<arma::uword>(k));

   if (!std::isfinite(gk) || gk <= 0.0) {
    throw std::runtime_error("sampleB_SBayesRC_ST_csr: active component has invalid gamma.");
   }

   ssb_scaled += b(iu) * b(iu) / gk;
   dfb += 1.0;
  }
 }

 const double scale = ssb_scaled + nub * ssb_prior;

 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleB_SBayesRC_ST_csr: invalid scale.");
 }

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);

 vb = std::max(scale / chi2, 1e-12);
}

inline void sampleB_SBayesRC_ST_csr_scaled(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  const arma::vec& gamma,
  const arma::rowvec& prior_scale,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb_scaled = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const int k = comp(iu);

  if (k > 0) {
   if (k >= static_cast<int>(gamma.n_elem)) {
    throw std::runtime_error("sampleB_SBayesRC_ST_csr_scaled: component index out of range.");
   }

   const double gk = gamma(static_cast<arma::uword>(k));

   if (!std::isfinite(gk) || gk <= 0.0) {
    throw std::runtime_error("sampleB_SBayesRC_ST_csr_scaled: active component has invalid gamma.");
   }

   const double scale_i = prior_scale(iu);
   if (!std::isfinite(scale_i) || scale_i <= 0.0) {
    throw std::runtime_error("sampleB_SBayesRC_ST_csr_scaled: invalid prior_scale value.");
   }

   // vb is the global variance; fixed selection_s enters the sufficient
   // statistic as b_j^2 / (gamma_m * prior_scale_j).
   ssb_scaled += b(iu) * b(iu) / (gk * scale_i);
   dfb += 1.0;
  }
 }

 const double scale = ssb_scaled + nub * ssb_prior;

 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleB_SBayesRC_ST_csr_scaled: invalid scale.");
 }

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);

 vb = std::max(scale / chi2, 1e-12);
}

inline void fill_selection_s_prior_scale_sbayesrc(
  int m,
  double selection_s,
  const arma::rowvec& log_h,
  arma::rowvec& prior_scale
) {
 if (static_cast<int>(prior_scale.n_elem) != m) {
  prior_scale.set_size(static_cast<arma::uword>(m));
 }

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double scale_i = std::exp((selection_s + 1.0) * log_h(iu));
  if (!std::isfinite(scale_i) || scale_i <= 0.0) {
   throw std::runtime_error("fill_selection_s_prior_scale_sbayesrc: invalid dynamic prior scale.");
  }
  prior_scale(iu) = scale_i;
 }
}

inline double logpost_selection_s_sbayesrc(
  double selection_s,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  double vb,
  const arma::vec& gamma,
  const arma::rowvec& log_h,
  double prior_lower,
  double prior_upper
) {
 if (!std::isfinite(selection_s) ||
     selection_s < prior_lower ||
     selection_s > prior_upper) {
  return -std::numeric_limits<double>::infinity();
 }
 if (!std::isfinite(vb) || vb <= 0.0) {
  return -std::numeric_limits<double>::infinity();
 }

 const int m = static_cast<int>(b.n_elem);
 double lp = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const int k = comp(iu);
  if (k <= 0) continue;
  if (k >= static_cast<int>(gamma.n_elem)) {
   return -std::numeric_limits<double>::infinity();
  }

  const double gk = gamma(static_cast<arma::uword>(k));
  if (!std::isfinite(gk) || gk <= 0.0) {
   return -std::numeric_limits<double>::infinity();
  }

  const double log_q = (selection_s + 1.0) * log_h(iu);
  const double q = std::exp(log_q);
  if (!std::isfinite(q) || q <= 0.0) {
   return -std::numeric_limits<double>::infinity();
  }

  const double bi = b(iu);
  lp += -0.5 * (log_q + bi * bi / (vb * gk * q));
 }

 return lp;
}

inline bool update_selection_s_sbayesrc(
  double& selection_s,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  double vb,
  const arma::vec& gamma,
  const arma::rowvec& log_h,
  double prior_lower,
  double prior_upper,
  double proposal_sd,
  std::mt19937& gen
) {
 std::normal_distribution<double> norm(0.0, proposal_sd);
 const double prop = selection_s + norm(gen);
 if (!std::isfinite(prop) || prop < prior_lower || prop > prior_upper) {
  return false;
 }

 const double lp_current = logpost_selection_s_sbayesrc(
  selection_s, b, comp, vb, gamma, log_h, prior_lower, prior_upper
 );
 const double lp_prop = logpost_selection_s_sbayesrc(
  prop, b, comp, vb, gamma, log_h, prior_lower, prior_upper
 );
 const double log_alpha = lp_prop - lp_current;
 if (!std::isfinite(log_alpha)) return false;

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 if (std::log(std::max(runif(gen), 1e-300)) < log_alpha) {
  selection_s = prop;
  return true;
 }

 return false;
}

inline double computeLE_SBayesRC_ST_csr(
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

#include "blr_csr_sbayesrc_core_impl.h"

struct CsrSBayesRCBindingMetadata {
 const arma::mat& wy_mat;
 const arma::vec& gamma;
 const std::vector<int>& sample_size;
 const Rcpp::List& operator_diagnostics;
 int marker_count;
 int trait_count;
 int annotation_count;
 int component_count;
 int step_count;
 int iterations;
 int burnin;
 int thinning;
 int chains;
 bool keep_chains;
 bool update_ld_swap;
 bool estimate_selection_s;
 bool use_selection_s_prior_scale;
};

static Rcpp::List stblr_csr_sbayesrc_result_to_raw(
 const sblr::core::CsrSBayesRCExecutionResult& execution_result,
 const CsrSBayesRCBindingMetadata& metadata
);

template <class MakeOperator>
Rcpp::List stblr_cpg_omp_csr_sbayesrc_impl(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> comp_init,
  bool use_comp_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  std::string ld_prefix,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  arma::mat A,
  arma::vec gamma,
  arma::mat alpha_init,
  arma::vec sigmaSqAlpha_init,
  bool intercept_flat,
  double sigmaSqAlpha_a,
  double sigmaSqAlpha_b,
  double pi_floor,
  double nub,
  double nue,
  bool updateAlpha,
  bool updateB,
  bool updateE,
  int alpha_update_every,
  double adjE,
  std::vector<int> n,
  int nit,
  int nburn,
  int nthin,
  int ncores,
  int seed,
  int nchains,
  bool keep_chains,
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds,
  bool updateLDswap,
  double ld_swap_prob,
  double ld_swap_r2,
  int ld_swap_max_friends,
  int ld_swap_moves,
  Rcpp::Nullable<Rcpp::NumericVector> selection_s_prior_scale,
  bool estimate_selection_s,
  double selection_s_init,
  Rcpp::NumericVector selection_s_prior,
  double selection_s_proposal_sd,
  Rcpp::Nullable<Rcpp::NumericVector> selection_s_log_h,
  MakeOperator make_operator
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: nt must be positive.");
 }

 const int m = static_cast<int>(wy[0].size());

 if (m <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: m must be positive.");
 }

 if (nit <= 0 || nburn < 0 || nthin <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: invalid nit/nburn/nthin.");
 }

 if (alpha_update_every <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_update_every must be positive.");
 }

 if (nchains <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: nchains must be positive.");
 }

 std::vector<int> chain_seeds_vec;
 if (chain_seeds.isNotNull()) {
  Rcpp::IntegerVector chain_seeds_r(chain_seeds);
  chain_seeds_vec = Rcpp::as<std::vector<int>>(chain_seeds_r);
 }

 if (!chain_seeds_vec.empty() && static_cast<int>(chain_seeds_vec.size()) != nchains) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: chain_seeds must have length nchains.");
 }

 if (!std::isfinite(ld_swap_prob) || ld_swap_prob < 0.0 || ld_swap_prob > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: ld_swap_prob must be in [0, 1].");
 }
 if (!std::isfinite(ld_swap_r2) || ld_swap_r2 < 0.0 || ld_swap_r2 > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: ld_swap_r2 must be in [0, 1].");
 }
 if (ld_swap_max_friends <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: ld_swap_max_friends must be positive.");
 }
 if (ld_swap_moves < 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: ld_swap_moves must be non-negative.");
 }

 bool use_selection_s_prior_scale = selection_s_prior_scale.isNotNull();
 Rcpp::NumericVector selection_s_prior_scale_vec;
 if (use_selection_s_prior_scale) {
  selection_s_prior_scale_vec = Rcpp::NumericVector(selection_s_prior_scale);
  use_selection_s_prior_scale = selection_s_prior_scale_vec.size() > 0;
 }

 if (use_selection_s_prior_scale &&
     static_cast<int>(selection_s_prior_scale_vec.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_prior_scale must have length m.");
 }

 if (estimate_selection_s && use_selection_s_prior_scale) {
  throw std::runtime_error(
   "stblr_cpg_omp_csr_sbayesrc: fixed selection_s_prior_scale and estimate_selection_s cannot both be used."
  );
 }

 if (estimate_selection_s) {
  if (!std::isfinite(selection_s_init)) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_init must be finite.");
  }
  if (selection_s_prior.size() != 2 ||
      !std::isfinite(selection_s_prior[0]) ||
      !std::isfinite(selection_s_prior[1])) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_prior must be finite length 2.");
  }
  if (selection_s_prior[0] >= selection_s_prior[1]) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_prior lower bound must be less than upper bound.");
  }
  if (selection_s_init < selection_s_prior[0] ||
      selection_s_init > selection_s_prior[1]) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_init must lie within selection_s_prior.");
  }
  if (!std::isfinite(selection_s_proposal_sd) || selection_s_proposal_sd <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_proposal_sd must be positive finite.");
  }
  if (selection_s_log_h.isNull()) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_log_h is required when estimate_selection_s = true.");
  }
 }

 if ((int)ww.size() != nt || (int)b_init.size() != nt ||
     (int)yy.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: inconsistent trait dimensions.");
 }

 if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: priors must be nt x nt.");
 }

 if ((int)A.n_rows != m) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: A must have m rows.");
 }

 const int nAnno = static_cast<int>(A.n_cols);
 if (nAnno <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: A must have at least one column.");
 }

 const int Kgamma = static_cast<int>(gamma.n_elem);
 if (Kgamma < 2) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma must have at least two components, including zero.");
 }

 if (!std::isfinite(gamma(0)) || gamma(0) != 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma[0] must be exactly 0.0.");
 }

 for (int k = 1; k < Kgamma; ++k) {
  const double g = gamma(static_cast<arma::uword>(k));
  if (!std::isfinite(g) || g <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma[k] must be positive for k > 0.");
  }
 }

 const int nstep = Kgamma - 1;

 if ((int)alpha_init.n_rows != nAnno || (int)alpha_init.n_cols != nstep) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_init must be ncol(A) x (length(gamma)-1).");
 }

 if ((int)sigmaSqAlpha_init.n_elem != nstep) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: sigmaSqAlpha_init must have length length(gamma)-1.");
 }

 if (!std::isfinite(pi_floor) || pi_floor <= 0.0 || pi_floor >= 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: pi_floor must be in (0,1).");
 }

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m ||
      (int)ww[t].size() != m ||
      (int)b_init[t].size() != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: inconsistent marker dimensions.");
  }
 }

 if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: B must be nt x nt.");
 }

 if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: E must be nt x nt.");
 }

 if (use_r_init) {
  if (static_cast<int>(r_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: r_init must have length nt when use_r_init = true.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(r_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: each r_init[t] must have length m.");
   }
  }
 }

 if (use_comp_init) {
  if (static_cast<int>(comp_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: comp_init must have length nt when enabled.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(comp_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: comp_init[t] must have length m.");
   }
  }
 }

 arma::rowvec prior_scale;
 if (use_selection_s_prior_scale) {
  prior_scale.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double scale_i = selection_s_prior_scale_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(scale_i) || scale_i <= 0.0) {
    throw std::runtime_error(
     "stblr_cpg_omp_csr_sbayesrc: selection_s_prior_scale must contain positive finite values."
    );
   }
   prior_scale(static_cast<arma::uword>(i)) = scale_i;
  }
 }

 arma::rowvec selection_s_log_h_row;
 if (estimate_selection_s) {
  Rcpp::NumericVector selection_s_log_h_vec(selection_s_log_h);
  if (static_cast<int>(selection_s_log_h_vec.size()) != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_log_h must have length m.");
  }
  selection_s_log_h_row.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double log_h_i = selection_s_log_h_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(log_h_i)) {
    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: selection_s_log_h must contain finite values.");
   }
   selection_s_log_h_row(static_cast<arma::uword>(i)) = log_h_i;
  }
 }

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);
 arma::mat comp_mat(nt, m, arma::fill::zeros);

 arma::vec yy_vec(nt, arma::fill::zeros);
 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  yy_vec(static_cast<arma::uword>(t)) = yy[t];

  for (int i = 0; i < m; ++i) {
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword iu = static_cast<arma::uword>(i);

   wy_mat(tu, iu) = wy[t][i];
   ww_mat(tu, iu) = ww[t][i];
   b_mat(tu, iu)  = b_init[t][i];
  }

  if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: priors must be nt x nt.");
  }

  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
   sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
  }
 }

 for (int t = 1; t < nt; ++t) {
  if (n[t] != n[0]) {
   throw std::runtime_error(
     "stblr_cpg_omp_csr_sbayesrc: current shared-LD scaling assumes equal n across traits."
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
      "stblr_cpg_omp_csr_sbayesrc: ww contains invalid value before LD pre-scaling."
    );
   }

   if (std::abs(w0 - wt) > tol) {
    throw std::runtime_error(
      "stblr_cpg_omp_csr_sbayesrc: ww differs across traits; pre-scaled shared ST LD is invalid."
    );
   }
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);

 for (int i = 0; i < m; ++i) {
  const double wi = ww_mat(0, static_cast<arma::uword>(i));
  if (!std::isfinite(wi) || wi <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: ww contains invalid value in trait 0.");
  }
  xx[static_cast<std::size_t>(i)] = wi;
 }

 arma::rowvec xx_row(static_cast<arma::uword>(m));
 for (int i = 0; i < m; ++i) {
  xx_row(static_cast<arma::uword>(i)) = xx[static_cast<std::size_t>(i)];
 }
 auto operator_context = make_operator(
  m, xx, xx_row, wy_mat, updateLDswap, ld_swap_r2, ld_swap_max_friends
 );
 const auto& op = operator_context.op;
 const SBayesRCLDLDFriends& ld_swap_friends = operator_context.ld_swap_friends;

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
           [&](int a, int b) {
            return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
           });

 sblr::core::SBayesRCAnnotationDesignView annotation_contract;
 annotation_contract.marker_count = static_cast<std::size_t>(m);
 annotation_contract.annotation_count = static_cast<std::size_t>(nAnno);
 annotation_contract.values = A.memptr();
 annotation_contract.value_count = A.n_elem;
 annotation_contract.includes_intercept = true;

 sblr::core::SBayesRCComponentSpec component_contract;
 component_contract.scales.assign(gamma.begin(), gamma.end());

 sblr::core::SBayesRCAlphaSpec alpha_contract;
 alpha_contract.annotation_count = static_cast<std::size_t>(nAnno);
 alpha_contract.step_count = static_cast<std::size_t>(nstep);
 alpha_contract.intercept_flat = intercept_flat;
 alpha_contract.update = updateAlpha;
 alpha_contract.variance_prior_a = sigmaSqAlpha_a;
 alpha_contract.variance_prior_b = sigmaSqAlpha_b;
 alpha_contract.update_every = alpha_update_every;

 sblr::core::SBayesRCProbabilityPolicy probability_contract;
 probability_contract.probability_floor = pi_floor;
 probability_contract.stick_order.resize(static_cast<std::size_t>(nstep));
 for (int j = 0; j < nstep; ++j) {
  probability_contract.stick_order[static_cast<std::size_t>(j)] =
   static_cast<std::size_t>(j + 1);
 }

 sblr::core::SBayesRCPriors prior_contract;
 prior_contract.marker_df = nub;
 prior_contract.residual_df = nue;

 sblr::core::SBayesRCControls control_contract;
 control_contract.iterations = nit;
 control_contract.burnin = nburn;
 control_contract.thinning = nthin;
 control_contract.chains = nchains;
 control_contract.cores = ncores;
 control_contract.seed = seed;
 control_contract.chain_seeds = chain_seeds_vec;
 control_contract.keep_chains = keep_chains;
 control_contract.update_marker_variance = updateB;
 control_contract.update_residual_variance = updateE;
 control_contract.rebuild_residual_before_update = rebuild_r_before_updateE;
 control_contract.update_ld_swap = updateLDswap;
 control_contract.ld_swap_probability = ld_swap_prob;
 control_contract.ld_swap_r2 = ld_swap_r2;
 control_contract.ld_swap_max_friends = ld_swap_max_friends;
 control_contract.ld_swap_moves = ld_swap_moves;

 sblr::core::SBayesRCOutputControl output_contract;
 output_contract.keep_chains = keep_chains;

 using Operator = typename std::remove_reference<decltype(op)>::type;
 sblr::core::CsrSBayesRCExecutionContext<Operator> execution_context{
  op,
  ld_swap_friends,
  wy_mat,
  ww_mat,
  b_mat,
  r_mat,
  comp_mat,
  yy_vec,
  ssb_prior_mat,
  sse_prior_mat,
  B,
  E,
  A,
  gamma,
  alpha_init,
  sigmaSqAlpha_init,
  comp_init,
  r_init,
  n,
  chain_seeds_vec,
  order,
  prior_scale,
  selection_s_log_h_row,
  annotation_contract,
  component_contract,
  alpha_contract,
  probability_contract,
  prior_contract,
  control_contract,
  output_contract,
  m,
  nt,
  nAnno,
  Kgamma,
  nstep,
  nit,
  nburn,
  nthin,
  ncores,
  seed,
  nchains,
  use_comp_init,
  use_r_init,
  rebuild_r_before_updateE,
  intercept_flat,
  sigmaSqAlpha_a,
  sigmaSqAlpha_b,
  pi_floor,
  nub,
  nue,
  updateAlpha,
  updateB,
  updateE,
  alpha_update_every,
  adjE,
  updateLDswap,
  ld_swap_prob,
  ld_swap_moves,
  use_selection_s_prior_scale,
  estimate_selection_s,
  selection_s_init,
  selection_s_prior[0],
  selection_s_prior[1],
  selection_s_proposal_sd
 };
 auto execution_result = sblr::core::run_csr_sbayesrc(execution_context);

 CsrSBayesRCBindingMetadata binding_metadata{
  wy_mat,
  gamma,
  n,
  operator_context.diagnostics,
  m,
  nt,
  nAnno,
  Kgamma,
  nstep,
  nit,
  nburn,
  nthin,
  nchains,
  keep_chains,
  updateLDswap,
  estimate_selection_s,
  use_selection_s_prior_scale
 };
 return stblr_csr_sbayesrc_result_to_raw(execution_result, binding_metadata);
}

static Rcpp::List stblr_csr_sbayesrc_result_to_raw(
 const sblr::core::CsrSBayesRCExecutionResult& execution_result,
 const CsrSBayesRCBindingMetadata& metadata
) {
 const arma::mat& wy_mat = metadata.wy_mat;
 const arma::vec& gamma = metadata.gamma;
 const std::vector<int>& n = metadata.sample_size;
 const int m = metadata.marker_count;
 const int nt = metadata.trait_count;
 const int nAnno = metadata.annotation_count;
 const int Kgamma = metadata.component_count;
 const int nstep = metadata.step_count;
 const int nit = metadata.iterations;
 const int nburn = metadata.burnin;
 const int nthin = metadata.thinning;
 const int nchains = metadata.chains;
 const bool keep_chains = metadata.keep_chains;
 const bool updateLDswap = metadata.update_ld_swap;
 const bool estimate_selection_s = metadata.estimate_selection_s;
 const bool use_selection_s_prior_scale = metadata.use_selection_s_prior_scale;

 const auto& bm_task = execution_result.bm_task;
 const auto& dm_task = execution_result.dm_task;
 const auto& comp_task_double = execution_result.comp_task_double;
 const auto& vbs_task = execution_result.vbs_task;
 const auto& vgs_task = execution_result.vgs_task;
 const auto& ves_task = execution_result.ves_task;
 const auto& pis_task = execution_result.pis_task;
 const auto& vles_task = execution_result.vles_task;
 const auto& vlds_task = execution_result.vlds_task;
 const auto& selection_s_task = execution_result.selection_s_task;
 const auto& final_pi_component_task = execution_result.final_pi_component_task;
 const auto& selection_s_attempted_task = execution_result.selection_s_attempted_task;
 const auto& selection_s_accepted_task = execution_result.selection_s_accepted_task;
 const auto& alpha_mean_task = execution_result.alpha_mean_task;
 const auto& sigmaSqAlpha_mean_task = execution_result.sigmaSqAlpha_mean_task;
 const auto& comp_prob_mean_task = execution_result.comp_prob_mean_task;
 const auto& ncomp_mean_task = execution_result.ncomp_mean_task;
 const auto& bm_mat = execution_result.bm_mat;
 const auto& dm_mat = execution_result.dm_mat;
 const auto& bm_sd_mat = execution_result.bm_sd_mat;
 const auto& dm_sd_mat = execution_result.dm_sd_mat;
 const auto& bm_min_mat = execution_result.bm_min_mat;
 const auto& dm_min_mat = execution_result.dm_min_mat;
 const auto& bm_max_mat = execution_result.bm_max_mat;
 const auto& dm_max_mat = execution_result.dm_max_mat;
 const auto& b_result_mat = execution_result.b_mat;
 const auto& r_result_mat = execution_result.r_mat;
 const auto& comp_result_mat = execution_result.comp_mat;
 const auto& vbs_mat = execution_result.vbs_mat;
 const auto& vgs_mat = execution_result.vgs_mat;
 const auto& ves_mat = execution_result.ves_mat;
 const auto& pis_mat = execution_result.pis_mat;
 const auto& vles_mat = execution_result.vles_mat;
 const auto& vlds_mat = execution_result.vlds_mat;
 const auto& selection_s_mat = execution_result.selection_s_mat;
 const auto& final_vb = execution_result.final_vb;
 const auto& final_vg = execution_result.final_vg;
 const auto& final_ve = execution_result.final_ve;
 const auto& final_pi_component = execution_result.final_pi_component;
 const auto& nsamples_vec = execution_result.nsamples_vec;
 const auto& selection_s_attempted_vec = execution_result.selection_s_attempted_vec;
 const auto& selection_s_accepted_vec = execution_result.selection_s_accepted_vec;
 const auto& alpha_mean = execution_result.alpha_mean;
 const auto& sigmaSqAlpha_mean = execution_result.sigmaSqAlpha_mean;
 const auto& comp_prob_mean = execution_result.comp_prob_mean;
 const auto& ncomp_mean = execution_result.ncomp_mean;
 const auto& task_seconds = execution_result.task_seconds;
 const auto& ld_swap_diagnostics = execution_result.ld_swap_diagnostics;
 const auto& ld_swap_chain_diagnostics = execution_result.ld_swap_chain_diagnostics;

 const bool return_chain_summaries = (nchains > 1) || keep_chains;
 const int n_trace = nit + nburn;

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
  for (int t = 0; t < nt; ++t) out(t, t) = x(static_cast<arma::uword>(t));
  return out;
 };

 Rcpp::CharacterVector component_names(Kgamma);
 for (int k = 0; k < Kgamma; ++k) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "gamma_%.2f", gamma(static_cast<arma::uword>(k)));
  component_names[k] = std::string(buf);
 }
 if (Kgamma > 0) component_names[0] = "gamma_0.00";

 Rcpp::CharacterVector step_names(nstep);
 for (int j = 0; j < nstep; ++j) {
  step_names[j] = "step_" + std::to_string(j + 1);
 }

 Rcpp::List comp_prob_out(nt);
 Rcpp::List alpha_mean_out(nt);
 Rcpp::List alpha_final_out(nt);
 arma::mat dm_component_mean_mat(nt, m, arma::fill::zeros);
 Rcpp::NumericMatrix sigma_mean_out(nstep, nt);
 Rcpp::NumericMatrix sigma_final_out(nstep, nt);
 Rcpp::NumericMatrix ncomp_out(nt, Kgamma);
 Rcpp::NumericMatrix pi_mean_out(nt, Kgamma);
 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  Rcpp::NumericMatrix cp = Rcpp::wrap(comp_prob_mean[static_cast<std::size_t>(t)]);
  Rcpp::colnames(cp) = component_names;
  comp_prob_out[t] = cp;
  alpha_mean_out[t] = alpha_mean[static_cast<std::size_t>(t)];
  // The SBayesRC sampler does not retain final alpha draws separately.
  // Expose the posterior mean in both fields to keep the raw namespace complete
  // without adding sampler-state mutations or RNG-dependent behavior.
  alpha_final_out[t] = alpha_mean[static_cast<std::size_t>(t)];
  for (int k = 0; k < Kgamma; ++k) {
   const arma::uword ku = static_cast<arma::uword>(k);
   ncomp_out(t, k) = ncomp_mean[static_cast<std::size_t>(t)](ku);
   pi_mean_out(t, k) = ncomp_mean[static_cast<std::size_t>(t)](ku) /
    static_cast<double>(m);
   for (int i = 0; i < m; ++i) {
    dm_component_mean_mat(tu, static_cast<arma::uword>(i)) +=
     static_cast<double>(k) *
     comp_prob_mean[static_cast<std::size_t>(t)](static_cast<arma::uword>(i), ku);
   }
  }
  for (int j = 0; j < nstep; ++j) {
   const arma::uword ju = static_cast<arma::uword>(j);
   sigma_mean_out(j, t) = sigmaSqAlpha_mean[static_cast<std::size_t>(t)](ju);
   sigma_final_out(j, t) = sigmaSqAlpha_mean[static_cast<std::size_t>(t)](ju);
  }
 }
 Rcpp::colnames(ncomp_out) = component_names;
 Rcpp::colnames(pi_mean_out) = component_names;
 Rcpp::NumericVector selection_mean(nt);
 Rcpp::NumericVector selection_sd(nt);
 Rcpp::NumericVector selection_min(nt);
 Rcpp::NumericVector selection_max(nt);
 Rcpp::NumericVector selection_acceptance(nt);
 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  selection_acceptance[t] =
   selection_s_attempted_vec(tu) > 0.0
   ? selection_s_accepted_vec(tu) / selection_s_attempted_vec(tu)
   : 0.0;
  if (estimate_selection_s) {
   double mean_s = 0.0;
   double min_s = std::numeric_limits<double>::infinity();
   double max_s = -std::numeric_limits<double>::infinity();
   int ns = 0;
   for (int it = nburn; it < n_trace; ++it) {
    const double val = selection_s_mat(tu, static_cast<arma::uword>(it));
    mean_s += val;
    min_s = std::min(min_s, val);
    max_s = std::max(max_s, val);
    ++ns;
   }
   if (ns > 0) mean_s /= static_cast<double>(ns);
   double ss = 0.0;
   if (ns > 1) {
    for (int it = nburn; it < n_trace; ++it) {
     const double diff = selection_s_mat(tu, static_cast<arma::uword>(it)) - mean_s;
     ss += diff * diff;
    }
    selection_sd[t] = std::sqrt(ss / static_cast<double>(ns - 1));
   } else {
    selection_sd[t] = NA_REAL;
   }
   selection_mean[t] = mean_s;
   selection_min[t] = ns > 0 ? min_s : NA_REAL;
   selection_max[t] = ns > 0 ? max_s : NA_REAL;
  } else {
   selection_mean[t] = NA_REAL;
   selection_sd[t] = NA_REAL;
   selection_min[t] = NA_REAL;
   selection_max[t] = NA_REAL;
  }
 }

 Rcpp::NumericVector nsamples(nt);
 Rcpp::IntegerVector n_used(nt);
 Rcpp::NumericVector seconds_mean(nt);
 Rcpp::NumericVector seconds_max(nt);
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
 }

 Rcpp::List marker = Rcpp::List::create(
  Rcpp::Named("bm") = marker_matrix(bm_mat),
  Rcpp::Named("dm") = marker_matrix(dm_mat),
  Rcpp::Named("wy") = marker_matrix(wy_mat),
  Rcpp::Named("r") = marker_matrix(r_result_mat),
  Rcpp::Named("b") = marker_matrix(b_result_mat),
  Rcpp::Named("state") = marker_matrix(comp_result_mat)
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
 Rcpp::Named("ld_swap") = updateLDswap ? Rcpp::wrap(ld_swap_diagnostics) : R_NilValue
 );
 if (metadata.operator_diagnostics.size() > 0) {
  diagnostics["block_eigen"] = metadata.operator_diagnostics;
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
     chain_state[i] = comp_task_double(task_u, iu);
    }
    Rcpp::NumericMatrix chain_ld(1, 3);
    chain_ld(0, 0) = ld_swap_chain_diagnostics(task_u, 2);
    chain_ld(0, 1) = ld_swap_chain_diagnostics(task_u, 3);
    chain_ld(0, 2) = ld_swap_chain_diagnostics(task_u, 4);
    Rcpp::NumericVector chain_pi_mean(Kgamma);
    for (int k = 0; k < Kgamma; ++k) {
     chain_pi_mean[k] =
      ncomp_mean_task[static_cast<std::size_t>(task)](static_cast<arma::uword>(k)) /
      static_cast<double>(m);
    }
    Rcpp::List chain_selection = Rcpp::List::create(
     Rcpp::Named("trace") = R_NilValue,
     Rcpp::Named("acceptance") = R_NilValue
    );
    if (estimate_selection_s) {
     const double attempted = selection_s_attempted_task(task_u);
     const double accepted = selection_s_accepted_task(task_u);
     chain_selection["trace"] = selection_s_task.row(task_u).t();
     chain_selection["acceptance"] =
      attempted > 0.0 ? accepted / attempted : 0.0;
    }
    trait_chains[chain] = Rcpp::List::create(
     Rcpp::Named("marker") = Rcpp::List::create(
      Rcpp::Named("bm") = chain_bm,
      Rcpp::Named("dm") = chain_dm,
      Rcpp::Named("state") = chain_state
     ),
     Rcpp::Named("trace") = Rcpp::List::create(
      Rcpp::Named("vbs") = vbs_task.row(task_u).t(),
      Rcpp::Named("vgs") = vgs_task.row(task_u).t(),
      Rcpp::Named("ves") = ves_task.row(task_u).t(),
      Rcpp::Named("vle") = vles_task.row(task_u).t(),
      Rcpp::Named("vld") = vlds_task.row(task_u).t(),
      Rcpp::Named("pis") = pis_task.row(task_u).t()
     ),
     Rcpp::Named("pi") = Rcpp::List::create(
      Rcpp::Named("final") = final_pi_component_task.row(task_u).t(),
      Rcpp::Named("mean") = chain_pi_mean
     ),
     Rcpp::Named("component") = Rcpp::List::create(
      Rcpp::Named("prob") = comp_prob_mean_task[static_cast<std::size_t>(task)],
      Rcpp::Named("ncomp") = ncomp_mean_task[static_cast<std::size_t>(task)].t()
     ),
     Rcpp::Named("annotation") = Rcpp::List::create(
      Rcpp::Named("alpha") = alpha_mean_task[static_cast<std::size_t>(task)],
      Rcpp::Named("sigmaSqAlpha") = sigmaSqAlpha_mean_task[static_cast<std::size_t>(task)]
     ),
     Rcpp::Named("selection") = chain_selection,
     Rcpp::Named("diagnostics") = Rcpp::List::create(
      Rcpp::Named("ld_swap") = updateLDswap ? Rcpp::wrap(chain_ld) : R_NilValue
     )
    );
   }
   trait_chains.attr("names") = chain_names;
   chains[t] = trait_chains;
  }
  chains.attr("names") = trait_names;
 }

 Rcpp::List selection = Rcpp::List::create(
  Rcpp::Named("enabled") = estimate_selection_s || use_selection_s_prior_scale,
  Rcpp::Named("fixed") = use_selection_s_prior_scale,
  Rcpp::Named("scale") = "standardized_genotype_effect",
  Rcpp::Named("trace") = estimate_selection_s ? Rcpp::wrap(trace_matrix(selection_s_mat)) : R_NilValue,
  Rcpp::Named("mean") = estimate_selection_s ? Rcpp::wrap(selection_mean) : R_NilValue,
  Rcpp::Named("sd") = estimate_selection_s ? Rcpp::wrap(selection_sd) : R_NilValue,
  Rcpp::Named("min") = estimate_selection_s ? Rcpp::wrap(selection_min) : R_NilValue,
  Rcpp::Named("max") = estimate_selection_s ? Rcpp::wrap(selection_max) : R_NilValue,
  Rcpp::Named("acceptance") = estimate_selection_s ? Rcpp::wrap(selection_acceptance) : R_NilValue
 );

 Rcpp::List raw = Rcpp::List::create(
  Rcpp::Named("schema") = Rcpp::List::create(
   Rcpp::Named("class") = "stblr_raw",
   Rcpp::Named("version") = 1
  ),
  Rcpp::Named("meta") = Rcpp::List::create(
   Rcpp::Named("model") = "sbayesrc",
   Rcpp::Named("backend") = "csr_sbayesrc",
   Rcpp::Named("data_level") = "summary",
   Rcpp::Named("prior_type") = "annotation_component",
   Rcpp::Named("m") = m,
   Rcpp::Named("nt") = nt,
   Rcpp::Named("n_trace") = n_trace,
   Rcpp::Named("nit") = nit,
   Rcpp::Named("nburn") = nburn,
   Rcpp::Named("nthin") = nthin,
   Rcpp::Named("nchains") = nchains,
   Rcpp::Named("keep_chains") = keep_chains,
   Rcpp::Named("n_components") = Kgamma,
   Rcpp::Named("n_annotations") = nAnno,
   Rcpp::Named("n_groups") = 0
  ),
  Rcpp::Named("marker") = marker,
  Rcpp::Named("trace") = trace,
  Rcpp::Named("variance") = variance,
  Rcpp::Named("pi") = Rcpp::List::create(
   Rcpp::Named("final") = final_pi_component,
   Rcpp::Named("mean") = pi_mean_out,
   Rcpp::Named("names") = component_names
  ),
  Rcpp::Named("diagnostics") = diagnostics,
  Rcpp::Named("chains") = chains,
  Rcpp::Named("prior") = Rcpp::List::create(),
  Rcpp::Named("group") = Rcpp::List::create(),
  Rcpp::Named("annotation") = Rcpp::List::create(
   Rcpp::Named("annotation_names") = R_NilValue,
   Rcpp::Named("alpha_mean") = alpha_mean_out,
   Rcpp::Named("alpha_final") = alpha_final_out,
   Rcpp::Named("sigmaSqAlpha_mean") = sigma_mean_out,
   Rcpp::Named("sigmaSqAlpha_final") = sigma_final_out
  ),
  Rcpp::Named("component") = Rcpp::List::create(
   Rcpp::Named("names") = component_names,
   Rcpp::Named("gamma") = gamma,
   Rcpp::Named("mixture_var") = gamma,
   Rcpp::Named("prob") = comp_prob_out,
   Rcpp::Named("ncomp") = ncomp_out,
   Rcpp::Named("dm_component_mean") = marker_matrix(dm_component_mean_mat)
  ),
  Rcpp::Named("selection") = selection
 );
 raw.attr("class") = Rcpp::CharacterVector::create("stblr_raw_v1", "stblr_raw", "list");
 return raw;
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_sbayesrc(
  std::vector<std::vector<double>> wy, std::vector<std::vector<double>> ww,
  std::vector<double> yy, std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> comp_init, bool use_comp_init,
  std::vector<std::vector<double>> r_init, bool use_r_init,
  bool rebuild_r_before_updateE, std::string ld_prefix, arma::mat B,
  arma::mat E, std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior, arma::mat A, arma::vec gamma,
  arma::mat alpha_init, arma::vec sigmaSqAlpha_init, bool intercept_flat,
  double sigmaSqAlpha_a, double sigmaSqAlpha_b, double pi_floor, double nub,
  double nue, bool updateAlpha, bool updateB, bool updateE,
  int alpha_update_every, double adjE, std::vector<int> n, int nit,
  int nburn, int nthin, int ncores, int seed, int nchains = 1,
  bool keep_chains = false,
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds = R_NilValue,
  bool updateLDswap = false, double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8, int ld_swap_max_friends = 50,
  int ld_swap_moves = 1,
  Rcpp::Nullable<Rcpp::NumericVector> selection_s_prior_scale = R_NilValue,
  bool estimate_selection_s = false, double selection_s_init = 0.0,
  Rcpp::NumericVector selection_s_prior = Rcpp::NumericVector::create(-3.0, 2.0),
  double selection_s_proposal_sd = 0.35,
  Rcpp::Nullable<Rcpp::NumericVector> selection_s_log_h = R_NilValue
) {
 auto make_csr_operator = [&](int m, const std::vector<double>& xx,
                              const arma::rowvec& xx_row, arma::mat& wy_mat,
                              bool update_ld_swap, double r2, int max_friends) {
  (void)wy_mat;
  STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
  CsrOperator op(ld, xx_row);
  SBayesRCLDLDFriends friends;
  if (update_ld_swap) {
   friends = build_ld_swap_friends_sbayesrc_ST_csr(m, ld, xx, r2, max_friends);
  } else {
   friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  }
  return SBayesRCOperatorContext<CsrOperator>(
   std::move(op), std::move(friends), Rcpp::List::create()
  );
 };
 return stblr_cpg_omp_csr_sbayesrc_impl(
  wy, ww, yy, b_init, comp_init, use_comp_init, r_init, use_r_init,
  rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, A, gamma,
  alpha_init, sigmaSqAlpha_init, intercept_flat, sigmaSqAlpha_a,
  sigmaSqAlpha_b, pi_floor, nub, nue, updateAlpha, updateB, updateE,
  alpha_update_every, adjE, n, nit, nburn, nthin, ncores, seed, nchains,
  keep_chains, chain_seeds, updateLDswap, ld_swap_prob, ld_swap_r2,
  ld_swap_max_friends, ld_swap_moves, selection_s_prior_scale,
  estimate_selection_s, selection_s_init, selection_s_prior,
  selection_s_proposal_sd, selection_s_log_h, make_csr_operator
 );
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_sbayesrc_block_eigen(
  std::vector<std::vector<double>> wy, std::vector<std::vector<double>> ww,
  std::vector<double> yy, std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> comp_init, bool use_comp_init,
  std::vector<std::vector<double>> r_init, bool use_r_init,
  bool rebuild_r_before_updateE, std::string ld_prefix, arma::mat B,
  arma::mat E, std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior, arma::mat A, arma::vec gamma,
  arma::mat alpha_init, arma::vec sigmaSqAlpha_init, bool intercept_flat,
  double sigmaSqAlpha_a, double sigmaSqAlpha_b, double pi_floor, double nub,
  double nue, bool updateAlpha, bool updateB, bool updateE,
  int alpha_update_every, double adjE, std::vector<int> n, int nit,
  int nburn, int nthin, int ncores, int seed, int nchains = 1,
  bool keep_chains = false,
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds = R_NilValue,
  bool updateLDswap = false, double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8, int ld_swap_max_friends = 50,
  int ld_swap_moves = 1,
  Rcpp::Nullable<Rcpp::NumericVector> selection_s_prior_scale = R_NilValue,
  bool estimate_selection_s = false, double selection_s_init = 0.0,
  Rcpp::NumericVector selection_s_prior = Rcpp::NumericVector::create(-3.0, 2.0),
  double selection_s_proposal_sd = 0.35,
  Rcpp::Nullable<Rcpp::NumericVector> selection_s_log_h = R_NilValue,
  Rcpp::CharacterVector bed_files = Rcpp::CharacterVector::create(),
  int n_bed = 0, Rcpp::List cls = R_NilValue,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::NumericVector af = Rcpp::NumericVector::create(),
  Rcpp::IntegerVector block_start = Rcpp::IntegerVector::create(),
  std::string eigen_filter = "hard_truncate", double eigen_tau = 0.01,
  double eigen_eta = 0.0
) {
 if (updateLDswap) {
  throw std::runtime_error(
   "LD-swap is not yet supported with the experimental block-eigen operator."
  );
 }
 if (bed_files.size() <= 0) throw std::runtime_error("bed_files must not be empty.");
 if (n_bed <= 0) throw std::runtime_error("n_bed must be positive.");
 if (cls.size() != bed_files.size()) {
  throw std::runtime_error("cls must have one element per BED file.");
 }
 const std::vector<std::string> bed_cpp =
  stblr_sbayesrc_copy_character_vector(bed_files, "bed_files");
 const std::vector<std::vector<int>> cls_cpp =
  stblr_sbayesrc_copy_int_list(cls, "cls");
 const std::vector<int> rows0 = stblr_sbayesrc_copy_rows0_or_empty(rows, n_bed);
 const std::vector<double> af_cpp = Rcpp::as<std::vector<double>>(af);
 const std::vector<int> starts = Rcpp::as<std::vector<int>>(block_start);
 const EigenFilterMode mode = parse_block_eigen_filter_mode(eigen_filter);

 auto make_block_operator = [&](int m, const std::vector<double>& xx,
                                const arma::rowvec& xx_row, arma::mat& wy_mat,
                                bool update_ld_swap, double r2, int max_friends) {
  (void)xx; (void)xx_row; (void)r2; (void)max_friends;
  if (update_ld_swap) {
   throw std::runtime_error(
    "LD-swap is not yet supported with the experimental block-eigen operator."
   );
  }
  if (static_cast<int>(af_cpp.size()) != m) {
   throw std::runtime_error("af length must equal m for block-eigen SBayesRC.");
  }
  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
   bed_cpp, n_bed, rows0.empty() ? nullptr : rows0.data(),
   static_cast<int>(rows0.size()), cls_cpp
  );
  if (G.m != m) throw std::runtime_error("BED marker count does not match m.");
  std::vector<BlockEigenDiag> block_diag;
  BlockEigenOperator op = build_block_eigen(
   G, af_cpp, starts, mode, eigen_tau, eigen_eta, wy_mat, ncores, &block_diag
  );
  SBayesRCLDLDFriends friends;
  friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  Rcpp::List diag = Rcpp::List::create(
   Rcpp::Named("blocks") = block_eigen_diagnostics_to_data_frame(block_diag)
  );
  return SBayesRCOperatorContext<BlockEigenOperator>(
   std::move(op), std::move(friends), diag
  );
 };
 return stblr_cpg_omp_csr_sbayesrc_impl(
  wy, ww, yy, b_init, comp_init, use_comp_init, r_init, use_r_init,
  rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, A, gamma,
  alpha_init, sigmaSqAlpha_init, intercept_flat, sigmaSqAlpha_a,
  sigmaSqAlpha_b, pi_floor, nub, nue, updateAlpha, updateB, updateE,
  alpha_update_every, adjE, n, nit, nburn, nthin, ncores, seed, nchains,
  keep_chains, chain_seeds, updateLDswap, ld_swap_prob, ld_swap_r2,
  ld_swap_max_friends, ld_swap_moves, selection_s_prior_scale,
  estimate_selection_s, selection_s_init, selection_s_prior,
  selection_s_proposal_sd, selection_s_log_h, make_block_operator
 );
}

// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
//
// #include "st_csr_common.h"
//
// #include <algorithm>
// #include <cmath>
// #include <cstdint>
// #include <cstdlib>
// #include <fstream>
// #include <limits>
// #include <numeric>
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // =============================================================================
// // STBLR summary-stat CSR: SBayesRC-style sampler with overlapping annotations
// // =============================================================================
// //
// // Exported function:
// //   stblr_cpg_omp_csr_sbayesrc(...)
// //
// // This is intentionally a full-sweep CSR sampler. No sparse scheduling is used.
// //
// // Model:
// //   gamma = c(0, 0.01, 0.1, 1)             fixed mixture variance multipliers
// //   b_i | z_i = k ~ N(0, vb * gamma[k])    for k > 0
// //   b_i | z_i = 0 = 0
// //
// // Annotation model:
// //   p_ij = Phi(A_i alpha_j), j = 1,...,Kgamma-1
// //
// // where p_ij is the conditional probability of stepping from component j-1
// // to a larger component.
// //
// // For Kgamma = 4 components:
// //   pi_i0 = 1 - p_i1
// //   pi_i1 = p_i1 * (1 - p_i2)
// //   pi_i2 = p_i1 * p_i2 * (1 - p_i3)
// //   pi_i3 = p_i1 * p_i2 * p_i3
// //
// // This handles overlapping annotations directly through dense A (m x nAnno).
// // If an intercept is desired, include a first column of 1s in A. The first
// // annotation coefficient is updated with a flat prior when intercept_flat=true.
// //
// // Return structure:
// //   0  bm
// //   1  dm = P(component > 0)
// //   2  wy
// //   3  r
// //   4  b
// //   5  component index, 0-based
// //   6  marker index
// //   7  vbs
// //   8  vgs
// //   9  ves
// //   10 covb
// //   11 covg
// //   12 cove
// //   13 vb
// //   14 vg
// //   15 ve
// //   16 final active pi, length 2: c(pi0, 1-pi0)
// //   17 posterior mean active pi, length 2
// //   18 posterior mean alpha, flattened as nAnno x (Kgamma-1), column-major
// //   19 posterior mean sigmaSqAlpha, length Kgamma-1
// //   20 vle = linkage-equilibrium variance trace
// //   21 vld = linkage-disequilibrium contribution trace, vg - vle
// //
// // Recommended R extraction for slot 18:
// //   alpha <- matrix(fit[[18]][[t]], nrow=ncol(A), ncol=length(gamma)-1)
// //
// // =============================================================================
//
// inline double logsumexp_vec(const std::vector<double>& x) {
//  double mx = -std::numeric_limits<double>::infinity();
//  for (double v : x) mx = std::max(mx, v);
//  if (!std::isfinite(mx)) return mx;
//
//  double s = 0.0;
//  for (double v : x) s += std::exp(v - mx);
//  return mx + std::log(s);
// }
//
// inline int sample_categorical_logprob(
//   const std::vector<double>& logp,
//   std::mt19937& gen
// ) {
//  const double lse = logsumexp_vec(logp);
//  if (!std::isfinite(lse)) {
//   throw std::runtime_error("sample_categorical_logprob: invalid log probability vector.");
//  }
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  const double u = runif(gen);
//  double cum = 0.0;
//
//  for (std::size_t k = 0; k < logp.size(); ++k) {
//   cum += std::exp(logp[k] - lse);
//   if (u <= cum) return static_cast<int>(k);
//  }
//
//  return static_cast<int>(logp.size() - 1);
// }
//
// inline double safe_pnorm(double x) {
//  double p = R::pnorm(x, 0.0, 1.0, 1, 0);
//  if (!std::isfinite(p)) p = (x > 0.0) ? 1.0 : 0.0;
//  return std::min(std::max(p, 1e-12), 1.0 - 1e-12);
// }
//
// inline double safe_qnorm(double p) {
//  p = std::min(std::max(p, 1e-12), 1.0 - 1e-12);
//  return R::qnorm(p, 0.0, 1.0, 1, 0);
// }
//
// inline double rtruncnorm_std(
//   double mu,
//   bool positive,
//   std::mt19937& gen
// ) {
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  const double u = runif(gen);
//
//  if (positive) {
//   // X ~ N(mu,1), truncated X > 0.
//   const double a = safe_pnorm(-mu);
//   const double p = a + u * (1.0 - a);
//   return mu + safe_qnorm(p);
//  }
//
//  // X ~ N(mu,1), truncated X <= 0.
//  const double b = safe_pnorm(-mu);
//  const double p = u * b;
//  return mu + safe_qnorm(p);
// }
//
// inline arma::mat compute_snp_pi_from_alpha(
//   const arma::mat& A,
//   const arma::mat& alpha,
//   double pi_floor
// ) {
//  const int m = static_cast<int>(A.n_rows);
//  const int nstep = static_cast<int>(alpha.n_cols);
//  const int Kgamma = nstep + 1;
//
//  arma::mat p(m, nstep, arma::fill::zeros);
//  arma::mat snpPi(m, Kgamma, arma::fill::zeros);
//
//  for (int j = 0; j < nstep; ++j) {
//   arma::vec eta = A * alpha.col(static_cast<arma::uword>(j));
//   for (int i = 0; i < m; ++i) {
//    p(static_cast<arma::uword>(i), static_cast<arma::uword>(j)) =
//     safe_pnorm(eta(static_cast<arma::uword>(i)));
//   }
//  }
//
//  for (int i = 0; i < m; ++i) {
//   double prod_prev = 1.0;
//
//   for (int k = 0; k < Kgamma; ++k) {
//    double val;
//
//    if (k < nstep) {
//     const double pk = p(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
//     val = prod_prev * (1.0 - pk);
//     prod_prev *= pk;
//    } else {
//     val = prod_prev;
//    }
//
//    snpPi(static_cast<arma::uword>(i), static_cast<arma::uword>(k)) =
//     std::max(val, pi_floor);
//   }
//
//   const double s = arma::accu(snpPi.row(static_cast<arma::uword>(i)));
//   if (!std::isfinite(s) || s <= 0.0) {
//    throw std::runtime_error("compute_snp_pi_from_alpha: invalid row probability sum.");
//   }
//   snpPi.row(static_cast<arma::uword>(i)) /= s;
//  }
//
//  return snpPi;
// }
//
// inline void sampleBeta_SBayesRC_ST_csr(
//   int i,
//   const arma::rowvec& pi_i,
//   const arma::vec& gamma,
//   double vb_t,
//   double vei_i,
//   const arma::rowvec& ww,
//   arma::rowvec& r,
//   arma::rowvec& b,
//   arma::Row<int>& comp,
//   const STLDCSR& ld,
//   std::mt19937& gen
// ) {
//  const arma::uword iu = static_cast<arma::uword>(i);
//  const double wi = ww(iu);
//
//  if (!std::isfinite(wi) || wi <= 0.0) {
//   throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: invalid ww value.");
//  }
//
//  if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//   throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: invalid vb.");
//  }
//
//  const int Kgamma = static_cast<int>(gamma.n_elem);
//  if (static_cast<int>(pi_i.n_elem) != Kgamma) {
//   throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: pi_i/gamma length mismatch.");
//  }
//
//  const double vei_safe = std::max(vei_i, 1e-300);
//  const double score = r(iu) + wi * b(iu);
//
//  std::vector<double> logp(static_cast<std::size_t>(Kgamma));
//
//  for (int k = 0; k < Kgamma; ++k) {
//   const double pik = std::max(static_cast<double>(pi_i(static_cast<arma::uword>(k))), 1e-300);
//   const double gk = gamma(static_cast<arma::uword>(k));
//
//   if (!std::isfinite(gk) || gk < 0.0) {
//    throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: gamma must be non-negative.");
//   }
//
//   if (gk <= 0.0) {
//    logp[static_cast<std::size_t>(k)] = std::log(pik);
//   } else {
//    const double vbk = std::max(vb_t * gk, 1e-300);
//    const double denom = std::max(vei_safe + wi * vbk, 1e-300);
//
//    const double logBF =
//     0.5 * std::log(vei_safe / denom)
//     + 0.5 * score * score * vbk / (vei_safe * denom);
//
//    logp[static_cast<std::size_t>(k)] = std::log(pik) + logBF;
//   }
//  }
//
//  const int k_new = sample_categorical_logprob(logp, gen);
//
//  double b_new = 0.0;
//  const double gamma_new = gamma(static_cast<arma::uword>(k_new));
//
//  if (gamma_new > 0.0) {
//   std::normal_distribution<double> norm01(0.0, 1.0);
//   const double vbk = std::max(vb_t * gamma_new, 1e-300);
//   const double lhs = wi + vei_safe / vbk;
//   const double mean = score / lhs;
//   const double sd = std::sqrt(vei_safe / lhs);
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
//  comp(iu) = k_new;
// }
//
// inline void sampleB_SBayesRC_ST_csr(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& comp,
//   const arma::vec& gamma,
//   double ssb_prior,
//   std::mt19937& gen
// ) {
//  double ssb_scaled = 0.0;
//  double dfb = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//   const int k = comp(iu);
//
//   if (k > 0) {
//    if (k >= static_cast<int>(gamma.n_elem)) {
//     throw std::runtime_error("sampleB_SBayesRC_ST_csr: component index out of range.");
//    }
//
//    const double gk = gamma(static_cast<arma::uword>(k));
//
//    if (!std::isfinite(gk) || gk <= 0.0) {
//     throw std::runtime_error("sampleB_SBayesRC_ST_csr: active component has invalid gamma.");
//    }
//
//    ssb_scaled += b(iu) * b(iu) / gk;
//    dfb += 1.0;
//   }
//  }
//
//  const double scale = ssb_scaled + nub * ssb_prior;
//
//  if (!std::isfinite(scale) || scale <= 0.0) {
//   throw std::runtime_error("sampleB_SBayesRC_ST_csr: invalid scale.");
//  }
//
//  std::chi_squared_distribution<double> rchisq(dfb + nub);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  vb = std::max(scale / chi2, 1e-12);
// }
//
// inline double computeLE_SBayesRC_ST_csr(
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
// inline void build_step_indicators(
//   const arma::Row<int>& comp,
//   arma::Mat<int>& zstep
// ) {
//  const int m = static_cast<int>(comp.n_elem);
//  const int nstep = static_cast<int>(zstep.n_cols);
//
//  for (int i = 0; i < m; ++i) {
//   const int ci = comp(static_cast<arma::uword>(i));
//   for (int j = 0; j < nstep; ++j) {
//    zstep(static_cast<arma::uword>(i), static_cast<arma::uword>(j)) = (ci > j) ? 1 : 0;
//   }
//  }
// }
//
// inline void update_alpha_sbayesrc(
//   const arma::mat& A,
//   const arma::Row<int>& comp,
//   arma::mat& alpha,
//   arma::vec& sigmaSqAlpha,
//   bool intercept_flat,
//   double sigmaSqAlpha_a,
//   double sigmaSqAlpha_b,
//   std::mt19937& gen
// ) {
//  const int m = static_cast<int>(A.n_rows);
//  const int nAnno = static_cast<int>(A.n_cols);
//  const int nstep = static_cast<int>(alpha.n_cols);
//
//  arma::Mat<int> zstep(m, nstep, arma::fill::zeros);
//  build_step_indicators(comp, zstep);
//
//  for (int j = 0; j < nstep; ++j) {
//   std::vector<int> idx;
//   idx.reserve(static_cast<std::size_t>(m));
//
//   for (int i = 0; i < m; ++i) {
//    if (j == 0 || zstep(static_cast<arma::uword>(i), static_cast<arma::uword>(j - 1)) > 0) {
//     idx.push_back(i);
//    }
//   }
//
//   if (idx.empty()) continue;
//
//   const int nj = static_cast<int>(idx.size());
//   arma::vec mu(nj, arma::fill::zeros);
//
//   for (int ii = 0; ii < nj; ++ii) {
//    const int i = idx[static_cast<std::size_t>(ii)];
//    double s = 0.0;
//    for (int k = 0; k < nAnno; ++k) {
//     s += A(static_cast<arma::uword>(i), static_cast<arma::uword>(k)) *
//      alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j));
//    }
//    mu(static_cast<arma::uword>(ii)) = s;
//   }
//
//   arma::vec latent(nj, arma::fill::zeros);
//
//   for (int ii = 0; ii < nj; ++ii) {
//    const int i = idx[static_cast<std::size_t>(ii)];
//    const bool positive = zstep(static_cast<arma::uword>(i), static_cast<arma::uword>(j)) > 0;
//    latent(static_cast<arma::uword>(ii)) = rtruncnorm_std(mu(static_cast<arma::uword>(ii)), positive, gen);
//   }
//
//   // Residualized latent variable, analogous to lj <- lj - mu in the R code.
//   arma::vec resid = latent - mu;
//
//   for (int k = 0; k < nAnno; ++k) {
//    const bool flat_prior = intercept_flat && (k == 0);
//
//    double diag_k = 0.0;
//    double rhs = 0.0;
//
//    const double old = alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j));
//
//    for (int ii = 0; ii < nj; ++ii) {
//     const int i = idx[static_cast<std::size_t>(ii)];
//     const double x = A(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
//     diag_k += x * x;
//     rhs += x * resid(static_cast<arma::uword>(ii));
//    }
//
//    rhs += diag_k * old;
//
//    if (diag_k <= 0.0) continue;
//
//    double invLhs;
//    if (flat_prior) {
//     invLhs = 1.0 / diag_k;
//    } else {
//     const double sig = std::max(sigmaSqAlpha(static_cast<arma::uword>(j)), 1e-12);
//     invLhs = 1.0 / (diag_k + 1.0 / sig);
//    }
//
//    const double mean = invLhs * rhs;
//    const double sd = std::sqrt(invLhs);
//
//    std::normal_distribution<double> norm(mean, sd);
//    const double anew = norm(gen);
//
//    alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j)) = anew;
//
//    const double diff_old_new = old - anew;
//    if (diff_old_new != 0.0) {
//     for (int ii = 0; ii < nj; ++ii) {
//      const int i = idx[static_cast<std::size_t>(ii)];
//      const double x = A(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
//      resid(static_cast<arma::uword>(ii)) += x * diff_old_new;
//     }
//    }
//   }
//
//   // Scaled inverse-chi-square / inverse-gamma style update.
//   // R code: sigmaSqAlpha[j] = (sum(alpha[-1,j]^2) + 2) / rchisq(nAnno-1+2)
//   // Here sigmaSqAlpha_a and sigmaSqAlpha_b generalize that default.
//   double ss = 0.0;
//   int ncoef = 0;
//   for (int k = 0; k < nAnno; ++k) {
//    if (intercept_flat && k == 0) continue;
//    const double ak = alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j));
//    ss += ak * ak;
//    ++ncoef;
//   }
//
//   if (ncoef > 0) {
//    const double df = static_cast<double>(ncoef) + sigmaSqAlpha_a;
//    const double scale = ss + sigmaSqAlpha_b;
//    std::chi_squared_distribution<double> rchisq(df);
//    const double chi2 = std::max(rchisq(gen), 1e-300);
//    sigmaSqAlpha(static_cast<arma::uword>(j)) = std::max(scale / chi2, 1e-12);
//   }
//  }
// }
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_sbayesrc(
//   std::vector<std::vector<double>> wy,
//   std::vector<std::vector<double>> ww,
//   std::vector<double> yy,
//   std::vector<std::vector<double>> b_init,
//   std::vector<std::vector<double>> comp_init,
//   bool use_comp_init,
//   std::vector<std::vector<double>> r_init,
//   bool use_r_init,
//   bool rebuild_r_before_updateE,
//   std::string ld_prefix,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   arma::mat A,
//   arma::vec gamma,
//   arma::mat alpha_init,
//   arma::vec sigmaSqAlpha_init,
//   bool intercept_flat,
//   double sigmaSqAlpha_a,
//   double sigmaSqAlpha_b,
//   double pi_floor,
//   double nub,
//   double nue,
//   bool updateAlpha,
//   bool updateB,
//   bool updateE,
//   int alpha_update_every,
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
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: nt must be positive.");
//  }
//
//  const int m = static_cast<int>(wy[0].size());
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: m must be positive.");
//  }
//
//  if (nit <= 0 || nburn < 0 || nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: invalid nit/nburn/nthin.");
//  }
//
//  if (alpha_update_every <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_update_every must be positive.");
//  }
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt ||
//      (int)yy.size() != nt || (int)n.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: inconsistent trait dimensions.");
//  }
//
//  if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: priors must be nt x nt.");
//  }
//
//  if ((int)A.n_rows != m) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: A must have m rows.");
//  }
//
//  const int nAnno = static_cast<int>(A.n_cols);
//  if (nAnno <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: A must have at least one column.");
//  }
//
//  const int Kgamma = static_cast<int>(gamma.n_elem);
//  if (Kgamma < 2) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma must have at least two components, including zero.");
//  }
//
//  if (!std::isfinite(gamma(0)) || gamma(0) != 0.0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma[0] must be exactly 0.0.");
//  }
//
//  for (int k = 1; k < Kgamma; ++k) {
//   const double g = gamma(static_cast<arma::uword>(k));
//   if (!std::isfinite(g) || g <= 0.0) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma[k] must be positive for k > 0.");
//   }
//  }
//
//  const int nstep = Kgamma - 1;
//
//  if ((int)alpha_init.n_rows != nAnno || (int)alpha_init.n_cols != nstep) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_init must be ncol(A) x (length(gamma)-1).");
//  }
//
//  if ((int)sigmaSqAlpha_init.n_elem != nstep) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: sigmaSqAlpha_init must have length length(gamma)-1.");
//  }
//
//  if (!std::isfinite(pi_floor) || pi_floor <= 0.0 || pi_floor >= 1.0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: pi_floor must be in (0,1).");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m ||
//       (int)ww[t].size() != m ||
//       (int)b_init[t].size() != m) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: inconsistent marker dimensions.");
//   }
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: E must be nt x nt.");
//  }
//
//  if (use_r_init) {
//   if (static_cast<int>(r_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: r_init must have length nt when use_r_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(r_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: each r_init[t] must have length m.");
//    }
//   }
//  }
//
//  if (use_comp_init) {
//   if (static_cast<int>(comp_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: comp_init must have length nt when enabled.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(comp_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: comp_init[t] must have length m.");
//    }
//   }
//  }
//
//  arma::mat wy_mat(nt, m, arma::fill::zeros);
//  arma::mat ww_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::mat r_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> comp_mat(nt, m, arma::fill::zeros);
//
//  arma::vec yy_vec(nt, arma::fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   yy_vec(static_cast<arma::uword>(t)) = yy[t];
//
//   for (int i = 0; i < m; ++i) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword iu = static_cast<arma::uword>(i);
//
//    wy_mat(tu, iu) = wy[t][i];
//    ww_mat(tu, iu) = ww[t][i];
//    b_mat(tu, iu)  = b_init[t][i];
//   }
//
//   if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
//    sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
//   }
//  }
//
//  for (int t = 1; t < nt; ++t) {
//   if (n[t] != n[0]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_csr_sbayesrc: current shared-LD scaling assumes equal n across traits."
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
//       "stblr_cpg_omp_csr_sbayesrc: ww contains invalid value before LD pre-scaling."
//     );
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr_sbayesrc: ww differs across traits; pre-scaled shared ST LD is invalid."
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
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: ww contains invalid value in trait 0.");
//   }
//   xx[static_cast<std::size_t>(i)] = wi;
//  }
//
//  STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
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
//            [&](int a, int b) {
//             return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
//            });
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
//  arma::vec final_pi_active(nt, arma::fill::zeros);
//  arma::vec final_vle(nt, arma::fill::zeros);
//  arma::vec final_vld(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  std::vector<arma::mat> alpha_mean(static_cast<std::size_t>(nt));
//  std::vector<arma::vec> sigmaSqAlpha_mean(static_cast<std::size_t>(nt));
//  for (int t = 0; t < nt; ++t) {
//   alpha_mean[static_cast<std::size_t>(t)] = arma::mat(nAnno, nstep, arma::fill::zeros);
//   sigmaSqAlpha_mean[static_cast<std::size_t>(t)] = arma::vec(nstep, arma::fill::zeros);
//  }
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
//  nthreads = std::max(1, std::min(ncores, nt));
//  omp_set_num_threads(nthreads);
//
//  Rcpp::Rcout
//  << "STBLR SBayesRC CSR OpenMP requested threads = "
//  << nthreads
//  << ", omp_get_max_threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
//  Rcpp::Rcout
//  << "STBLR real-SBayesRC CSR: m=" << m
//  << ", nt=" << nt
//  << ", annotations=" << nAnno
//  << ", components=" << Kgamma
//  << ", updateAlpha=" << updateAlpha
//  << ", alpha_update_every=" << alpha_update_every
//  << "\n";
//
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(nthreads) schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
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
//    arma::Row<int> comp_t(m, arma::fill::zeros);
//
//    if (use_comp_init) {
//     for (int i = 0; i < m; ++i) {
//      const int k = static_cast<int>(std::round(comp_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)]));
//      if (k < 0 || k >= Kgamma) {
//       throw std::runtime_error("comp_init contains component outside 0..Kgamma-1.");
//      }
//      comp_t(static_cast<arma::uword>(i)) = k;
//     }
//    } else {
//     for (int i = 0; i < m; ++i) {
//      comp_t(static_cast<arma::uword>(i)) =
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
//      throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: r_init contains NaN/Inf.");
//     }
//    } else {
//     rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
//    }
//
//    double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
//    double vle_t = computeLE_SBayesRC_ST_csr(m, b_t, ww_t, n[t]);
//    double vld_t = vg_t - vle_t;
//    double vei_t = ve_t + adjE * vg_t;
//
//    arma::mat alpha_t = alpha_init;
//    arma::vec sigmaSqAlpha_t = sigmaSqAlpha_init;
//
//    for (int j = 0; j < nstep; ++j) {
//     if (!std::isfinite(sigmaSqAlpha_t(static_cast<arma::uword>(j))) ||
//         sigmaSqAlpha_t(static_cast<arma::uword>(j)) <= 0.0) {
//      throw std::runtime_error("sigmaSqAlpha_init contains invalid value.");
//     }
//    }
//
//    arma::mat snpPi_t = compute_snp_pi_from_alpha(A, alpha_t, pi_floor);
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
//    arma::mat alpha_accum(nAnno, nstep, arma::fill::zeros);
//    arma::vec sigmaSqAlpha_accum(nstep, arma::fill::zeros);
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//     for (int isort = 0; isort < m; ++isort) {
//      const int i = order[static_cast<std::size_t>(isort)];
//
//      sampleBeta_SBayesRC_ST_csr(
//       i,
//       snpPi_t.row(static_cast<arma::uword>(i)),
//       gamma,
//       vb_t,
//       vei_t,
//       ww_t,
//       r_t,
//       b_t,
//       comp_t,
//       ld,
//       gen_t
//      );
//     }
//
//     if (updateAlpha && ((it + 1) % alpha_update_every == 0)) {
//      update_alpha_sbayesrc(
//       A,
//       comp_t,
//       alpha_t,
//       sigmaSqAlpha_t,
//       intercept_flat,
//       sigmaSqAlpha_a,
//       sigmaSqAlpha_b,
//       gen_t
//      );
//
//      snpPi_t = compute_snp_pi_from_alpha(A, alpha_t, pi_floor);
//     }
//
//     if (updateB) {
//      sampleB_SBayesRC_ST_csr(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       comp_t,
//       gamma,
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
//       rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
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
//     vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
//     vle_t = computeLE_SBayesRC_ST_csr(m, b_t, ww_t, n[t]);
//     vld_t = vg_t - vle_t;
//
//     if (!std::isfinite(vg_t)) {
//      throw std::runtime_error("vg became NaN/Inf after computeG. iter=" + std::to_string(it));
//     }
//
//     if (!std::isfinite(vle_t)) {
//      throw std::runtime_error("vle became NaN/Inf after computeLE. iter=" + std::to_string(it));
//     }
//
//     if (!std::isfinite(vld_t)) {
//      throw std::runtime_error("vld became NaN/Inf after computeLE. iter=" + std::to_string(it));
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
//     double pi_active = 0.0;
//     for (int i = 0; i < m; ++i) {
//      double row_active = 0.0;
//      for (int k = 1; k < Kgamma; ++k) {
//       row_active += snpPi_t(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
//      }
//      pi_active += row_active;
//     }
//     pi_active /= static_cast<double>(m);
//
//     vbs_t(static_cast<arma::uword>(it)) = vb_t;
//     ves_t(static_cast<arma::uword>(it)) = ve_t;
//     vgs_t(static_cast<arma::uword>(it)) = vg_t;
//     pis_t(static_cast<arma::uword>(it)) = pi_active;
//     vles_t(static_cast<arma::uword>(it)) = vle_t;
//     vlds_t(static_cast<arma::uword>(it)) = vld_t;
//
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       const arma::uword iu = static_cast<arma::uword>(i);
//       bm_t(iu) += b_t(iu);
//       dm_t(iu) += (comp_t(iu) > 0) ? 1.0 : 0.0;
//      }
//
//      alpha_accum += alpha_t;
//      sigmaSqAlpha_accum += sigmaSqAlpha_t;
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    bm_t /= nsamples_t;
//    dm_t /= nsamples_t;
//    alpha_accum /= nsamples_t;
//    sigmaSqAlpha_accum /= nsamples_t;
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
//    comp_mat.row(static_cast<arma::uword>(t)) = comp_t;
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
//    final_pi_active(static_cast<arma::uword>(t)) = pis_t(static_cast<arma::uword>(nit + nburn - 1));
//    final_vle(static_cast<arma::uword>(t)) = vle_t;
//    final_vld(static_cast<arma::uword>(t)) = vld_t;
//    nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;
//
//    alpha_mean[static_cast<std::size_t>(t)] = alpha_accum;
//    sigmaSqAlpha_mean[static_cast<std::size_t>(t)] = sigmaSqAlpha_accum;
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
//      "stblr_cpg_omp_csr_sbayesrc failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[static_cast<std::size_t>(t)]
//    );
//   }
//  }
//
//  std::vector<std::vector<std::vector<double>>> result(22);
//
//  for (int k = 0; k < 22; ++k) {
//   result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[1][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[2][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[3][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[4][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[5][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[6][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//
//   result[7][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
//   result[8][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
//   result[9][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
//
//   result[10][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[11][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[12][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[13][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[14][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[15][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//
//   result[16][static_cast<std::size_t>(t)].resize(2);
//   result[17][static_cast<std::size_t>(t)].resize(2);
//
//   result[18][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nAnno * nstep));
//   result[19][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nstep));
//   result[20][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
//   result[21][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
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
//    result[5][ts][is] = static_cast<double>(comp_mat(tu, iu));
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
//   result[16][ts][0] = 1.0 - final_pi_active(static_cast<arma::uword>(t));
//   result[16][ts][1] = final_pi_active(static_cast<arma::uword>(t));
//
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi_active(static_cast<arma::uword>(t));
//
//   result[17][ts][0] = 1.0 - mean_pi;
//   result[17][ts][1] = mean_pi;
//
//   // Flatten nAnno x nstep alpha matrix in column-major order, matching R matrix().
//   for (int j = 0; j < nstep; ++j) {
//    for (int a = 0; a < nAnno; ++a) {
//     const std::size_t pos = static_cast<std::size_t>(j * nAnno + a);
//     result[18][ts][pos] = alpha_mean[ts](static_cast<arma::uword>(a), static_cast<arma::uword>(j));
//    }
//   }
//
//   for (int j = 0; j < nstep; ++j) {
//    result[19][ts][static_cast<std::size_t>(j)] =
//     sigmaSqAlpha_mean[ts](static_cast<arma::uword>(j));
//   }
//  }
//
//  return result;
// }

// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
//
// #include "st_csr_common.h"
//
// #include <algorithm>
// #include <cmath>
// #include <cstdio>
// #include <cstdint>
// #include <cstdlib>
// #include <fstream>
// #include <limits>
// #include <numeric>
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // =============================================================================
// // STBLR summary-stat CSR: SBayesRC-style sampler with overlapping annotations
// // =============================================================================
// //
// // Exported function:
// //
// //   stblr_cpg_omp_csr_sbayesrc(...)
// //
// // This implements the SBayesRC-style annotation model from the R code:
// //
// //   gamma = c(0, 0.01, 0.1, 1)             fixed mixture variance multipliers
// //   b_i | z_i = k ~ N(0, vb * gamma[k])    for k > 0
// //   b_i | z_i = 0 = 0
// //
// // Annotation model:
// //
// //   p_ij = Phi(A_i alpha_j), j = 1,...,Kgamma-1
// //
// // where p_ij is the conditional probability of stepping from component j-1
// // to a larger component.
// //
// // For Kgamma = 4 components:
// //
// //   pi_i0 = 1 - p_i1
// //   pi_i1 = p_i1 * (1 - p_i2)
// //   pi_i2 = p_i1 * p_i2 * (1 - p_i3)
// //   pi_i3 = p_i1 * p_i2 * p_i3
// //
// // This handles overlapping annotations directly through dense A (m x nAnno).
// // If an intercept is desired, include a first column of 1s in A. The first
// // annotation coefficient is updated with a flat prior when intercept_flat=true.
// //
// // Return structure keeps the usual 20 slots:
// //   0  bm
// //   1  dm = P(component > 0)
// //   2  wy
// //   3  r
// //   4  b
// //   5  component index, 0-based
// //   6  marker index
// //   7  vbs
// //   8  vgs
// //   9  ves
// //   10 covb
// //   11 covg
// //   12 cove
// //   13 vb
// //   14 vg
// //   15 ve
// //   16 final active pi, length 2: c(pi0, 1-pi0)
// //   17 posterior mean active pi, length 2
// //   18 posterior mean alpha, flattened as nAnno x (Kgamma-1), column-major
// //   19 posterior mean sigmaSqAlpha, length Kgamma-1
// //
// // Recommended R extraction for slot 18:
// //   alpha <- matrix(fit[[18]][[t]], nrow=ncol(A), ncol=length(gamma)-1)
// //
// // =============================================================================
//
// // struct STLDCSR {
// //  std::vector<uint64_t> ptr;
// //  std::vector<int> idx;
// //  std::vector<float> xij;
// // };
// //
// // inline void read_exact_file(
// //   const std::string& path,
// //   void* data,
// //   std::size_t nbytes
// // ) {
// //  FILE* fs = std::fopen(path.c_str(), "rb");
// //
// //  if (!fs) {
// //   throw std::runtime_error("Could not open file: " + path);
// //  }
// //
// //  const std::size_t got = std::fread(data, 1, nbytes, fs);
// //  std::fclose(fs);
// //
// //  if (got != nbytes) {
// //   throw std::runtime_error("Short read from file: " + path);
// //  }
// // }
// //
// // inline uint64_t parse_uint64_from_meta(
// //   const std::string& value,
// //   const std::string& key
// // ) {
// //  if (value.empty()) {
// //   throw std::runtime_error("Empty metadata value for key: " + key);
// //  }
// //
// //  char* endptr = nullptr;
// //  const unsigned long long out = std::strtoull(value.c_str(), &endptr, 10);
// //
// //  if (endptr == value.c_str() || *endptr != '\0') {
// //   throw std::runtime_error("Invalid unsigned integer metadata value for key: " + key);
// //  }
// //
// //  return static_cast<uint64_t>(out);
// // }
// //
// // inline STLDCSR read_and_build_st_ld_csr(
// //   const std::string& prefix,
// //   int m,
// //   const std::vector<double>& xx
// // ) {
// //  const std::string row_file  = prefix + ".row_ptr.u64.bin";
// //  const std::string col_file  = prefix + ".col_idx.u32.0based.bin";
// //  const std::string val_file  = prefix + ".values.f32.bin";
// //  const std::string meta_file = prefix + ".meta.txt";
// //
// //  if (m <= 0) {
// //   throw std::runtime_error("read_and_build_st_ld_csr: m must be positive.");
// //  }
// //
// //  if (static_cast<int>(xx.size()) != m) {
// //   throw std::runtime_error("read_and_build_st_ld_csr: xx must have length m.");
// //  }
// //
// //  std::ifstream meta(meta_file.c_str());
// //  if (!meta.is_open()) {
// //   throw std::runtime_error("Could not open metadata file: " + meta_file);
// //  }
// //
// //  int m_meta = -1;
// //  uint64_t nnz_u64 = 0;
// //  bool have_nnz = false;
// //
// //  std::string line;
// //  while (std::getline(meta, line)) {
// //   const std::string key_m   = "n_variants=";
// //   const std::string key_nnz = "nnz=";
// //
// //   if (line.rfind(key_m, 0) == 0) {
// //    m_meta = std::stoi(line.substr(key_m.size()));
// //   } else if (line.rfind(key_nnz, 0) == 0) {
// //    nnz_u64 = parse_uint64_from_meta(line.substr(key_nnz.size()), "nnz");
// //    have_nnz = true;
// //   }
// //  }
// //  meta.close();
// //
// //  if (m_meta <= 0) {
// //   throw std::runtime_error("Could not read n_variants from metadata.");
// //  }
// //
// //  if (m_meta != m) {
// //   throw std::runtime_error("LD metadata n_variants does not match marker dimension.");
// //  }
// //
// //  if (!have_nnz) {
// //   throw std::runtime_error("Could not read nnz from metadata.");
// //  }
// //
// //  const std::size_t nnz = static_cast<std::size_t>(nnz_u64);
// //
// //  std::vector<uint64_t> row_ptr(static_cast<std::size_t>(m) + 1);
// //  std::vector<uint32_t> col_idx_u32(nnz);
// //  std::vector<float> values_r(nnz);
// //
// //  read_exact_file(row_file, row_ptr.data(), row_ptr.size() * sizeof(uint64_t));
// //  read_exact_file(col_file, col_idx_u32.data(), col_idx_u32.size() * sizeof(uint32_t));
// //  read_exact_file(val_file, values_r.data(), values_r.size() * sizeof(float));
// //
// //  if (row_ptr[0] != 0 || row_ptr[static_cast<std::size_t>(m)] != nnz_u64) {
// //   throw std::runtime_error("Invalid LD row_ptr: expected 0-based row_ptr ending at nnz.");
// //  }
// //
// //  for (int i = 0; i < m; ++i) {
// //   if (row_ptr[static_cast<std::size_t>(i + 1)] < row_ptr[static_cast<std::size_t>(i)]) {
// //    throw std::runtime_error("Invalid LD row_ptr: row pointers are not nondecreasing.");
// //   }
// //
// //   if (!std::isfinite(xx[static_cast<std::size_t>(i)]) || xx[static_cast<std::size_t>(i)] <= 0.0) {
// //    throw std::runtime_error(
// //      "read_and_build_st_ld_csr: xx contains invalid value at marker " +
// //       std::to_string(i)
// //    );
// //   }
// //  }
// //
// //  std::vector<uint64_t> degree(static_cast<std::size_t>(m), 0);
// //
// //  for (int i = 0; i < m; ++i) {
// //   const uint64_t start = row_ptr[static_cast<std::size_t>(i)];
// //   const uint64_t end   = row_ptr[static_cast<std::size_t>(i + 1)];
// //
// //   if (end > nnz_u64) {
// //    throw std::runtime_error("Invalid LD row_ptr: row end exceeds nnz.");
// //   }
// //
// //   for (uint64_t p = start; p < end; ++p) {
// //    const uint32_t j_u32 = col_idx_u32[static_cast<std::size_t>(p)];
// //
// //    if (j_u32 >= static_cast<uint32_t>(m)) {
// //     throw std::runtime_error("LD column index out of range.");
// //    }
// //
// //    const int j = static_cast<int>(j_u32);
// //
// //    if (j == i) continue;
// //
// //    ++degree[static_cast<std::size_t>(i)];
// //    ++degree[static_cast<std::size_t>(j)];
// //   }
// //  }
// //
// //  STLDCSR ld;
// //  ld.ptr.resize(static_cast<std::size_t>(m) + 1);
// //  ld.ptr[0] = 0;
// //
// //  for (int i = 0; i < m; ++i) {
// //   ld.ptr[static_cast<std::size_t>(i + 1)] =
// //    ld.ptr[static_cast<std::size_t>(i)] + degree[static_cast<std::size_t>(i)];
// //  }
// //
// //  const uint64_t nnz_sym = ld.ptr[static_cast<std::size_t>(m)];
// //
// //  ld.idx.resize(static_cast<std::size_t>(nnz_sym));
// //  ld.xij.resize(static_cast<std::size_t>(nnz_sym));
// //
// //  std::vector<uint64_t> offset = ld.ptr;
// //
// //  double max_abs_rij = 0.0;
// //  double max_abs_xij = 0.0;
// //
// //  for (int i = 0; i < m; ++i) {
// //   const uint64_t start = row_ptr[static_cast<std::size_t>(i)];
// //   const uint64_t end   = row_ptr[static_cast<std::size_t>(i + 1)];
// //
// //   for (uint64_t p = start; p < end; ++p) {
// //    const int j = static_cast<int>(col_idx_u32[static_cast<std::size_t>(p)]);
// //
// //    if (j == i) continue;
// //
// //    const double rij = static_cast<double>(values_r[static_cast<std::size_t>(p)]);
// //
// //    if (!std::isfinite(rij)) {
// //     throw std::runtime_error("LD value contains NaN/Inf.");
// //    }
// //
// //    max_abs_rij = std::max(max_abs_rij, std::abs(rij));
// //
// //    if (std::abs(rij) > 1.0001) {
// //     throw std::runtime_error(
// //       "LD value is not a correlation. Did you pass X_i'X_j instead of r_ij?"
// //     );
// //    }
// //
// //    const double xij =
// //     rij * std::sqrt(xx[static_cast<std::size_t>(i)] * xx[static_cast<std::size_t>(j)]);
// //
// //    if (!std::isfinite(xij)) {
// //     throw std::runtime_error("Computed X_i'X_j contains NaN/Inf.");
// //    }
// //
// //    max_abs_xij = std::max(max_abs_xij, std::abs(xij));
// //
// //    const float xij_f = static_cast<float>(xij);
// //
// //    const uint64_t pos_i = offset[static_cast<std::size_t>(i)]++;
// //    ld.idx[static_cast<std::size_t>(pos_i)] = j;
// //    ld.xij[static_cast<std::size_t>(pos_i)] = xij_f;
// //
// //    const uint64_t pos_j = offset[static_cast<std::size_t>(j)]++;
// //    ld.idx[static_cast<std::size_t>(pos_j)] = i;
// //    ld.xij[static_cast<std::size_t>(pos_j)] = xij_f;
// //   }
// //  }
// //
// //  for (int i = 0; i < m; ++i) {
// //   if (offset[static_cast<std::size_t>(i)] != ld.ptr[static_cast<std::size_t>(i + 1)]) {
// //    throw std::runtime_error("Internal LD CSR fill-count mismatch.");
// //   }
// //  }
// //
// //  Rcpp::Rcout
// //  << "ST SBayesRC flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
// //  << ", symmetric nnz=" << static_cast<double>(nnz_sym)
// //  << ", max_abs_rij=" << max_abs_rij
// //  << ", max_abs_xij=" << max_abs_xij
// //  << "\n";
// //
// //  return ld;
// // }
// //
// // inline void rebuild_residual_st_csr(
// //   int m,
// //   const arma::rowvec& wy,
// //   const arma::rowvec& ww,
// //   const arma::rowvec& b,
// //   arma::rowvec& r,
// //   const STLDCSR& ld
// // ) {
// //  r = wy;
// //
// //  for (int i = 0; i < m; ++i) {
// //   const arma::uword iu = static_cast<arma::uword>(i);
// //   const double bi = b(iu);
// //
// //   if (bi == 0.0) continue;
// //
// //   r(iu) -= ww(iu) * bi;
// //
// //   const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
// //   const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];
// //
// //   for (uint64_t p = start; p < end; ++p) {
// //    const int j = ld.idx[static_cast<std::size_t>(p)];
// //    r(static_cast<arma::uword>(j)) -=
// //     static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * bi;
// //   }
// //  }
// // }
//
// inline double logsumexp_vec(const std::vector<double>& x) {
//  double mx = -std::numeric_limits<double>::infinity();
//  for (double v : x) mx = std::max(mx, v);
//  if (!std::isfinite(mx)) return mx;
//  double s = 0.0;
//  for (double v : x) s += std::exp(v - mx);
//  return mx + std::log(s);
// }
//
// inline int sample_categorical_logprob(
//   const std::vector<double>& logp,
//   std::mt19937& gen
// ) {
//  const double lse = logsumexp_vec(logp);
//  if (!std::isfinite(lse)) {
//   throw std::runtime_error("sample_categorical_logprob: invalid log probability vector.");
//  }
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  const double u = runif(gen);
//  double cum = 0.0;
//
//  for (std::size_t k = 0; k < logp.size(); ++k) {
//   cum += std::exp(logp[k] - lse);
//   if (u <= cum) return static_cast<int>(k);
//  }
//
//  return static_cast<int>(logp.size() - 1);
// }
//
// inline double safe_pnorm(double x) {
//  double p = R::pnorm(x, 0.0, 1.0, 1, 0);
//  if (!std::isfinite(p)) p = (x > 0.0) ? 1.0 : 0.0;
//  return std::min(std::max(p, 1e-12), 1.0 - 1e-12);
// }
//
// inline double safe_qnorm(double p) {
//  p = std::min(std::max(p, 1e-12), 1.0 - 1e-12);
//  return R::qnorm(p, 0.0, 1.0, 1, 0);
// }
//
// inline double rtruncnorm_std(
//   double mu,
//   bool positive,
//   std::mt19937& gen
// ) {
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  const double u = runif(gen);
//
//  if (positive) {
//   // X ~ N(mu,1), truncated X > 0.
//   const double a = safe_pnorm(-mu);
//   const double p = a + u * (1.0 - a);
//   return mu + safe_qnorm(p);
//  }
//
//  // X ~ N(mu,1), truncated X <= 0.
//  const double b = safe_pnorm(-mu);
//  const double p = u * b;
//  return mu + safe_qnorm(p);
// }
//
// inline arma::mat compute_snp_pi_from_alpha(
//   const arma::mat& A,
//   const arma::mat& alpha,
//   double pi_floor
// ) {
//  const int m = static_cast<int>(A.n_rows);
//  const int nstep = static_cast<int>(alpha.n_cols);
//  const int Kgamma = nstep + 1;
//
//  arma::mat p(m, nstep, arma::fill::zeros);
//  arma::mat snpPi(m, Kgamma, arma::fill::zeros);
//
//  for (int j = 0; j < nstep; ++j) {
//   arma::vec eta = A * alpha.col(static_cast<arma::uword>(j));
//   for (int i = 0; i < m; ++i) {
//    p(static_cast<arma::uword>(i), static_cast<arma::uword>(j)) =
//     safe_pnorm(eta(static_cast<arma::uword>(i)));
//   }
//  }
//
//  for (int i = 0; i < m; ++i) {
//   double prod_prev = 1.0;
//
//   for (int k = 0; k < Kgamma; ++k) {
//    double val;
//
//    if (k < nstep) {
//     const double pk = p(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
//     val = prod_prev * (1.0 - pk);
//     prod_prev *= pk;
//    } else {
//     val = prod_prev;
//    }
//
//    snpPi(static_cast<arma::uword>(i), static_cast<arma::uword>(k)) =
//     std::max(val, pi_floor);
//   }
//
//   const double s = arma::accu(snpPi.row(static_cast<arma::uword>(i)));
//   if (!std::isfinite(s) || s <= 0.0) {
//    throw std::runtime_error("compute_snp_pi_from_alpha: invalid row probability sum.");
//   }
//   snpPi.row(static_cast<arma::uword>(i)) /= s;
//  }
//
//  return snpPi;
// }
//
// inline void sampleBeta_SBayesRC_ST_csr(
//   int i,
//   const arma::rowvec& pi_i,
//   const arma::vec& gamma,
//   double vb_t,
//   double vei_i,
//   const arma::rowvec& ww,
//   arma::rowvec& r,
//   arma::rowvec& b,
//   arma::Row<int>& comp,
//   const STLDCSR& ld,
//   std::mt19937& gen
// ) {
//  const arma::uword iu = static_cast<arma::uword>(i);
//  const double wi = ww(iu);
//
//  if (!std::isfinite(wi) || wi <= 0.0) {
//   throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: invalid ww value.");
//  }
//
//  if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//   throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: invalid vb.");
//  }
//
//  const int Kgamma = static_cast<int>(gamma.n_elem);
//  if (static_cast<int>(pi_i.n_elem) != Kgamma) {
//   throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: pi_i/gamma length mismatch.");
//  }
//
//  const double vei_safe = std::max(vei_i, 1e-300);
//  const double score = r(iu) + wi * b(iu);
//
//  std::vector<double> logp(static_cast<std::size_t>(Kgamma));
//
//  for (int k = 0; k < Kgamma; ++k) {
//   const double pik = std::max(static_cast<double>(pi_i(static_cast<arma::uword>(k))), 1e-300);
//   const double gk = gamma(static_cast<arma::uword>(k));
//
//   if (!std::isfinite(gk) || gk < 0.0) {
//    throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: gamma must be non-negative.");
//   }
//
//   if (gk <= 0.0) {
//    logp[static_cast<std::size_t>(k)] = std::log(pik);
//   } else {
//    const double vbk = std::max(vb_t * gk, 1e-300);
//    const double denom = std::max(vei_safe + wi * vbk, 1e-300);
//
//    const double logBF =
//     0.5 * std::log(vei_safe / denom)
//     + 0.5 * score * score * vbk / (vei_safe * denom);
//
//    logp[static_cast<std::size_t>(k)] = std::log(pik) + logBF;
//   }
//  }
//
//  const int k_new = sample_categorical_logprob(logp, gen);
//
//  double b_new = 0.0;
//  const double gamma_new = gamma(static_cast<arma::uword>(k_new));
//
//  if (gamma_new > 0.0) {
//   std::normal_distribution<double> norm01(0.0, 1.0);
//   const double vbk = std::max(vb_t * gamma_new, 1e-300);
//   const double lhs = wi + vei_safe / vbk;
//   const double mean = score / lhs;
//   const double sd = std::sqrt(vei_safe / lhs);
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
//  comp(iu) = k_new;
// }
//
// inline void sampleB_SBayesRC_ST_csr(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& comp,
//   const arma::vec& gamma,
//   double ssb_prior,
//   std::mt19937& gen
// ) {
//  double ssb_scaled = 0.0;
//  double dfb = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//   const int k = comp(iu);
//
//   if (k > 0) {
//    if (k >= static_cast<int>(gamma.n_elem)) {
//     throw std::runtime_error("sampleB_SBayesRC_ST_csr: component index out of range.");
//    }
//
//    const double gk = gamma(static_cast<arma::uword>(k));
//
//    if (!std::isfinite(gk) || gk <= 0.0) {
//     throw std::runtime_error("sampleB_SBayesRC_ST_csr: active component has invalid gamma.");
//    }
//
//    ssb_scaled += b(iu) * b(iu) / gk;
//    dfb += 1.0;
//   }
//  }
//
//  const double scale = ssb_scaled + nub * ssb_prior;
//
//  if (!std::isfinite(scale) || scale <= 0.0) {
//   throw std::runtime_error("sampleB_SBayesRC_ST_csr: invalid scale.");
//  }
//
//  std::chi_squared_distribution<double> rchisq(dfb + nub);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  vb = std::max(scale / chi2, 1e-12);
// }
//
// // inline void sampleE_ST_csr(
// //   int m,
// //   double nue,
// //   double& ve,
// //   const arma::rowvec& b,
// //   const arma::rowvec& wy,
// //   const arma::rowvec& r,
// //   double sse_prior,
// //   double yy,
// //   int n,
// //   std::mt19937& gen
// // ) {
// //  double b_dot_r_plus_wy = 0.0;
// //
// //  for (int i = 0; i < m; ++i) {
// //   const arma::uword iu = static_cast<arma::uword>(i);
// //   b_dot_r_plus_wy += b(iu) * (r(iu) + wy(iu));
// //  }
// //
// //  const double sse = yy - b_dot_r_plus_wy;
// //  const double scale = sse + nue * sse_prior;
// //
// //  if (!std::isfinite(scale) || scale <= 0.0) {
// //   throw std::runtime_error("sampleE_ST_csr: invalid residual scale.");
// //  }
// //
// //  std::chi_squared_distribution<double> rchisq(n + nue);
// //  const double chi2 = std::max(rchisq(gen), 1e-300);
// //
// //  const double ve_new = scale / chi2;
// //
// //  if (!std::isfinite(ve_new) || ve_new <= 0.0) {
// //   throw std::runtime_error("sampleE_ST_csr: sampled ve is invalid.");
// //  }
// //
// //  ve = std::max(ve_new, 1e-12);
// // }
//
// // inline double computeG_ST_csr(
// //   const arma::rowvec& b,
// //   const arma::rowvec& wy,
// //   const arma::rowvec& r,
// //   int n
// // ) {
// //  double ssg = 0.0;
// //  const arma::uword m = b.n_elem;
// //
// //  for (arma::uword i = 0; i < m; ++i) {
// //   ssg += b(i) * (wy(i) - r(i));
// //  }
// //
// //  return ssg / static_cast<double>(n);
// // }
// //
// inline void build_step_indicators(
//   const arma::Row<int>& comp,
//   arma::Mat<int>& zstep
// ) {
//  const int m = static_cast<int>(comp.n_elem);
//  const int nstep = static_cast<int>(zstep.n_cols);
//
//  for (int i = 0; i < m; ++i) {
//   const int ci = comp(static_cast<arma::uword>(i));
//   for (int j = 0; j < nstep; ++j) {
//    zstep(static_cast<arma::uword>(i), static_cast<arma::uword>(j)) = (ci > j) ? 1 : 0;
//   }
//  }
// }
//
// inline void update_alpha_sbayesrc(
//   const arma::mat& A,
//   const arma::Row<int>& comp,
//   arma::mat& alpha,
//   arma::vec& sigmaSqAlpha,
//   bool intercept_flat,
//   double sigmaSqAlpha_a,
//   double sigmaSqAlpha_b,
//   std::mt19937& gen
// ) {
//  const int m = static_cast<int>(A.n_rows);
//  const int nAnno = static_cast<int>(A.n_cols);
//  const int nstep = static_cast<int>(alpha.n_cols);
//
//  arma::Mat<int> zstep(m, nstep, arma::fill::zeros);
//  build_step_indicators(comp, zstep);
//
//  for (int j = 0; j < nstep; ++j) {
//   std::vector<int> idx;
//   idx.reserve(static_cast<std::size_t>(m));
//
//   for (int i = 0; i < m; ++i) {
//    if (j == 0 || zstep(static_cast<arma::uword>(i), static_cast<arma::uword>(j - 1)) > 0) {
//     idx.push_back(i);
//    }
//   }
//
//   if (idx.empty()) continue;
//
//   const int nj = static_cast<int>(idx.size());
//   arma::vec mu(nj, arma::fill::zeros);
//
//   for (int ii = 0; ii < nj; ++ii) {
//    const int i = idx[static_cast<std::size_t>(ii)];
//    double s = 0.0;
//    for (int k = 0; k < nAnno; ++k) {
//     s += A(static_cast<arma::uword>(i), static_cast<arma::uword>(k)) *
//      alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j));
//    }
//    mu(static_cast<arma::uword>(ii)) = s;
//   }
//
//   arma::vec latent(nj, arma::fill::zeros);
//
//   for (int ii = 0; ii < nj; ++ii) {
//    const int i = idx[static_cast<std::size_t>(ii)];
//    const bool positive = zstep(static_cast<arma::uword>(i), static_cast<arma::uword>(j)) > 0;
//    latent(static_cast<arma::uword>(ii)) = rtruncnorm_std(mu(static_cast<arma::uword>(ii)), positive, gen);
//   }
//
//   // residualized latent variable, analogous to lj <- lj - mu in the R code.
//   arma::vec resid = latent - mu;
//
//   for (int k = 0; k < nAnno; ++k) {
//    const bool flat_prior = intercept_flat && (k == 0);
//
//    double diag_k = 0.0;
//    double rhs = 0.0;
//
//    const double old = alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j));
//
//    for (int ii = 0; ii < nj; ++ii) {
//     const int i = idx[static_cast<std::size_t>(ii)];
//     const double x = A(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
//     diag_k += x * x;
//     rhs += x * resid(static_cast<arma::uword>(ii));
//    }
//
//    rhs += diag_k * old;
//
//    if (diag_k <= 0.0) continue;
//
//    double invLhs;
//    if (flat_prior) {
//     invLhs = 1.0 / diag_k;
//    } else {
//     const double sig = std::max(sigmaSqAlpha(static_cast<arma::uword>(j)), 1e-12);
//     invLhs = 1.0 / (diag_k + 1.0 / sig);
//    }
//
//    const double mean = invLhs * rhs;
//    const double sd = std::sqrt(invLhs);
//
//    std::normal_distribution<double> norm(mean, sd);
//    const double anew = norm(gen);
//
//    alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j)) = anew;
//
//    const double diff_old_new = old - anew;
//    if (diff_old_new != 0.0) {
//     for (int ii = 0; ii < nj; ++ii) {
//      const int i = idx[static_cast<std::size_t>(ii)];
//      const double x = A(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
//      resid(static_cast<arma::uword>(ii)) += x * diff_old_new;
//     }
//    }
//   }
//
//   // Scaled inverse-chi-square / inverse-gamma style update.
//   // R code: sigmaSqAlpha[j] = (sum(alpha[-1,j]^2) + 2) / rchisq(nAnno-1+2)
//   // Here sigmaSqAlpha_a and sigmaSqAlpha_b generalize that default.
//   double ss = 0.0;
//   int ncoef = 0;
//   for (int k = 0; k < nAnno; ++k) {
//    if (intercept_flat && k == 0) continue;
//    const double ak = alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j));
//    ss += ak * ak;
//    ++ncoef;
//   }
//
//   if (ncoef > 0) {
//    const double df = static_cast<double>(ncoef) + sigmaSqAlpha_a;
//    const double scale = ss + sigmaSqAlpha_b;
//    std::chi_squared_distribution<double> rchisq(df);
//    const double chi2 = std::max(rchisq(gen), 1e-300);
//    sigmaSqAlpha(static_cast<arma::uword>(j)) = std::max(scale / chi2, 1e-12);
//   }
//  }
// }
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_sbayesrc(
//   std::vector<std::vector<double>> wy,
//   std::vector<std::vector<double>> ww,
//   std::vector<double> yy,
//   std::vector<std::vector<double>> b_init,
//   std::vector<std::vector<double>> comp_init,
//   bool use_comp_init,
//   std::vector<std::vector<double>> r_init,
//   bool use_r_init,
//   bool rebuild_r_before_updateE,
//   std::string ld_prefix,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   arma::mat A,
//   arma::vec gamma,
//   arma::mat alpha_init,
//   arma::vec sigmaSqAlpha_init,
//   bool intercept_flat,
//   double sigmaSqAlpha_a,
//   double sigmaSqAlpha_b,
//   double pi_floor,
//   double nub,
//   double nue,
//   bool updateAlpha,
//   bool updateB,
//   bool updateE,
//   int alpha_update_every,
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
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: nt must be positive.");
//  }
//
//  const int m = static_cast<int>(wy[0].size());
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: m must be positive.");
//  }
//
//  if (nit <= 0 || nburn < 0 || nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: invalid nit/nburn/nthin.");
//  }
//
//  if (alpha_update_every <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_update_every must be positive.");
//  }
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt ||
//      (int)yy.size() != nt || (int)n.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: inconsistent trait dimensions.");
//  }
//
//  if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: priors must be nt x nt.");
//  }
//
//  if ((int)A.n_rows != m) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: A must have m rows.");
//  }
//
//  const int nAnno = static_cast<int>(A.n_cols);
//  if (nAnno <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: A must have at least one column.");
//  }
//
//  const int Kgamma = static_cast<int>(gamma.n_elem);
//  if (Kgamma < 2) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma must have at least two components, including zero.");
//  }
//
//  if (!std::isfinite(gamma(0)) || gamma(0) != 0.0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma[0] must be exactly 0.0.");
//  }
//
//  for (int k = 1; k < Kgamma; ++k) {
//   const double g = gamma(static_cast<arma::uword>(k));
//   if (!std::isfinite(g) || g <= 0.0) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma[k] must be positive for k > 0.");
//   }
//  }
//
//  const int nstep = Kgamma - 1;
//
//  if ((int)alpha_init.n_rows != nAnno || (int)alpha_init.n_cols != nstep) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_init must be ncol(A) x (length(gamma)-1).");
//  }
//
//  if ((int)sigmaSqAlpha_init.n_elem != nstep) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: sigmaSqAlpha_init must have length length(gamma)-1.");
//  }
//
//  if (!std::isfinite(pi_floor) || pi_floor <= 0.0 || pi_floor >= 1.0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: pi_floor must be in (0,1).");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m ||
//       (int)ww[t].size() != m ||
//       (int)b_init[t].size() != m) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: inconsistent marker dimensions.");
//   }
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: E must be nt x nt.");
//  }
//
//  if (use_r_init) {
//   if (static_cast<int>(r_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: r_init must have length nt when use_r_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(r_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: each r_init[t] must have length m.");
//    }
//   }
//  }
//
//  if (use_comp_init) {
//   if (static_cast<int>(comp_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: comp_init must have length nt when enabled.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(comp_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: comp_init[t] must have length m.");
//    }
//   }
//  }
//
//  arma::mat wy_mat(nt, m, arma::fill::zeros);
//  arma::mat ww_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::mat r_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> comp_mat(nt, m, arma::fill::zeros);
//
//  arma::vec yy_vec(nt, arma::fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   yy_vec(static_cast<arma::uword>(t)) = yy[t];
//
//   for (int i = 0; i < m; ++i) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword iu = static_cast<arma::uword>(i);
//
//    wy_mat(tu, iu) = wy[t][i];
//    ww_mat(tu, iu) = ww[t][i];
//    b_mat(tu, iu)  = b_init[t][i];
//   }
//
//   if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
//    sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
//   }
//  }
//
//  for (int t = 1; t < nt; ++t) {
//   if (n[t] != n[0]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_csr_sbayesrc: current shared-LD scaling assumes equal n across traits."
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
//       "stblr_cpg_omp_csr_sbayesrc: ww contains invalid value before LD pre-scaling."
//     );
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr_sbayesrc: ww differs across traits; pre-scaled shared ST LD is invalid."
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
//    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: ww contains invalid value in trait 0.");
//   }
//   xx[static_cast<std::size_t>(i)] = wi;
//  }
//
//  STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
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
//            [&](int a, int b) {
//             return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
//            });
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
//  arma::vec final_pi_active(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  std::vector<arma::mat> alpha_mean(static_cast<std::size_t>(nt));
//  std::vector<arma::vec> sigmaSqAlpha_mean(static_cast<std::size_t>(nt));
//  for (int t = 0; t < nt; ++t) {
//   alpha_mean[static_cast<std::size_t>(t)] = arma::mat(nAnno, nstep, arma::fill::zeros);
//   sigmaSqAlpha_mean[static_cast<std::size_t>(t)] = arma::vec(nstep, arma::fill::zeros);
//  }
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
//  nthreads = std::max(1, std::min(ncores, nt));
//  omp_set_num_threads(nthreads);
//
//  Rcpp::Rcout
//  << "STBLR SBayesRC CSR OpenMP requested threads = "
//  << nthreads
//  << ", omp_get_max_threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
//  Rcpp::Rcout
//  << "STBLR real-SBayesRC CSR: m=" << m
//  << ", nt=" << nt
//  << ", annotations=" << nAnno
//  << ", components=" << Kgamma
//  << ", updateAlpha=" << updateAlpha
//  << ", alpha_update_every=" << alpha_update_every
//  << "\n";
//
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(nthreads) schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
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
//    arma::Row<int> comp_t(m, arma::fill::zeros);
//
//    if (use_comp_init) {
//     for (int i = 0; i < m; ++i) {
//      const int k = static_cast<int>(std::round(comp_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)]));
//      if (k < 0 || k >= Kgamma) {
//       throw std::runtime_error("comp_init contains component outside 0..Kgamma-1.");
//      }
//      comp_t(static_cast<arma::uword>(i)) = k;
//     }
//    } else {
//     for (int i = 0; i < m; ++i) {
//      comp_t(static_cast<arma::uword>(i)) =
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
//      throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: r_init contains NaN/Inf.");
//     }
//    } else {
//     rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
//    }
//
//    double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double vg_t = 0.0;
//    double vei_t = ve_t + adjE * vg_t;
//
//    arma::mat alpha_t = alpha_init;
//    arma::vec sigmaSqAlpha_t = sigmaSqAlpha_init;
//
//    for (int j = 0; j < nstep; ++j) {
//     if (!std::isfinite(sigmaSqAlpha_t(static_cast<arma::uword>(j))) ||
//         sigmaSqAlpha_t(static_cast<arma::uword>(j)) <= 0.0) {
//      throw std::runtime_error("sigmaSqAlpha_init contains invalid value.");
//     }
//    }
//
//    arma::mat snpPi_t = compute_snp_pi_from_alpha(A, alpha_t, pi_floor);
//
//    arma::rowvec bm_t(m, arma::fill::zeros);
//    arma::rowvec dm_t(m, arma::fill::zeros);
//    arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
//
//    arma::mat alpha_accum(nAnno, nstep, arma::fill::zeros);
//    arma::vec sigmaSqAlpha_accum(nstep, arma::fill::zeros);
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//     for (int isort = 0; isort < m; ++isort) {
//      const int i = order[static_cast<std::size_t>(isort)];
//
//      sampleBeta_SBayesRC_ST_csr(
//       i,
//       snpPi_t.row(static_cast<arma::uword>(i)),
//       gamma,
//       vb_t,
//       vei_t,
//       ww_t,
//       r_t,
//       b_t,
//       comp_t,
//       ld,
//       gen_t
//      );
//     }
//
//     if (updateAlpha && ((it + 1) % alpha_update_every == 0)) {
//      update_alpha_sbayesrc(
//       A,
//       comp_t,
//       alpha_t,
//       sigmaSqAlpha_t,
//       intercept_flat,
//       sigmaSqAlpha_a,
//       sigmaSqAlpha_b,
//       gen_t
//      );
//
//      snpPi_t = compute_snp_pi_from_alpha(A, alpha_t, pi_floor);
//     }
//
//     if (updateB) {
//      sampleB_SBayesRC_ST_csr(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       comp_t,
//       gamma,
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
//       rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
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
//     vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
//
//     if (!std::isfinite(vg_t)) {
//      throw std::runtime_error("vg became NaN/Inf after computeG. iter=" + std::to_string(it));
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
//     double pi_active = 0.0;
//     for (int i = 0; i < m; ++i) {
//      double row_active = 0.0;
//      for (int k = 1; k < Kgamma; ++k) {
//       row_active += snpPi_t(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
//      }
//      pi_active += row_active;
//     }
//     pi_active /= static_cast<double>(m);
//
//     vbs_t(static_cast<arma::uword>(it)) = vb_t;
//     ves_t(static_cast<arma::uword>(it)) = ve_t;
//     vgs_t(static_cast<arma::uword>(it)) = vg_t;
//     pis_t(static_cast<arma::uword>(it)) = pi_active;
//
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       const arma::uword iu = static_cast<arma::uword>(i);
//       bm_t(iu) += b_t(iu);
//       dm_t(iu) += (comp_t(iu) > 0) ? 1.0 : 0.0;
//      }
//
//      alpha_accum += alpha_t;
//      sigmaSqAlpha_accum += sigmaSqAlpha_t;
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    bm_t /= nsamples_t;
//    dm_t /= nsamples_t;
//    alpha_accum /= nsamples_t;
//    sigmaSqAlpha_accum /= nsamples_t;
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
//    comp_mat.row(static_cast<arma::uword>(t)) = comp_t;
//
//    vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
//    vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
//    ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
//    pis_mat.row(static_cast<arma::uword>(t)) = pis_t;
//
//    final_vb(static_cast<arma::uword>(t)) = vb_t;
//    final_ve(static_cast<arma::uword>(t)) = ve_t;
//    final_vg(static_cast<arma::uword>(t)) = vg_t;
//    final_pi_active(static_cast<arma::uword>(t)) = pis_t(static_cast<arma::uword>(nit + nburn - 1));
//    nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;
//
//    alpha_mean[static_cast<std::size_t>(t)] = alpha_accum;
//    sigmaSqAlpha_mean[static_cast<std::size_t>(t)] = sigmaSqAlpha_accum;
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
//      "stblr_cpg_omp_csr_sbayesrc failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[static_cast<std::size_t>(t)]
//    );
//   }
//  }
//
//  std::vector<std::vector<std::vector<double>>> result(20);
//
//  for (int k = 0; k < 20; ++k) {
//   result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[1][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[2][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[3][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[4][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[5][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//   result[6][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
//
//   result[7][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
//   result[8][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
//   result[9][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
//
//   result[10][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[11][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[12][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[13][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[14][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//   result[15][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
//
//   result[16][static_cast<std::size_t>(t)].resize(2);
//   result[17][static_cast<std::size_t>(t)].resize(2);
//
//   result[18][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nAnno * nstep));
//   result[19][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nstep));
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
//    result[5][ts][is] = static_cast<double>(comp_mat(tu, iu));
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
//   result[16][ts][0] = 1.0 - final_pi_active(static_cast<arma::uword>(t));
//   result[16][ts][1] = final_pi_active(static_cast<arma::uword>(t));
//
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi_active(static_cast<arma::uword>(t));
//
//   result[17][ts][0] = 1.0 - mean_pi;
//   result[17][ts][1] = mean_pi;
//
//   // Flatten nAnno x nstep alpha matrix in column-major order, matching R matrix().
//   for (int j = 0; j < nstep; ++j) {
//    for (int a = 0; a < nAnno; ++a) {
//     const std::size_t pos = static_cast<std::size_t>(j * nAnno + a);
//     result[18][ts][pos] = alpha_mean[ts](static_cast<arma::uword>(a), static_cast<arma::uword>(j));
//    }
//   }
//
//   for (int j = 0; j < nstep; ++j) {
//    result[19][ts][static_cast<std::size_t>(j)] =
//     sigmaSqAlpha_mean[ts](static_cast<arma::uword>(j));
//   }
//  }
//
//  return result;
// }
