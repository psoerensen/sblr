// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "st_chain_utils.h"
#include "st_bayesrc_annotation_prior.h"
#include "st_bayesrc_annotation_selection.h"
#include "st_bayesrc_pairwise_allocation.h"
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
// intercept coefficient uses the resolved proper prior (or explicit legacy mode).
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
//   37 maf_effect_s trace, optional when estimate_maf_effect_s = true
//   38 maf_effect_s acceptance, optional when estimate_maf_effect_s = true
//   39 chain maf_effect_s trace, optional when estimate_maf_effect_s and keep_chains
//   40 chain maf_effect_s acceptance, optional when estimate_maf_effect_s and keep_chains
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
 const double score = op.corrected_rhs(i, b(iu), r);

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
  op.apply_difference(i, diff, r);
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
 const double score = op.corrected_rhs(i, b(iu), r);

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
  op.apply_difference(i, diff, r);
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
  op.apply_difference(i, diff, r);
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

   // vb is the global variance; fixed maf_effect_s enters the sufficient
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

inline void fill_maf_effect_s_prior_scale_sbayesrc(
  int m,
  double maf_effect_s,
  const arma::rowvec& log_h,
  arma::rowvec& prior_scale
) {
 if (static_cast<int>(prior_scale.n_elem) != m) {
  prior_scale.set_size(static_cast<arma::uword>(m));
 }

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double scale_i = std::exp((maf_effect_s + 1.0) * log_h(iu));
  if (!std::isfinite(scale_i) || scale_i <= 0.0) {
   throw std::runtime_error("fill_maf_effect_s_prior_scale_sbayesrc: invalid dynamic prior scale.");
  }
  prior_scale(iu) = scale_i;
 }
}

inline double logpost_maf_effect_s_sbayesrc(
  double maf_effect_s,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  double vb,
  const arma::vec& gamma,
  const arma::rowvec& log_h,
  double prior_lower,
  double prior_upper
) {
 if (!std::isfinite(maf_effect_s) ||
     maf_effect_s < prior_lower ||
     maf_effect_s > prior_upper) {
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

  const double log_q = (maf_effect_s + 1.0) * log_h(iu);
  const double q = std::exp(log_q);
  if (!std::isfinite(q) || q <= 0.0) {
   return -std::numeric_limits<double>::infinity();
  }

  const double bi = b(iu);
  lp += -0.5 * (log_q + bi * bi / (vb * gk * q));
 }

 return lp;
}

inline bool update_maf_effect_s_sbayesrc(
  double& maf_effect_s,
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
 const double prop = maf_effect_s + norm(gen);
 if (!std::isfinite(prop) || prop < prior_lower || prop > prior_upper) {
  return false;
 }

 const double lp_current = logpost_maf_effect_s_sbayesrc(
  maf_effect_s, b, comp, vb, gamma, log_h, prior_lower, prior_upper
 );
 const double lp_prop = logpost_maf_effect_s_sbayesrc(
  prop, b, comp, vb, gamma, log_h, prior_lower, prior_upper
 );
 const double log_alpha = lp_prop - lp_current;
 if (!std::isfinite(log_alpha)) return false;

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 if (std::log(std::max(runif(gen), 1e-300)) < log_alpha) {
  maf_effect_s = prop;
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
 bool estimate_maf_effect_s;
 bool use_maf_effect_s_prior_scale;
 bool uses_retained_low_rank;
 bool selection_enabled;
 const std::vector<int>& convergence_markers;
 const std::vector<double>& phenotype_variance;
 sblr::core::BlockResidualControl block_residual_control;
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
  arma::mat intercept_prior_resolved,
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
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale,
  bool estimate_maf_effect_s,
  double maf_effect_s_init,
  Rcpp::NumericVector maf_effect_s_prior,
  double maf_effect_s_proposal_sd,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h,
  Rcpp::IntegerVector convergence_markers,
  bool convergence_annotations,
  bool convergence_b,
  bool convergence_d,
  bool convergence_component,
  int low_rank_residual_rebuild_every,
  std::vector<double> phenotype_variance,
  std::string residual_policy,
  std::string block_ve_mode,
  double resam_thresh,
  double minimum_ve_ratio,
  bool block_ve_keep_history,
  StBayesRCSelectionGenomicConfig selection_config,
  MakeOperator make_operator
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: nt must be positive.");
 }

 const int m = static_cast<int>(wy[0].size());
 const std::vector<int> convergence_markers_cpp=
  Rcpp::as<std::vector<int>>(convergence_markers);

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

 bool use_maf_effect_s_prior_scale = maf_effect_s_prior_scale.isNotNull();
 Rcpp::NumericVector maf_effect_s_prior_scale_vec;
 if (use_maf_effect_s_prior_scale) {
  maf_effect_s_prior_scale_vec = Rcpp::NumericVector(maf_effect_s_prior_scale);
  use_maf_effect_s_prior_scale = maf_effect_s_prior_scale_vec.size() > 0;
 }

 if (use_maf_effect_s_prior_scale &&
     static_cast<int>(maf_effect_s_prior_scale_vec.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_prior_scale must have length m.");
 }

 if (estimate_maf_effect_s && use_maf_effect_s_prior_scale) {
  throw std::runtime_error(
   "stblr_cpg_omp_csr_sbayesrc: fixed maf_effect_s_prior_scale and estimate_maf_effect_s cannot both be used."
  );
 }

 if (estimate_maf_effect_s) {
  if (!std::isfinite(maf_effect_s_init)) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_init must be finite.");
  }
  if (maf_effect_s_prior.size() != 2 ||
      !std::isfinite(maf_effect_s_prior[0]) ||
      !std::isfinite(maf_effect_s_prior[1])) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_prior must be finite length 2.");
  }
  if (maf_effect_s_prior[0] >= maf_effect_s_prior[1]) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_prior lower bound must be less than upper bound.");
  }
  if (maf_effect_s_init < maf_effect_s_prior[0] ||
      maf_effect_s_init > maf_effect_s_prior[1]) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_init must lie within maf_effect_s_prior.");
  }
  if (!std::isfinite(maf_effect_s_proposal_sd) || maf_effect_s_proposal_sd <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_proposal_sd must be positive finite.");
  }
  if (maf_effect_s_log_h.isNull()) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_log_h is required when estimate_maf_effect_s = true.");
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
 if (use_maf_effect_s_prior_scale) {
  prior_scale.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double scale_i = maf_effect_s_prior_scale_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(scale_i) || scale_i <= 0.0) {
    throw std::runtime_error(
     "stblr_cpg_omp_csr_sbayesrc: maf_effect_s_prior_scale must contain positive finite values."
    );
   }
   prior_scale(static_cast<arma::uword>(i)) = scale_i;
  }
 }

 arma::rowvec maf_effect_s_log_h_row;
 if (estimate_maf_effect_s) {
  Rcpp::NumericVector maf_effect_s_log_h_vec(maf_effect_s_log_h);
  if (static_cast<int>(maf_effect_s_log_h_vec.size()) != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_log_h must have length m.");
  }
  maf_effect_s_log_h_row.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double log_h_i = maf_effect_s_log_h_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(log_h_i)) {
    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: maf_effect_s_log_h must contain finite values.");
   }
   maf_effect_s_log_h_row(static_cast<arma::uword>(i)) = log_h_i;
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
 sblr::core::BlockResidualControl block_residual_control;
 block_residual_control.policy =
  sblr::core::parse_block_residual_policy(residual_policy);
 block_residual_control.mode = sblr::core::parse_block_ve_mode(block_ve_mode);
 block_residual_control.resam_threshold = resam_thresh;
 block_residual_control.minimum_ve_ratio = minimum_ve_ratio;
 block_residual_control.keep_history = block_ve_keep_history;
 sblr::core::validate_block_residual_control(block_residual_control);
 if (block_residual_control.uses_block_variance() && !op.uses_retained_low_rank())
  throw std::runtime_error(
   "block residual policies require representation = 'low_rank'.");
 if (low_rank_residual_rebuild_every < 0) {
  throw std::runtime_error(
   "low_rank_residual_rebuild_every must be non-negative.");
 }
 if (op.uses_retained_low_rank() && rebuild_r_before_updateE) {
  throw std::runtime_error(
   "rebuild_r_before_updateE is incompatible with retained low rank; use "
   "low_rank_residual_rebuild_every instead.");
 }
 if (!op.uses_retained_low_rank() &&
     low_rank_residual_rebuild_every != 0) {
  throw std::runtime_error(
   "low_rank_residual_rebuild_every is only supported by retained low rank.");
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
 const StBayesRCInterceptPrior intercept_prior =
  st_bayesrc_parse_intercept_prior(intercept_prior_resolved, nstep);
 alpha_contract.intercept_prior_legacy_flat = intercept_prior.legacy_flat;
 alpha_contract.intercept_prior_mean.assign(
  intercept_prior.mean.begin(), intercept_prior.mean.end());
 alpha_contract.intercept_prior_precision.assign(
  intercept_prior.precision.begin(), intercept_prior.precision.end());
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
  maf_effect_s_log_h_row,
  annotation_contract,
  component_contract,
  alpha_contract,
  probability_contract,
  prior_contract,
  control_contract,
  output_contract,
  selection_config,
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
  low_rank_residual_rebuild_every,
  intercept_prior,
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
  use_maf_effect_s_prior_scale,
  estimate_maf_effect_s,
  maf_effect_s_init,
  maf_effect_s_prior[0],
  maf_effect_s_prior[1],
  maf_effect_s_proposal_sd,
  convergence_markers_cpp,
  phenotype_variance,
  block_residual_control,
  convergence_annotations,
  convergence_b,
  convergence_d,
  convergence_component
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
  estimate_maf_effect_s,
  use_maf_effect_s_prior_scale,
  op.uses_retained_low_rank(),
  selection_config.enabled,
  convergence_markers_cpp,
  phenotype_variance,
  block_residual_control
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
 const bool estimate_maf_effect_s = metadata.estimate_maf_effect_s;
 const bool use_maf_effect_s_prior_scale = metadata.use_maf_effect_s_prior_scale;
 const bool selection_enabled = metadata.selection_enabled;
 const std::vector<int>& convergence_markers=metadata.convergence_markers;
 const std::vector<double>& phenotype_variance=metadata.phenotype_variance;
 const auto& block_residual_control=metadata.block_residual_control;

 const auto& bm_task = execution_result.bm_task;
 const auto& dm_task = execution_result.dm_task;
 const auto& comp_task_double = execution_result.comp_task_double;
 const auto& vbs_task = execution_result.vbs_task;
 const auto& vgs_task = execution_result.vgs_task;
 const auto& ves_task = execution_result.ves_task;
 const auto& pis_task = execution_result.pis_task;
 const auto& vles_task = execution_result.vles_task;
 const auto& vlds_task = execution_result.vlds_task;
 const auto& maf_effect_s_task = execution_result.maf_effect_s_task;
 const auto& final_pi_component_task = execution_result.final_pi_component_task;
 const auto& maf_effect_s_attempted_task = execution_result.maf_effect_s_attempted_task;
 const auto& convergence_alpha_task=execution_result.convergence_alpha_task;
 const auto& convergence_sigma_task=execution_result.convergence_sigma_task;
 const auto& convergence_b_task=execution_result.convergence_b_task;
 const auto& convergence_d_task=execution_result.convergence_d_task;
 const auto& convergence_component_task=execution_result.convergence_component_task;
 const auto& convergence_aggregate_task=execution_result.convergence_aggregate_task;
 const auto& selection_pip_task = execution_result.selection_pip_task;
 const auto& selection_pi_a_mean_task = execution_result.selection_pi_a_mean_task;
 const auto& selection_tau2_mean_task = execution_result.selection_tau2_mean_task;
 const auto& selection_included_mean_task = execution_result.selection_included_mean_task;
 const auto& selection_switches_task = execution_result.selection_switches_task;
 const auto& selection_alpha_conditional_mean_task =
  execution_result.selection_alpha_conditional_mean_task;
 const auto& selection_delta_final_task = execution_result.selection_delta_final_task;
 const auto& maf_effect_s_accepted_task = execution_result.maf_effect_s_accepted_task;
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
 const auto& maf_effect_s_mat = execution_result.maf_effect_s_mat;
 const auto& final_vb = execution_result.final_vb;
 const auto& final_vg = execution_result.final_vg;
 const auto& final_ve = execution_result.final_ve;
 const auto& final_pi_component = execution_result.final_pi_component;
 const auto& nsamples_vec = execution_result.nsamples_vec;
 const auto& maf_effect_s_attempted_vec = execution_result.maf_effect_s_attempted_vec;
 const auto& maf_effect_s_accepted_vec = execution_result.maf_effect_s_accepted_vec;
 const auto& alpha_mean = execution_result.alpha_mean;
 const auto& sigmaSqAlpha_mean = execution_result.sigmaSqAlpha_mean;
 const auto& comp_prob_mean = execution_result.comp_prob_mean;
 const auto& ncomp_mean = execution_result.ncomp_mean;
 const auto& task_seconds = execution_result.task_seconds;
 const auto& ld_swap_diagnostics = execution_result.ld_swap_diagnostics;
 const auto& ld_swap_chain_diagnostics = execution_result.ld_swap_chain_diagnostics;
 const auto& block_ve_posterior_mean_task =
  execution_result.block_ve_posterior_mean_task;
 const auto& block_ve_final_task = execution_result.block_ve_final_task;
 const auto& block_ve_resampled_task = execution_result.block_ve_resampled_task;
 const auto& block_ve_reset_task = execution_result.block_ve_reset_task;
 const auto& summary_heritability_task =
  execution_result.summary_heritability_task;
 const auto& block_ve_history_task = execution_result.block_ve_history_task;

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
 Rcpp::NumericVector maf_effect_sd(nt);
 Rcpp::NumericVector selection_min(nt);
 Rcpp::NumericVector selection_max(nt);
 Rcpp::NumericVector selection_acceptance(nt);
 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
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
 if (block_residual_control.uses_block_variance()) {
  arma::mat summary_h2(nt, n_trace, arma::fill::zeros);
  for (int t = 0; t < nt; ++t)
   for (int chain = 0; chain < nchains; ++chain)
    summary_h2.row(static_cast<arma::uword>(t)) +=
     summary_heritability_task.row(
      static_cast<arma::uword>(t * nchains + chain));
  summary_h2 /= static_cast<double>(nchains);
  trace["heritability_summary"] = trace_matrix(summary_h2);
 }

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
 if (metadata.uses_retained_low_rank) {
  diagnostics["low_rank_residual"] = Rcpp::List::create(
   Rcpp::Named("low_rank_residual_rebuild_every") =
    execution_result.low_rank_residual_diagnostics.col(0),
   Rcpp::Named("low_rank_residual_rebuild_count") =
    execution_result.low_rank_residual_diagnostics.col(1),
   Rcpp::Named("low_rank_residual_max_abs_drift") =
    execution_result.low_rank_residual_diagnostics.col(2)
  );
 }

 if (block_residual_control.uses_block_variance()) {
  Rcpp::RObject history = R_NilValue;
  if (block_residual_control.keep_history) {
   Rcpp::List history_list(block_ve_history_task.size());
   for (std::size_t task = 0; task < block_ve_history_task.size(); ++task)
    history_list[task] = block_ve_history_task[task];
   history = history_list;
  }
  diagnostics["block_residual"] = Rcpp::List::create(
   Rcpp::Named("residual_policy") =
    sblr::core::block_residual_policy_name(block_residual_control.policy),
   Rcpp::Named("block_ve_mode") =
    sblr::core::block_ve_mode_name(block_residual_control.mode),
   Rcpp::Named("resam_thresh") = block_residual_control.resam_threshold,
   Rcpp::Named("minimum_ve_ratio") = block_residual_control.minimum_ve_ratio,
   Rcpp::Named("phenotype_variance") = phenotype_variance,
   Rcpp::Named("posterior_mean_per_chain_block") =
    block_ve_posterior_mean_task,
   Rcpp::Named("final_per_chain_block") = block_ve_final_task,
   Rcpp::Named("resampled_per_chain_block") = block_ve_resampled_task,
   Rcpp::Named("minimum_ratio_resets_per_chain_block") = block_ve_reset_task,
   Rcpp::Named("history") = history
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
    if (estimate_maf_effect_s) {
     const double attempted = maf_effect_s_attempted_task(task_u);
     const double accepted = maf_effect_s_accepted_task(task_u);
     chain_selection["trace"] = maf_effect_s_task.row(task_u).t();
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
      Rcpp::Named("sigmaSqAlpha") = sigmaSqAlpha_mean_task[static_cast<std::size_t>(task)],
      Rcpp::Named("annotation_pip") = selection_enabled ?
       Rcpp::wrap(selection_pip_task[static_cast<std::size_t>(task)]) : R_NilValue,
      Rcpp::Named("annotation_delta_final") = selection_enabled ?
       Rcpp::wrap(selection_delta_final_task[static_cast<std::size_t>(task)]) : R_NilValue,
      Rcpp::Named("alpha_mean_given_inclusion") = selection_enabled ?
       Rcpp::wrap(selection_alpha_conditional_mean_task[static_cast<std::size_t>(task)]) : R_NilValue,
      Rcpp::Named("annotation_pi_A") = selection_enabled ?
       Rcpp::wrap(selection_pi_a_mean_task[static_cast<std::size_t>(task)]) : R_NilValue,
      Rcpp::Named("annotation_tau2") = selection_enabled ?
       Rcpp::wrap(selection_tau2_mean_task[static_cast<std::size_t>(task)]) : R_NilValue,
      Rcpp::Named("annotation_included_mean") = selection_enabled ?
       Rcpp::wrap(selection_included_mean_task[static_cast<std::size_t>(task)]) : R_NilValue,
      Rcpp::Named("annotation_switches") = selection_enabled ?
       Rcpp::wrap(selection_switches_task[static_cast<std::size_t>(task)]) : R_NilValue
     ),
     Rcpp::Named("convergence_trace")=Rcpp::List::create(
      Rcpp::Named("marker_index")=Rcpp::wrap(convergence_markers),
      Rcpp::Named("b")=convergence_b_task[static_cast<std::size_t>(task)],
      Rcpp::Named("d")=convergence_d_task[static_cast<std::size_t>(task)],
      Rcpp::Named("component")=convergence_component_task[static_cast<std::size_t>(task)],
      Rcpp::Named("component_count")=convergence_aggregate_task[static_cast<std::size_t>(task)].component_count,
      Rcpp::Named("realized_active_count")=convergence_aggregate_task[static_cast<std::size_t>(task)].realized_active_count,
      Rcpp::Named("stick_eligible_count")=convergence_aggregate_task[static_cast<std::size_t>(task)].stick_eligible_count,
      Rcpp::Named("stick_continue_count")=convergence_aggregate_task[static_cast<std::size_t>(task)].stick_continue_count,
      Rcpp::Named("stick_stop_count")=convergence_aggregate_task[static_cast<std::size_t>(task)].stick_stop_count,
      Rcpp::Named("alpha")=convergence_alpha_task[static_cast<std::size_t>(task)],
      Rcpp::Named("sigmaSqAlpha")=convergence_sigma_task[static_cast<std::size_t>(task)]),
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
   Rcpp::Named("model") = selection_enabled ? "sbayesrc_selection" : "sbayesrc",
   Rcpp::Named("backend") = selection_enabled ?
    "csr_sbayesrc_selection_internal" : "csr_sbayesrc",
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
   Rcpp::Named("n_groups") = 0,
   Rcpp::Named("residual_policy") =
    sblr::core::block_residual_policy_name(block_residual_control.policy),
   Rcpp::Named("block_ve_mode") =
    sblr::core::block_ve_mode_name(block_residual_control.mode),
   Rcpp::Named("block_ve_definition") =
    block_residual_control.uses_block_variance() ?
     "mean_block_residual_variance" : "global_projected_residual_variance",
   Rcpp::Named("heritability_definition") =
    block_residual_control.uses_block_variance() ?
     "sum_block_genetic_variance_over_phenotype_variance" :
     "genetic_over_genetic_plus_residual"
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

// Internal deterministic regression hook for the shared BED/CSR BayesRC
// latent-variable kernel.
// [[Rcpp::export(name = ".st_bayesrc_truncated_normal_draws")]]
Rcpp::NumericMatrix st_bayesrc_truncated_normal_draws_test(
 Rcpp::NumericVector location,
 int draws,
 int seed
) {
 if (location.size() == 0 || draws <= 0)
  throw std::runtime_error("location and draws must be non-empty and positive.");
 std::mt19937 generator(static_cast<unsigned int>(seed));
 Rcpp::NumericMatrix result(draws, 2 * location.size());
 for (R_xlen_t j = 0; j < location.size(); ++j) {
  for (int i = 0; i < draws; ++i) {
   result(i, j) = st_bayesrc_sample_truncated_normal_std(
    location[j], true, generator);
   result(i, j + location.size()) = st_bayesrc_sample_truncated_normal_std(
    location[j], false, generator);
  }
 }
 return result;
}

// Deterministic internal test hook for the shared BED/CSR/block-eigen/MT
// annotation kernel. It is deliberately not exported from the R namespace.
// [[Rcpp::export(name = ".st_bayesrc_annotation_update")]]
Rcpp::List st_bayesrc_annotation_update_test(
  arma::mat annotation,
  arma::rowvec component_numeric,
  arma::mat alpha,
  arma::vec sigma_sq_alpha,
  arma::mat intercept_prior_resolved,
  double sigma_alpha_a,
  double sigma_alpha_b,
  int seed
) {
 arma::Row<int> component(component_numeric.n_elem);
 for (arma::uword i = 0; i < component_numeric.n_elem; ++i) {
  const double value = component_numeric(i);
  if (!std::isfinite(value) || value < 0.0 || value != std::floor(value)) {
   throw std::invalid_argument("component must contain non-negative integer indices.");
  }
  component(i) = static_cast<int>(value);
 }
 const int nstep = static_cast<int>(alpha.n_cols);
 const auto prior = st_bayesrc_parse_intercept_prior(
  intercept_prior_resolved, nstep);
 std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));
 const auto diagnostics = st_bayesrc_update_annotation_prior(
  annotation, component, alpha, sigma_sq_alpha, prior,
  sigma_alpha_a, sigma_alpha_b, gen);
 return Rcpp::List::create(
  Rcpp::Named("alpha") = alpha,
  Rcpp::Named("sigmaSqAlpha") = sigma_sq_alpha,
  Rcpp::Named("eligible") = diagnostics.eligible,
  Rcpp::Named("continuation") = diagnostics.continuation,
  Rcpp::Named("prior_only") = diagnostics.prior_only,
  Rcpp::Named("px_attempted") = diagnostics.px_attempted,
  Rcpp::Named("px_accepted") = diagnostics.px_accepted,
  Rcpp::Named("px_abs_log_scale") = diagnostics.px_abs_log_scale,
  Rcpp::Named("px_alpha_jump") = diagnostics.px_alpha_jump
 );
}

// Deterministic development hook for alpha conditional moments with fixed
// Albert--Chib latent responses. It does not mutate a production chain.
// [[Rcpp::export(name = ".st_bayesrc_alpha_conditional_moments")]]
arma::mat st_bayesrc_alpha_conditional_moments_test(
 arma::mat annotation,
 arma::vec latent,
 arma::vec alpha,
 double sigma_sq_alpha,
 double intercept_mean,
 double intercept_sd
) {
 if (annotation.n_rows != latent.n_elem || annotation.n_cols != alpha.n_elem ||
     !annotation.is_finite() || !latent.is_finite() || !alpha.is_finite() ||
     !std::isfinite(sigma_sq_alpha) || sigma_sq_alpha <= 0.0 ||
     !std::isfinite(intercept_mean) || !std::isfinite(intercept_sd) ||
     intercept_sd <= 0.0) {
  throw std::invalid_argument("invalid fixed-latent alpha conditional inputs.");
 }
 arma::vec residual = latent - annotation * alpha;
 arma::mat result(annotation.n_cols, 2, arma::fill::zeros);
 for (arma::uword k = 0; k < annotation.n_cols; ++k) {
  const arma::vec x = annotation.col(k);
  const double diagonal = arma::dot(x, x);
  const double rhs = arma::dot(x, residual) + diagonal * alpha(k);
  const double prior_mean = k == 0 ? intercept_mean : 0.0;
  const double prior_precision = k == 0 ?
   1.0 / (intercept_sd * intercept_sd) : 1.0 / sigma_sq_alpha;
  const auto conditional = st_bayesrc_scalar_conditional(
   diagonal, rhs, prior_mean, prior_precision);
  result(k, 0) = conditional.mean;
  result(k, 1) = conditional.variance;
 }
 return result;
}

// Stochastic development hook for the exact sigmaSqAlpha conditional.
// [[Rcpp::export(name = ".st_bayesrc_sigma_sq_alpha_draws")]]
Rcpp::NumericVector st_bayesrc_sigma_sq_alpha_draws_test(
 arma::vec alpha_non_intercept,
 double prior_a,
 double prior_b,
 int draws,
 int seed
) {
 if (!alpha_non_intercept.is_finite() || draws <= 0)
  throw std::invalid_argument("invalid sigmaSqAlpha draw inputs.");
 std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));
 Rcpp::NumericVector result(draws);
 const double ss = arma::dot(alpha_non_intercept, alpha_non_intercept);
 for (int i = 0; i < draws; ++i) {
  result[i] = st_bayesrc_sample_sigma_sq_alpha(
   ss, static_cast<int>(alpha_non_intercept.n_elem), prior_a, prior_b, gen);
 }
 return result;
}

// Development-only hierarchy sampler with allocations fixed. This isolates
// the shared probit-stick update from every likelihood and marker-effect
// transition and returns compact chain traces.
// [[Rcpp::export(name = ".st_bayesrc_frozen_hierarchy_chains")]]
Rcpp::List st_bayesrc_frozen_hierarchy_chains_test(
 arma::mat annotation,
 arma::rowvec component_numeric,
 arma::mat alpha_initial,
 arma::vec sigma_initial,
 arma::mat intercept_prior_resolved,
 double sigma_alpha_a,
 double sigma_alpha_b,
 double probability_floor,
 int iterations,
 Rcpp::IntegerVector chain_seeds,
 int cores
) {
 if (iterations <= 0 || chain_seeds.size() <= 0 || cores <= 0 ||
     annotation.n_rows != component_numeric.n_elem ||
     alpha_initial.n_rows != annotation.n_cols ||
     sigma_initial.n_elem != alpha_initial.n_cols ||
     !annotation.is_finite() || !alpha_initial.is_finite() ||
     !sigma_initial.is_finite() || arma::any(sigma_initial <= 0.0)) {
  throw std::invalid_argument("invalid frozen alpha-hierarchy inputs.");
 }
 arma::Row<int> component(component_numeric.n_elem);
 for (arma::uword i = 0; i < component_numeric.n_elem; ++i) {
  const double value = component_numeric(i);
  if (!std::isfinite(value) || value < 0.0 || value != std::floor(value) ||
      value > static_cast<double>(alpha_initial.n_cols)) {
   throw std::invalid_argument("component contains an invalid state.");
  }
  component(i) = static_cast<int>(value);
 }
 const int chain_count = chain_seeds.size();
 const std::vector<int> chain_seed_values =
  Rcpp::as<std::vector<int>>(chain_seeds);
 const int coefficient_count = static_cast<int>(alpha_initial.n_elem);
 const int step_count = static_cast<int>(alpha_initial.n_cols);
 const int component_count = step_count + 1;
 std::vector<arma::mat> alpha_trace(static_cast<std::size_t>(chain_count));
 std::vector<arma::mat> sigma_trace(static_cast<std::size_t>(chain_count));
 std::vector<arma::mat> probability_trace(static_cast<std::size_t>(chain_count));
 std::vector<std::string> error(static_cast<std::size_t>(chain_count));

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
 for (int chain = 0; chain < chain_count; ++chain) {
  try {
   arma::mat alpha = alpha_initial;
   arma::vec sigma = sigma_initial;
   const auto prior = st_bayesrc_parse_intercept_prior(
    intercept_prior_resolved, step_count);
   std::mt19937 generator(static_cast<std::mt19937::result_type>(
    chain_seed_values[static_cast<std::size_t>(chain)]));
   arma::mat alpha_out(iterations, coefficient_count, arma::fill::zeros);
   arma::mat sigma_out(iterations, step_count, arma::fill::zeros);
   arma::mat probability_out(
    iterations, 1 + 3 * component_count, arma::fill::zeros);
   for (int iteration = 0; iteration < iterations; ++iteration) {
    st_bayesrc_update_annotation_prior(
     annotation, component, alpha, sigma, prior,
     sigma_alpha_a, sigma_alpha_b, generator);
    alpha_out.row(static_cast<arma::uword>(iteration)) =
     arma::vectorise(alpha).t();
    sigma_out.row(static_cast<arma::uword>(iteration)) = sigma.t();
    const arma::mat marker_probability = st_bayesrc_compute_snp_pi(
     annotation, alpha, probability_floor);
    probability_out(iteration, 0) =
     arma::accu(1.0 - marker_probability.col(0));
    for (int k = 0; k < component_count; ++k) {
     const arma::vec column = marker_probability.col(static_cast<arma::uword>(k));
     probability_out(iteration, 1 + k) = arma::mean(column);
     probability_out(iteration, 1 + component_count + k) = column.min();
     probability_out(iteration, 1 + 2 * component_count + k) = column.max();
    }
   }
   alpha_trace[static_cast<std::size_t>(chain)] = std::move(alpha_out);
   sigma_trace[static_cast<std::size_t>(chain)] = std::move(sigma_out);
   probability_trace[static_cast<std::size_t>(chain)] =
    std::move(probability_out);
  } catch (const std::exception& exception) {
   error[static_cast<std::size_t>(chain)] = exception.what();
  }
 }
 for (const std::string& message : error) {
  if (!message.empty()) throw std::runtime_error(message);
 }
 Rcpp::List chains(chain_count);
 for (int chain = 0; chain < chain_count; ++chain) {
  chains[chain] = Rcpp::List::create(
   Rcpp::Named("alpha") = alpha_trace[static_cast<std::size_t>(chain)],
   Rcpp::Named("sigmaSqAlpha") = sigma_trace[static_cast<std::size_t>(chain)],
   Rcpp::Named("probability_summary") =
    probability_trace[static_cast<std::size_t>(chain)]);
 }
 return chains;
}

static std::vector<arma::uvec> st_bayesrc_selection_rows_from_r(
 Rcpp::List eligible,
 arma::uword marker_count
) {
 std::vector<arma::uvec> result(static_cast<std::size_t>(eligible.size()));
 for (R_xlen_t stick = 0; stick < eligible.size(); ++stick) {
  Rcpp::IntegerVector input = eligible[stick];
  arma::uvec rows(static_cast<arma::uword>(input.size()));
  for (R_xlen_t i = 0; i < input.size(); ++i) {
   if (input[i] == NA_INTEGER || input[i] <= 0 ||
       input[i] > static_cast<int>(marker_count)) {
    throw std::invalid_argument("SBayesRC-S eligible rows must be one-based valid indices.");
   }
   rows(static_cast<arma::uword>(i)) = static_cast<arma::uword>(input[i] - 1);
  }
  result[static_cast<std::size_t>(stick)] = std::move(rows);
 }
 return result;
}

static std::vector<arma::ivec> st_bayesrc_selection_outcome_from_r(
 Rcpp::List outcome
) {
 std::vector<arma::ivec> result(static_cast<std::size_t>(outcome.size()));
 for (R_xlen_t stick = 0; stick < outcome.size(); ++stick) {
  Rcpp::IntegerVector input = outcome[stick];
  arma::ivec values(static_cast<arma::uword>(input.size()));
  for (R_xlen_t i = 0; i < input.size(); ++i) {
   if (input[i] == NA_INTEGER || (input[i] != 0 && input[i] != 1)) {
    throw std::invalid_argument("SBayesRC-S outcomes must be zero or one.");
   }
   values(static_cast<arma::uword>(i)) = input[i];
  }
  result[static_cast<std::size_t>(stick)] = std::move(values);
 }
 return result;
}

static arma::uvec st_bayesrc_selection_delta_from_r(
 Rcpp::IntegerVector delta,
 arma::uword annotation_count
) {
 if (delta.size() != static_cast<R_xlen_t>(annotation_count)) {
  throw std::invalid_argument("SBayesRC-S delta length mismatch.");
 }
 arma::uvec result(annotation_count);
 for (R_xlen_t j = 0; j < delta.size(); ++j) {
  if (delta[j] == NA_INTEGER || (delta[j] != 0 && delta[j] != 1)) {
   throw std::invalid_argument("SBayesRC-S delta must contain zero or one.");
  }
  result(static_cast<arma::uword>(j)) = static_cast<arma::uword>(delta[j]);
 }
 return result;
}

// Deterministic internal oracle for every Phase-4A hierarchy calculation.
// It is deliberately absent from NAMESPACE and cannot run a genomic sampler.
// [[Rcpp::export(name = ".st_bayesrc_selection_math")]]
Rcpp::List st_bayesrc_selection_math_test(
 arma::mat annotation,
 Rcpp::List eligible,
 Rcpp::List latent,
 arma::mat alpha,
 Rcpp::IntegerVector delta,
 double pi_a,
 arma::vec tau2,
 double a_pi,
 double b_pi,
 double a_tau,
 double b_tau,
 double probability_floor
) {
 const auto rows = st_bayesrc_selection_rows_from_r(
  eligible, annotation.n_rows);
 if (latent.size() != eligible.size()) {
  throw std::invalid_argument("SBayesRC-S latent stick count mismatch.");
 }
 std::vector<arma::vec> latent_values(static_cast<std::size_t>(latent.size()));
 for (R_xlen_t stick = 0; stick < latent.size(); ++stick) {
  latent_values[static_cast<std::size_t>(stick)] =
   Rcpp::as<arma::vec>(latent[stick]);
 }
 StBayesRCSelectionState state{
  st_bayesrc_selection_delta_from_r(delta, annotation.n_cols),
  alpha, pi_a, tau2};
 const StBayesRCSelectionHyperParameters hyper{a_pi, b_pi, a_tau, b_tau};
 st_bayesrc_selection_validate(annotation, rows, nullptr, state, hyper);

 const arma::uword annotation_count = annotation.n_cols;
 const arma::uword stick_count = rows.size();
 arma::mat s(annotation_count, stick_count, arma::fill::zeros);
 arma::mat t(annotation_count, stick_count, arma::fill::zeros);
 arma::mat log_bf(annotation_count, stick_count, arma::fill::zeros);
 arma::mat conditional_mean(annotation_count, stick_count, arma::fill::zeros);
 arma::mat conditional_variance(annotation_count, stick_count, arma::fill::zeros);
 arma::vec inclusion_logit(annotation_count, arma::fill::zeros);
 Rcpp::List eta(stick_count);
 for (arma::uword stick = 0; stick < stick_count; ++stick) {
  arma::vec eta_stick(rows[static_cast<std::size_t>(stick)].n_elem,
                      arma::fill::zeros);
  for (arma::uword ii = 0; ii < eta_stick.n_elem; ++ii) {
   const arma::uword row = rows[static_cast<std::size_t>(stick)](ii);
   eta_stick(ii) = state.alpha(0u, stick) +
    arma::dot(annotation.row(row), state.alpha.col(stick).subvec(
     1u, annotation_count));
  }
  eta[stick] = eta_stick;
 }
 for (arma::uword j = 0; j < annotation_count; ++j) {
  double logit = std::log(pi_a) - std::log1p(-pi_a);
  for (arma::uword stick = 0; stick < stick_count; ++stick) {
   const arma::uvec& index = rows[static_cast<std::size_t>(stick)];
   arma::vec residual = latent_values[static_cast<std::size_t>(stick)] -
    state.alpha(0u, stick);
   for (arma::uword other = 0; other < annotation_count; ++other) {
    if (other != j) residual -= st_bayesrc_selection_column_rows(
     annotation, other, index) *
     state.alpha(other + 1u, stick);
   }
   const auto moments = st_bayesrc_selection_moments(
    st_bayesrc_selection_column_rows(annotation, j, index),
    residual, tau2(stick));
   s(j, stick) = moments.s;
   t(j, stick) = moments.t;
   log_bf(j, stick) = moments.log_bf;
   conditional_mean(j, stick) = moments.mean;
   conditional_variance(j, stick) = moments.variance;
   logit += moments.log_bf;
  }
  inclusion_logit(j) = logit;
 }
 const arma::vec beta_parameters = st_bayesrc_selection_beta_parameters(
  state.delta, a_pi, b_pi);
 const arma::mat ig_parameters = st_bayesrc_selection_ig_parameters(
  state.alpha, state.delta, a_tau, b_tau);
 const arma::mat q = st_bayesrc_selection_compute_q(annotation, state.alpha);
 const arma::mat component_probability = st_bayesrc_selection_compute_pi(
  annotation, state.alpha, probability_floor);
 return Rcpp::List::create(
  Rcpp::Named("eta") = eta,
  Rcpp::Named("s") = s,
  Rcpp::Named("t") = t,
  Rcpp::Named("log_bf") = log_bf,
  Rcpp::Named("inclusion_logit") = inclusion_logit,
  Rcpp::Named("inclusion_probability") = 1.0 / (1.0 + arma::exp(-inclusion_logit)),
  Rcpp::Named("conditional_mean") = conditional_mean,
  Rcpp::Named("conditional_variance") = conditional_variance,
  Rcpp::Named("beta_parameters") = beta_parameters,
  Rcpp::Named("ig_parameters") = ig_parameters,
  Rcpp::Named("q") = q,
  Rcpp::Named("component_probability") = component_probability
 );
}

// Isolated C++ implementation of the validated Phase-3 hierarchy.  This hook
// has no beta, allocation likelihood, LD, residual, or genomic variance state.
// [[Rcpp::export(name = ".st_bayesrc_selection_hierarchy")]]
Rcpp::List st_bayesrc_selection_hierarchy_test(
 arma::mat annotation,
 Rcpp::List eligible,
 Rcpp::List outcome,
 arma::mat alpha_initial,
 Rcpp::IntegerVector delta_initial,
 double pi_a_initial,
 arma::vec tau2_initial,
 double a_pi,
 double b_pi,
 double a_tau,
 double b_tau,
 double probability_floor,
 int iterations,
 int burn,
 int seed,
 Rcpp::IntegerVector fixed_delta
) {
 if (iterations <= burn || burn < 0) {
  throw std::invalid_argument("SBayesRC-S iterations must exceed burn-in.");
 }
 const auto rows = st_bayesrc_selection_rows_from_r(
  eligible, annotation.n_rows);
 const auto outcomes = st_bayesrc_selection_outcome_from_r(outcome);
 StBayesRCSelectionState state{
  st_bayesrc_selection_delta_from_r(delta_initial, annotation.n_cols),
  alpha_initial, pi_a_initial, tau2_initial};
 const StBayesRCSelectionHyperParameters hyper{a_pi, b_pi, a_tau, b_tau};
 st_bayesrc_selection_validate(annotation, rows, &outcomes, state, hyper);
 for (arma::uword j = 0; j < state.delta.n_elem; ++j) {
  if (state.delta(j) == 0u) state.alpha.row(j + 1u).zeros();
 }
 arma::ivec fixed;
 const arma::ivec* fixed_pointer = nullptr;
 if (fixed_delta.size() > 0) {
  const arma::uvec fixed_unsigned = st_bayesrc_selection_delta_from_r(
   fixed_delta, annotation.n_cols);
  fixed = arma::conv_to<arma::ivec>::from(fixed_unsigned);
  state.delta = fixed_unsigned;
  fixed_pointer = &fixed;
 }

 const int retained = iterations - burn;
 arma::umat delta_draws(retained, annotation.n_cols, arma::fill::zeros);
 arma::cube alpha_draws(retained, annotation.n_cols + 1u,
                        rows.size(), arma::fill::zeros);
 arma::vec pi_a_draws(retained, arma::fill::zeros);
 arma::mat tau2_draws(retained, rows.size(), arma::fill::zeros);
 arma::uvec included_draws(retained, arma::fill::zeros);
 arma::mat q_sum(annotation.n_rows, rows.size(), arma::fill::zeros);
 arma::mat pi_sum(annotation.n_rows, rows.size() + 1u, arma::fill::zeros);
 arma::umat switches(2u, annotation.n_cols, arma::fill::zeros);
 arma::uvec previous = state.delta;
 std::mt19937 generator(static_cast<std::mt19937::result_type>(seed));
 for (int iteration = 0; iteration < iterations; ++iteration) {
  const auto latent = st_bayesrc_selection_sample_latent(
   annotation, rows, outcomes, state.alpha, generator);
  st_bayesrc_selection_delta_sweep(
   annotation, rows, latent, state, generator, fixed_pointer);
  st_bayesrc_selection_blocked_redraw(
   annotation, rows, latent, state, generator);
  st_bayesrc_selection_update_hyperparameters(state, hyper, generator);
  if (iteration >= burn) {
   const arma::uword draw = static_cast<arma::uword>(iteration - burn);
   delta_draws.row(draw) = state.delta.t();
   for (arma::uword stick = 0; stick < state.alpha.n_cols; ++stick) {
    alpha_draws.slice(stick).row(draw) = state.alpha.col(stick).t();
   }
   pi_a_draws(draw) = state.pi_a;
   tau2_draws.row(draw) = state.tau2.t();
   included_draws(draw) = arma::accu(state.delta);
   q_sum += st_bayesrc_selection_compute_q(annotation, state.alpha);
   pi_sum += st_bayesrc_selection_compute_pi(
    annotation, state.alpha, probability_floor);
   if (draw > 0u) {
    for (arma::uword j = 0; j < state.delta.n_elem; ++j) {
     if (previous(j) == 0u && state.delta(j) == 1u) switches(0u, j)++;
     if (previous(j) == 1u && state.delta(j) == 0u) switches(1u, j)++;
    }
   }
  }
  previous = state.delta;
 }
 return Rcpp::List::create(
  Rcpp::Named("delta_draws") = delta_draws,
  Rcpp::Named("alpha_draws") = alpha_draws,
  Rcpp::Named("pi_A_draws") = pi_a_draws,
  Rcpp::Named("tau2_draws") = tau2_draws,
  Rcpp::Named("included_draws") = included_draws,
  Rcpp::Named("q_mean") = q_sum / static_cast<double>(retained),
  Rcpp::Named("component_probability_mean") =
   pi_sum / static_cast<double>(retained),
  Rcpp::Named("switches") = switches,
  Rcpp::Named("delta_final") = state.delta,
  Rcpp::Named("alpha_final") = state.alpha,
  Rcpp::Named("pi_A_final") = state.pi_a,
  Rcpp::Named("tau2_final") = state.tau2
 );
}

// Deterministic development hook for validating the exact two-marker
// component conditional. It does not mutate a production chain.
// [[Rcpp::export(name = ".st_bayesrc_pairwise_conditional")]]
Rcpp::List st_bayesrc_pairwise_conditional_test(
 arma::rowvec prior_i, arma::rowvec prior_j, arma::vec gamma,
 double marker_variance, double prior_scale_i, double prior_scale_j,
 double residual_variance, double diagonal_i, double diagonal_j,
 double cross_product, double score_i, double score_j
) {
 const auto states = st_bayesrc_pairwise_conditional(
  prior_i, prior_j, gamma, marker_variance, prior_scale_i, prior_scale_j,
  residual_variance, diagonal_i, diagonal_j, cross_product, score_i, score_j);
 const arma::uword count = static_cast<arma::uword>(states.size());
 arma::imat component(count, 2, arma::fill::zeros);
 arma::vec log_weight(count), probability(count);
 arma::mat mean(count, 2, arma::fill::zeros);
 arma::cube covariance(2, 2, count, arma::fill::zeros);
 for (arma::uword index = 0; index < count; ++index) {
  const auto& state = states[static_cast<std::size_t>(index)];
  component(index, 0) = state.component_i;
  component(index, 1) = state.component_j;
  log_weight(index) = state.log_weight;
  probability(index) = state.probability;
  mean(index, 0) = state.mean_i;
  mean(index, 1) = state.mean_j;
  covariance(0, 0, index) = state.covariance_ii;
  covariance(0, 1, index) = state.covariance_ij;
  covariance(1, 0, index) = state.covariance_ij;
  covariance(1, 1, index) = state.covariance_jj;
 }
 return Rcpp::List::create(
  Rcpp::Named("component") = component,
  Rcpp::Named("log_weight") = log_weight,
  Rcpp::Named("probability") = probability,
  Rcpp::Named("mean") = mean,
  Rcpp::Named("covariance") = covariance);
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
  arma::mat alpha_init, arma::vec sigmaSqAlpha_init,
  arma::mat intercept_prior_resolved,
  double sigmaSqAlpha_a, double sigmaSqAlpha_b, double pi_floor, double nub,
  double nue, bool updateAlpha, bool updateB, bool updateE,
  int alpha_update_every, double adjE, std::vector<int> n, int nit,
  int nburn, int nthin, int ncores, int seed, int nchains = 1,
  bool keep_chains = false,
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds = R_NilValue,
  bool updateLDswap = false, double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8, int ld_swap_max_friends = 50,
  int ld_swap_moves = 1,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale = R_NilValue,
  bool estimate_maf_effect_s = false, double maf_effect_s_init = 0.0,
  Rcpp::NumericVector maf_effect_s_prior = Rcpp::NumericVector::create(-3.0, 2.0),
  double maf_effect_s_proposal_sd = 0.35,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h = R_NilValue,
  Rcpp::IntegerVector convergence_markers=Rcpp::IntegerVector::create(),
  bool convergence_annotations=false, bool convergence_b=false,
  bool convergence_d=false, bool convergence_component=false
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
  alpha_init, sigmaSqAlpha_init, intercept_prior_resolved, sigmaSqAlpha_a,
  sigmaSqAlpha_b, pi_floor, nub, nue, updateAlpha, updateB, updateE,
  alpha_update_every, adjE, n, nit, nburn, nthin, ncores, seed, nchains,
  keep_chains, chain_seeds, updateLDswap, ld_swap_prob, ld_swap_r2,
  ld_swap_max_friends, ld_swap_moves, maf_effect_s_prior_scale,
  estimate_maf_effect_s, maf_effect_s_init, maf_effect_s_prior,
  maf_effect_s_proposal_sd, maf_effect_s_log_h, convergence_markers,
  convergence_annotations, convergence_b, convergence_d,
  convergence_component, 0, std::vector<double>(),
  "global_projected_legacy", "fixVe", 1.1, 0.7, false,
  StBayesRCSelectionGenomicConfig{},
  make_csr_operator
 );
}

// Internal Phase-4B qualification binding.  This composes the validated CSR
// genomic engine with the SBayesRC-S annotation policy without exposing a
// supported public model identifier or wrapper.
// [[Rcpp::export(name = ".st_sbayesrc_selection_csr")]]
Rcpp::List st_sbayesrc_selection_csr_internal(
  std::vector<std::vector<double>> wy, std::vector<std::vector<double>> ww,
  std::vector<double> yy, std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> comp_init, bool use_comp_init,
  std::vector<std::vector<double>> r_init, bool use_r_init,
  std::string ld_prefix, arma::mat B, arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior, arma::mat A, arma::vec gamma,
  arma::mat alpha_init, arma::uvec delta_init, double pi_A_init,
  arma::vec tau2_init, double a_pi, double b_pi, double a_tau, double b_tau,
  Rcpp::IntegerVector fixed_delta, bool update_hierarchy,
  arma::mat intercept_prior_resolved,
  double pi_floor, double nub, double nue, bool updateB, bool updateE,
  double adjE, std::vector<int> n, int nit, int nburn, int nthin,
  int ncores, int seed, int nchains = 1,
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds = R_NilValue
) {
 if (A.n_cols < 2u || alpha_init.n_rows != A.n_cols ||
     alpha_init.n_cols + 1u != gamma.n_elem ||
     delta_init.n_elem + 1u != A.n_cols ||
     tau2_init.n_elem != alpha_init.n_cols) {
  throw std::invalid_argument("invalid internal genomic SBayesRC-S dimensions");
 }
 StBayesRCSelectionGenomicConfig selection_config;
 selection_config.enabled = true;
 selection_config.delta_init = delta_init;
 selection_config.pi_a_init = pi_A_init;
 selection_config.tau2_init = tau2_init;
 selection_config.hyper = StBayesRCSelectionHyperParameters{
  a_pi, b_pi, a_tau, b_tau};
 if (fixed_delta.size() > 0) {
  selection_config.fixed_delta = true;
  selection_config.fixed_delta_value = Rcpp::as<arma::ivec>(fixed_delta);
 }
 const arma::vec legacy_sigma = tau2_init;
 auto make_csr_operator = [&](int m, const std::vector<double>& xx,
                              const arma::rowvec& xx_row, arma::mat& wy_mat,
                              bool update_ld_swap, double r2, int max_friends) {
  (void)wy_mat; (void)r2; (void)max_friends;
  if (update_ld_swap) throw std::runtime_error(
   "LD-swap is unavailable in the internal SBayesRC-S qualification binding");
  STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
  CsrOperator op(ld, xx_row);
  SBayesRCLDLDFriends friends;
  friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  return SBayesRCOperatorContext<CsrOperator>(
   std::move(op), std::move(friends), Rcpp::List::create());
 };
 return stblr_cpg_omp_csr_sbayesrc_impl(
  wy, ww, yy, b_init, comp_init, use_comp_init, r_init, use_r_init,
  false, ld_prefix, B, E, ssb_prior, sse_prior, A, gamma, alpha_init,
  legacy_sigma, intercept_prior_resolved, 2.0, 2.0, pi_floor, nub, nue,
  update_hierarchy, updateB, updateE, 1, adjE, n, nit, nburn, nthin, ncores, seed,
  nchains, true, chain_seeds, false, 0.0, 0.8, 50, 0, R_NilValue,
  false, 0.0, Rcpp::NumericVector::create(-3.0, 2.0), 0.35, R_NilValue,
  Rcpp::IntegerVector::create(), true, false, false, true, 0,
  std::vector<double>(), "global_projected_legacy", "fixVe", 1.1, 0.7,
  false, selection_config, make_csr_operator);
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_sbayesrc_block_eigen(
  std::vector<std::vector<double>> wy, std::vector<std::vector<double>> ww,
  std::vector<double> yy, std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> comp_init, bool use_comp_init,
  std::vector<std::vector<double>> r_init, bool use_r_init,
  bool rebuild_r_before_updateE, arma::mat B,
  arma::mat E, std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior, arma::mat A, arma::vec gamma,
  arma::mat alpha_init, arma::vec sigmaSqAlpha_init,
  arma::mat intercept_prior_resolved,
  double sigmaSqAlpha_a, double sigmaSqAlpha_b, double pi_floor, double nub,
  double nue, bool updateAlpha, bool updateB, bool updateE,
  int alpha_update_every, double adjE, std::vector<int> n, int nit,
  int nburn, int nthin, int ncores, int seed, int nchains = 1,
  bool keep_chains = false,
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds = R_NilValue,
  bool updateLDswap = false, double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8, int ld_swap_max_friends = 50,
  int ld_swap_moves = 1,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale = R_NilValue,
  bool estimate_maf_effect_s = false, double maf_effect_s_init = 0.0,
  Rcpp::NumericVector maf_effect_s_prior = Rcpp::NumericVector::create(-3.0, 2.0),
  double maf_effect_s_proposal_sd = 0.35,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h = R_NilValue,
  Rcpp::IntegerVector convergence_markers=Rcpp::IntegerVector::create(),
  bool convergence_annotations=false, bool convergence_b=false,
  bool convergence_d=false, bool convergence_component=false,
  Rcpp::CharacterVector bed_files = Rcpp::CharacterVector::create(),
  int n_bed = 0, Rcpp::List cls = R_NilValue,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::NumericVector af = Rcpp::NumericVector::create(),
  Rcpp::IntegerVector block_start = Rcpp::IntegerVector::create(),
  std::string eigen_filter = "hard_truncate", double eigen_tau = 0.01,
  double eigen_eta = 0.0, std::string representation = "dense_reconstructed",
  double eigen_prop = 0.995, Rcpp::List low_rank_residual_config = R_NilValue
) {
 const int low_rank_residual_rebuild_every =
  low_rank_residual_config.size() > 0 ? Rcpp::as<int>(
   low_rank_residual_config["rebuild_every"]) : 100;
 Rcpp::List block_residual_config = low_rank_residual_config.size() > 0 ?
  Rcpp::as<Rcpp::List>(low_rank_residual_config["block_residual"]) :
  Rcpp::List::create();
 const bool has_block_residual_config = block_residual_config.size() > 0;
 Rcpp::NumericVector phenotype_variance = has_block_residual_config ?
  Rcpp::as<Rcpp::NumericVector>(block_residual_config["phenotype_variance"]) :
  Rcpp::NumericVector::create();
 const std::string residual_policy = has_block_residual_config ?
  Rcpp::as<std::string>(block_residual_config["residual_policy"]) :
  "global_projected_legacy";
 const std::string block_ve_mode = has_block_residual_config ?
  Rcpp::as<std::string>(block_residual_config["block_ve_mode"]) : "fixVe";
 const double resam_thresh = has_block_residual_config ?
  Rcpp::as<double>(block_residual_config["resam_thresh"]) : 1.1;
 const double minimum_ve_ratio = has_block_residual_config ?
  Rcpp::as<double>(block_residual_config["minimum_ve_ratio"]) : 0.7;
 const bool block_ve_keep_history = has_block_residual_config ?
  Rcpp::as<bool>(block_residual_config["block_ve_keep_history"]) : false;
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
 if (representation != "low_rank" && representation != "dense_reconstructed")
  throw std::runtime_error("unknown block-eigen representation.");
 if (representation == "low_rank" && use_r_init)
  throw std::runtime_error("marker-space r_init is unavailable for the low-rank representation.");

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
  BlockEigenDispatchOperator op;
  op.low_rank = representation == "low_rank";
  Rcpp::List diag;
  if (op.low_rank) {
   std::vector<BlockLowRankDiag> block_diag;
   op.retained = build_block_low_rank(
    G, af_cpp, starts, eigen_prop, wy_mat, ncores, &block_diag
   );
   diag = Rcpp::List::create(
    Rcpp::Named("blocks") = block_low_rank_diagnostics_to_data_frame(block_diag),
    Rcpp::Named("operator_contract") = "block_low_rank_v1",
    Rcpp::Named("operator_representation") = "low_rank",
    Rcpp::Named("operator_scale_contract") = "general_cross_product",
    Rcpp::Named("eigen_policy") = "cumulative_positive_mass",
    Rcpp::Named("eigen_prop") = eigen_prop,
    Rcpp::Named("build") = block_low_rank_build_metadata(op.retained)
   );
  } else {
   std::vector<BlockEigenDiag> block_diag;
   op.dense = build_block_eigen(
    G, af_cpp, starts, mode, eigen_tau, eigen_eta, wy_mat, ncores, &block_diag
   );
   diag = Rcpp::List::create(
    Rcpp::Named("blocks") = block_eigen_diagnostics_to_data_frame(block_diag),
    Rcpp::Named("operator_contract") = "block_dense_reconstructed_v1",
    Rcpp::Named("operator_representation") = "dense_reconstructed"
   );
  }
  SBayesRCLDLDFriends friends;
  friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  return SBayesRCOperatorContext<BlockEigenDispatchOperator>(
   std::move(op), std::move(friends), diag
  );
 };
 return stblr_cpg_omp_csr_sbayesrc_impl(
  wy, ww, yy, b_init, comp_init, use_comp_init, r_init, use_r_init,
  rebuild_r_before_updateE, "", B, E, ssb_prior, sse_prior, A, gamma,
  alpha_init, sigmaSqAlpha_init, intercept_prior_resolved, sigmaSqAlpha_a,
  sigmaSqAlpha_b, pi_floor, nub, nue, updateAlpha, updateB, updateE,
  alpha_update_every, adjE, n, nit, nburn, nthin, ncores, seed, nchains,
  keep_chains, chain_seeds, updateLDswap, ld_swap_prob, ld_swap_r2,
  ld_swap_max_friends, ld_swap_moves, maf_effect_s_prior_scale,
  estimate_maf_effect_s, maf_effect_s_init, maf_effect_s_prior,
  maf_effect_s_proposal_sd, maf_effect_s_log_h, convergence_markers,
  convergence_annotations, convergence_b, convergence_d,
  convergence_component, low_rank_residual_rebuild_every,
  Rcpp::as<std::vector<double>>(phenotype_variance),
  residual_policy, block_ve_mode, resam_thresh, minimum_ve_ratio,
  block_ve_keep_history, StBayesRCSelectionGenomicConfig{}, make_block_operator
 );
}

// [[Rcpp::export(.st_sbayesrc_invalid_scale_diagnostic_fixture)]]
void st_sbayesrc_invalid_scale_diagnostic_fixture(std::string path) {
 sblr::core::CsrSBayesRCFailureState state;
 state.captured = true;
 state.operator_name = "csr";
 state.trait = 0; state.chain = 0; state.iteration = 1;
 state.internal_iteration = 0; state.sample_size = 8;
 state.residual_df = 4.0; state.yy = 1.0; state.sse_prior = 0.1;
 state.prior_contribution = 0.4;
 state.maintained_sse = -1.4; state.rebuilt_sse = -1.4;
 state.quadratic_sse = -1.4; state.maintained_scale = -1.0;
 state.rebuilt_scale = -1.0; state.quadratic_scale = -1.0;
 state.effects = arma::rowvec(1, arma::fill::ones);
 state.residual = arma::rowvec(1, arma::fill::zeros);
 state.rebuilt_residual = state.residual;
 state.score = arma::rowvec(1, arma::fill::ones);
 state.diagonal = arma::rowvec(1, arma::fill::ones);
 state.component = arma::Row<int>(1, arma::fill::ones);
 state.effects_finite = true; state.residual_finite = true;
 state.rebuilt_residual_finite = true;
 sblr::core::write_csr_sbayesrc_failure_state(state, path);
 double ve = 1.0;
 std::mt19937 generator(718);
 try {
  sampleE_ST_operator_from_scale(4.0, ve, -1.0, 8, generator);
 } catch (const std::exception&) {
  throw std::runtime_error(
   "sampleE_ST_operator: invalid projected residual scale at iteration 1 "
   "(internal 0); maintained_scale=-1.000000");
 }
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
// #endif
//
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
// #endif
//
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
