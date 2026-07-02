// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "st_chain_utils.h"
#include "st_csr_common.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

inline double computeLE_bayesr_ST_csr(
  int m,
  const arma::rowvec& b,
  const arma::rowvec& ww,
  int n
) {
 double vle = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double bi = b(iu);
  if (bi != 0.0) vle += ww(iu) * bi * bi;
 }

 return vle / static_cast<double>(n);
}

inline void residual_sse_terms_bayesr_ST_csr(
  int m,
  const arma::rowvec& b,
  const arma::rowvec& wy,
  const arma::rowvec& r,
  double yy,
  double& bwy,
  double& br,
  double& bXb,
  double& sse
) {
 bwy = 0.0;
 br = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  bwy += b(iu) * wy(iu);
  br += b(iu) * r(iu);
 }

 bXb = bwy - br;
 sse = yy - bwy - br;
}

struct BayesRUpdateEDiagnostics {
 double bwy = 0.0;
 double br = 0.0;
 double bXb = 0.0;
 double sse = 0.0;
 double residual_scale = 0.0;
 int nonzero_components = 0;
 double max_abs_b = 0.0;
 double max_abs_r = 0.0;
};

inline BayesRUpdateEDiagnostics residual_diagnostics_bayesr_ST_csr(
  int m,
  double nue,
  const arma::rowvec& b,
  const arma::rowvec& r,
  const arma::Row<int>& comp,
  const arma::rowvec& wy,
  double sse_prior,
  double yy
) {
 BayesRUpdateEDiagnostics out;
 residual_sse_terms_bayesr_ST_csr(m, b, wy, r, yy, out.bwy, out.br, out.bXb, out.sse);
 out.residual_scale = out.sse + nue * sse_prior;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (comp(iu) > 0) ++out.nonzero_components;
  out.max_abs_b = std::max(out.max_abs_b, std::abs(b(iu)));
  out.max_abs_r = std::max(out.max_abs_r, std::abs(r(iu)));
 }

 return out;
}

inline void ensure_null_effects_bayesr_ST_csr(
  int m,
  int trait,
  int chain,
  int iter,
  const arma::rowvec& b,
  const arma::Row<int>& comp
) {
 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (comp(iu) == 0 && b(iu) != 0.0) {
   throw std::runtime_error(
    "BayesR CSR state mismatch: null component has nonzero effect. trait=" +
    std::to_string(trait) +
    ", chain=" + std::to_string(chain) +
    ", iter=" + std::to_string(iter) +
    ", marker=" + std::to_string(i) +
    ", b=" + std::to_string(b(iu))
   );
  }
 }
}

inline void check_residual_scale_bayesr_ST_csr(
  int m,
  int trait,
  int chain,
  int iter,
  double nue,
  double ve,
  double vb,
  const arma::vec& mixture_var,
  const arma::rowvec& b,
  const arma::rowvec& r,
  const arma::Row<int>& comp,
  const arma::rowvec& wy,
  double sse_prior,
  double yy,
  int n,
  double adjE
) {
 const BayesRUpdateEDiagnostics diag = residual_diagnostics_bayesr_ST_csr(
  m,
  nue,
  b,
  r,
  comp,
  wy,
  sse_prior,
  yy
 );

 if (std::isfinite(diag.residual_scale) && diag.residual_scale > 0.0) return;

 const double b_min = b.n_elem > 0 ? b.min() : std::numeric_limits<double>::quiet_NaN();
 const double b_max = b.n_elem > 0 ? b.max() : std::numeric_limits<double>::quiet_NaN();
 const double b_sum = b.n_elem > 0 ? arma::accu(b) : std::numeric_limits<double>::quiet_NaN();
 const double r_min = r.n_elem > 0 ? r.min() : std::numeric_limits<double>::quiet_NaN();
 const double r_max = r.n_elem > 0 ? r.max() : std::numeric_limits<double>::quiet_NaN();
 const double mix_min = mixture_var.n_elem > 0 ? mixture_var.min() : std::numeric_limits<double>::quiet_NaN();
 const double mix_max = mixture_var.n_elem > 0 ? mixture_var.max() : std::numeric_limits<double>::quiet_NaN();

 throw std::runtime_error(
  "sampleE_ST_csr: invalid residual scale in BayesR CSR. trait=" +
  std::to_string(trait) +
  ", chain=" + std::to_string(chain) +
  ", iter=" + std::to_string(iter) +
  ", yy=" + std::to_string(yy) +
  ", n=" + std::to_string(n) +
  ", adjE=" + std::to_string(adjE) +
  ", ve=" + std::to_string(ve) +
  ", vb=" + std::to_string(vb) +
  ", mixture_var_min=" + std::to_string(mix_min) +
  ", mixture_var_max=" + std::to_string(mix_max) +
  ", sse_prior=" + std::to_string(sse_prior) +
  ", bwy=" + std::to_string(diag.bwy) +
  ", br=" + std::to_string(diag.br) +
  ", bXb=" + std::to_string(diag.bXb) +
  ", sse=" + std::to_string(diag.sse) +
  ", residual_scale=" + std::to_string(diag.residual_scale) +
  ", r_finite=" + std::to_string(r.is_finite() ? 1 : 0) +
  ", r_min=" + std::to_string(r_min) +
  ", r_max=" + std::to_string(r_max) +
  ", nonzero_components=" + std::to_string(diag.nonzero_components) +
  ", b_finite=" + std::to_string(b.is_finite() ? 1 : 0) +
  ", b_min=" + std::to_string(b_min) +
  ", b_max=" + std::to_string(b_max) +
  ", b_sum=" + std::to_string(b_sum)
 );
}

inline double logsumexp_bayesr(const std::vector<double>& x) {
 double mx = -std::numeric_limits<double>::infinity();
 for (double v : x) mx = std::max(mx, v);
 if (!std::isfinite(mx)) return mx;

 double s = 0.0;
 for (double v : x) s += std::exp(v - mx);
 return mx + std::log(s);
}

inline int sample_categorical_logprob_bayesr(
  const std::vector<double>& logp,
  std::mt19937& gen
) {
 const double lse = logsumexp_bayesr(logp);
 if (!std::isfinite(lse)) {
  throw std::runtime_error("sample_categorical_logprob_bayesr: invalid log probabilities.");
 }

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 const double u = runif(gen);
 double cum = 0.0;

 for (std::size_t k = 0; k < logp.size(); ++k) {
  cum += std::exp(logp[k] - lse);
  if (u <= cum) return static_cast<int>(k);
 }

 return static_cast<int>(logp.size() - 1L);
}

inline void sampleBetaR_ST_csr(
  int i,
  const std::vector<double>& pi,
  const arma::vec& mixture_var,
  double vb,
  const arma::rowvec& prior_scale,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const STLDCSR& ld,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double wi = ww(iu);

 if (!std::isfinite(wi) || wi <= 0.0) {
  throw std::runtime_error("sampleBetaR_ST_csr: invalid ww value.");
 }
 if (!std::isfinite(vb) || vb <= 0.0) {
  throw std::runtime_error("sampleBetaR_ST_csr: invalid vb.");
 }

 const int K = static_cast<int>(mixture_var.n_elem);
 if (static_cast<int>(pi.size()) != K) {
  throw std::runtime_error("sampleBetaR_ST_csr: pi and mixture_var length mismatch.");
 }

 const double vei_safe = std::max(vei_i, 1e-300);
 const double score = r(iu) + wi * b(iu);
 const double scale_i = prior_scale(iu);
 if (!std::isfinite(scale_i) || scale_i <= 0.0) {
  throw std::runtime_error("sampleBetaR_ST_csr: invalid prior_scale value.");
 }
 std::vector<double> logp(static_cast<std::size_t>(K));

 for (int k = 0; k < K; ++k) {
  const double pik = std::max(pi[static_cast<std::size_t>(k)], 1e-300);
  const double ck = mixture_var(static_cast<arma::uword>(k));

  if (!std::isfinite(ck) || ck < 0.0) {
   throw std::runtime_error("sampleBetaR_ST_csr: mixture_var must be non-negative.");
  }

  if (ck <= 0.0) {
   logp[static_cast<std::size_t>(k)] = std::log(pik);
  } else {
   const double vbk = std::max(vb * ck * scale_i, 1e-300);
   const double denom = std::max(vei_safe + wi * vbk, 1e-300);
   const double logBF =
    0.5 * std::log(vei_safe / denom) +
    0.5 * score * score * vbk / (vei_safe * denom);
   logp[static_cast<std::size_t>(k)] = std::log(pik) + logBF;
  }
 }

 const int k_new = sample_categorical_logprob_bayesr(logp, gen);

 double b_new = 0.0;
 const double ck_new = mixture_var(static_cast<arma::uword>(k_new));

 if (ck_new > 0.0) {
  std::normal_distribution<double> norm01(0.0, 1.0);
  const double vbk = std::max(vb * ck_new * scale_i, 1e-300);
  const double lhs = wi + vei_safe / vbk;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  r(iu) -= wi * diff;

  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   r(static_cast<arma::uword>(j)) -=
    static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
  }
 }

 b(iu) = b_new;
 comp(iu) = k_new;
}

inline void sampleBetaR_ST_csr_unscaled(
  int i,
  const std::vector<double>& pi,
  const arma::vec& mixture_var,
  double vb,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const STLDCSR& ld,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double wi = ww(iu);

 if (!std::isfinite(wi) || wi <= 0.0) {
  throw std::runtime_error("sampleBetaR_ST_csr_unscaled: invalid ww value.");
 }
 if (!std::isfinite(vb) || vb <= 0.0) {
  throw std::runtime_error("sampleBetaR_ST_csr_unscaled: invalid vb.");
 }

 const int K = static_cast<int>(mixture_var.n_elem);
 if (static_cast<int>(pi.size()) != K) {
  throw std::runtime_error("sampleBetaR_ST_csr_unscaled: pi and mixture_var length mismatch.");
 }

 const double vei_safe = std::max(vei_i, 1e-300);
 const double score = r(iu) + wi * b(iu);
 std::vector<double> logp(static_cast<std::size_t>(K));

 for (int k = 0; k < K; ++k) {
  const double pik = std::max(pi[static_cast<std::size_t>(k)], 1e-300);
  const double ck = mixture_var(static_cast<arma::uword>(k));

  if (!std::isfinite(ck) || ck < 0.0) {
   throw std::runtime_error("sampleBetaR_ST_csr_unscaled: mixture_var must be non-negative.");
  }

  if (ck <= 0.0) {
   logp[static_cast<std::size_t>(k)] = std::log(pik);
  } else {
   const double vbk = std::max(vb * ck, 1e-300);
   const double denom = std::max(vei_safe + wi * vbk, 1e-300);
   const double logBF =
    0.5 * std::log(vei_safe / denom) +
    0.5 * score * score * vbk / (vei_safe * denom);
   logp[static_cast<std::size_t>(k)] = std::log(pik) + logBF;
  }
 }

 const int k_new = sample_categorical_logprob_bayesr(logp, gen);

 double b_new = 0.0;
 const double ck_new = mixture_var(static_cast<arma::uword>(k_new));

 if (ck_new > 0.0) {
  std::normal_distribution<double> norm01(0.0, 1.0);
  const double vbk = std::max(vb * ck_new, 1e-300);
  const double lhs = wi + vei_safe / vbk;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  r(iu) -= wi * diff;

  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   r(static_cast<arma::uword>(j)) -=
    static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
  }
 }

 b(iu) = b_new;
 comp(iu) = k_new;
}

struct BayesRLDLDFriends {
 std::vector<uint64_t> ptr;
 std::vector<int> idx;
 std::vector<double> r2;
};

inline BayesRLDLDFriends build_ld_swap_friends_bayesr_ST_csr(
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

 BayesRLDLDFriends friends;
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

inline void set_marker_state_bayesr_ST_csr(
  int i,
  double b_new,
  int comp_new,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const STLDCSR& ld
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  r(iu) -= ww(iu) * diff;

  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   r(static_cast<arma::uword>(j)) -=
    static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
  }
 }

 b(iu) = b_new;
 comp(iu) = comp_new;
}

inline int count_null_ld_friends_bayesr_ST_csr(
  int i,
  const arma::Row<int>& comp,
  const arma::rowvec& b,
  const BayesRLDLDFriends& friends
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

inline int collect_ld_swap_candidates_bayesr_ST_csr(
  int m,
  const arma::Row<int>& comp,
  const arma::rowvec& b,
  const BayesRLDLDFriends& friends,
  std::vector<int>& candidates,
  std::vector<int>& n_null
) {
 candidates.clear();
 n_null.clear();

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (comp(iu) <= 0 || b(iu) == 0.0) continue;

  const int nf = count_null_ld_friends_bayesr_ST_csr(i, comp, b, friends);
  if (nf > 0) {
   candidates.push_back(i);
   n_null.push_back(nf);
  }
 }

 return static_cast<int>(candidates.size());
}

inline void throw_ld_swap_error_bayesr_ST_csr(
  int trait,
  int chain,
  int iter,
  int source,
  int target,
  int component_moved,
  double sse_old,
  double sse_new,
  double vei,
  double log_q_forward,
  double log_q_reverse,
  const arma::rowvec& b,
  const arma::rowvec& r
) {
 throw std::runtime_error(
  "BayesR CSR LD-swap invalid proposal. trait=" +
  std::to_string(trait) +
  ", chain=" + std::to_string(chain) +
  ", iter=" + std::to_string(iter) +
  ", source=" + std::to_string(source) +
  ", target=" + std::to_string(target) +
  ", component_moved=" + std::to_string(component_moved) +
  ", sse_old=" + std::to_string(sse_old) +
  ", sse_new=" + std::to_string(sse_new) +
  ", vei=" + std::to_string(vei) +
  ", log_q_forward=" + std::to_string(log_q_forward) +
  ", log_q_reverse=" + std::to_string(log_q_reverse) +
  ", b_finite=" + std::to_string(b.is_finite() ? 1 : 0) +
  ", r_finite=" + std::to_string(r.is_finite() ? 1 : 0)
 );
}

inline bool attempt_ld_swap_bayesr_ST_csr(
  int m,
  int trait,
  int chain,
  int iter,
  double vei,
  double vb,
  double yy,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  const arma::vec& mixture_var,
  const arma::rowvec& prior_scale,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const STLDCSR& ld,
  const BayesRLDLDFriends& friends,
  std::mt19937& gen,
  bool& attempted
) {
 attempted = false;
 if (!std::isfinite(vei) || vei <= 0.0) return false;

 std::vector<int> candidates;
 std::vector<int> n_null;
 const int n_candidates =
  collect_ld_swap_candidates_bayesr_ST_csr(m, comp, b, friends, candidates, n_null);
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
 if (comp_j_old >= static_cast<int>(mixture_var.n_elem)) {
  throw std::runtime_error("attempt_ld_swap_bayesr_ST_csr: component index out of range.");
 }
 const double ck = mixture_var(static_cast<arma::uword>(comp_j_old));
 if (!std::isfinite(ck) || ck <= 0.0) {
  throw std::runtime_error("attempt_ld_swap_bayesr_ST_csr: active component has invalid mixture_var.");
 }
 const double vb_j = std::max(vb * ck * prior_scale(ju), 1e-300);
 const double vb_k = std::max(vb * ck * prior_scale(ku), 1e-300);

 attempted = true;
 const double sse_old = residual_diagnostics_bayesr_ST_csr(
  m, 0.0, b, r, comp, wy, 0.0, yy
 ).sse;
 if (!std::isfinite(sse_old)) {
  throw_ld_swap_error_bayesr_ST_csr(
   trait, chain, iter, j, k, comp_j_old, sse_old,
   std::numeric_limits<double>::quiet_NaN(), vei,
   std::numeric_limits<double>::quiet_NaN(),
   std::numeric_limits<double>::quiet_NaN(), b, r
  );
 }

 const arma::rowvec r_old = r;
 set_marker_state_bayesr_ST_csr(j, 0.0, 0, ww, r, b, comp, ld);
 set_marker_state_bayesr_ST_csr(k, b_j_old, comp_j_old, ww, r, b, comp, ld);

 const double sse_new = residual_diagnostics_bayesr_ST_csr(
  m, 0.0, b, r, comp, wy, 0.0, yy
 ).sse;

 std::vector<int> reverse_candidates;
 std::vector<int> reverse_n_null;
 const int n_reverse_candidates =
  collect_ld_swap_candidates_bayesr_ST_csr(
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

 if (!std::isfinite(sse_new) ||
     !std::isfinite(log_q_forward) ||
     !std::isfinite(log_q_reverse)) {
  r = r_old;
  b(ju) = b_j_old;
  b(ku) = b_k_old;
  comp(ju) = comp_j_old;
  comp(ku) = comp_k_old;
  throw_ld_swap_error_bayesr_ST_csr(
   trait, chain, iter, j, k, comp_j_old, sse_old, sse_new,
   vei, log_q_forward, log_q_reverse, b, r
  );
 }

 // With constant BayesR prior variance, active/null LD-swap relocations
 // cancel the prior density. Fixed selection_s makes the moved component's
 // N(0, vb * mixture_var_m * prior_scale_j) density marker-specific unless
 // prior_scale_j == prior_scale_k, so it must enter the MH ratio.
 const double log_prior_ratio =
  -0.5 * (std::log(vb_k / vb_j) +
          b_j_old * b_j_old * (1.0 / vb_k - 1.0 / vb_j));
 const double log_alpha =
  -0.5 * (sse_new - sse_old) / vei +
  log_prior_ratio + log_q_reverse - log_q_forward;

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

inline bool attempt_ld_swap_bayesr_ST_csr_unscaled(
  int m,
  int trait,
  int chain,
  int iter,
  double vei,
  double yy,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& comp,
  const STLDCSR& ld,
  const BayesRLDLDFriends& friends,
  std::mt19937& gen,
  bool& attempted
) {
 attempted = false;
 if (!std::isfinite(vei) || vei <= 0.0) return false;

 std::vector<int> candidates;
 std::vector<int> n_null;
 const int n_candidates =
  collect_ld_swap_candidates_bayesr_ST_csr(m, comp, b, friends, candidates, n_null);
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
 const double sse_old = residual_diagnostics_bayesr_ST_csr(
  m, 0.0, b, r, comp, wy, 0.0, yy
 ).sse;
 if (!std::isfinite(sse_old)) {
  throw_ld_swap_error_bayesr_ST_csr(
   trait, chain, iter, j, k, comp_j_old, sse_old,
   std::numeric_limits<double>::quiet_NaN(), vei,
   std::numeric_limits<double>::quiet_NaN(),
   std::numeric_limits<double>::quiet_NaN(), b, r
  );
 }

 const arma::rowvec r_old = r;
 set_marker_state_bayesr_ST_csr(j, 0.0, 0, ww, r, b, comp, ld);
 set_marker_state_bayesr_ST_csr(k, b_j_old, comp_j_old, ww, r, b, comp, ld);

 const double sse_new = residual_diagnostics_bayesr_ST_csr(
  m, 0.0, b, r, comp, wy, 0.0, yy
 ).sse;

 std::vector<int> reverse_candidates;
 std::vector<int> reverse_n_null;
 const int n_reverse_candidates =
  collect_ld_swap_candidates_bayesr_ST_csr(
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

 if (!std::isfinite(sse_new) ||
     !std::isfinite(log_q_forward) ||
     !std::isfinite(log_q_reverse)) {
  r = r_old;
  b(ju) = b_j_old;
  b(ku) = b_k_old;
  comp(ju) = comp_j_old;
  comp(ku) = comp_k_old;
  throw_ld_swap_error_bayesr_ST_csr(
   trait, chain, iter, j, k, comp_j_old, sse_old, sse_new,
   vei, log_q_forward, log_q_reverse, b, r
  );
 }

 const double log_alpha =
  -0.5 * (sse_new - sse_old) / vei +
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

inline void sampleB_bayesr_ST_csr(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  const arma::vec& mixture_var,
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
   if (k >= static_cast<int>(mixture_var.n_elem)) {
    throw std::runtime_error("sampleB_bayesr_ST_csr: component index out of range.");
   }
   const double ck = mixture_var(static_cast<arma::uword>(k));
   if (!std::isfinite(ck) || ck <= 0.0) {
    throw std::runtime_error("sampleB_bayesr_ST_csr: active component has invalid mixture_var.");
   }
   const double scale_i = prior_scale(iu);
   if (!std::isfinite(scale_i) || scale_i <= 0.0) {
    throw std::runtime_error("sampleB_bayesr_ST_csr: invalid prior_scale value.");
   }
   ssb_scaled += b(iu) * b(iu) / (ck * scale_i);
   dfb += 1.0;
  }
 }

 const double scale = ssb_scaled + nub * ssb_prior;
 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleB_bayesr_ST_csr: invalid scale.");
 }

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 vb = std::max(scale / chi2, 1e-12);
}

inline void sampleB_bayesr_ST_csr_unscaled(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  const arma::vec& mixture_var,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb_scaled = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const int k = comp(iu);

  if (k > 0) {
   if (k >= static_cast<int>(mixture_var.n_elem)) {
    throw std::runtime_error("sampleB_bayesr_ST_csr_unscaled: component index out of range.");
   }
   const double ck = mixture_var(static_cast<arma::uword>(k));
   if (!std::isfinite(ck) || ck <= 0.0) {
    throw std::runtime_error("sampleB_bayesr_ST_csr_unscaled: active component has invalid mixture_var.");
   }
   ssb_scaled += b(iu) * b(iu) / ck;
   dfb += 1.0;
  }
 }

 const double scale = ssb_scaled + nub * ssb_prior;
 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleB_bayesr_ST_csr_unscaled: invalid scale.");
 }

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 vb = std::max(scale / chi2, 1e-12);
}

inline void samplePi_bayesr_ST_csr(
  const arma::Row<int>& comp,
  std::vector<double>& pi,
  const std::vector<double>& alpha,
  std::mt19937& gen
) {
 const int K = static_cast<int>(pi.size());
 std::vector<double> counts(static_cast<std::size_t>(K), 0.0);

 for (arma::uword i = 0; i < comp.n_elem; ++i) {
  const int k = comp(i);
  if (k < 0 || k >= K) {
   throw std::runtime_error("samplePi_bayesr_ST_csr: component index out of range.");
  }
  counts[static_cast<std::size_t>(k)] += 1.0;
 }

 double sum_g = 0.0;
 for (int k = 0; k < K; ++k) {
  const double shape = alpha[static_cast<std::size_t>(k)] + counts[static_cast<std::size_t>(k)];
  if (!std::isfinite(shape) || shape <= 0.0) {
   throw std::runtime_error("samplePi_bayesr_ST_csr: invalid Dirichlet shape.");
  }
  std::gamma_distribution<double> rg(shape, 1.0);
  pi[static_cast<std::size_t>(k)] = std::max(rg(gen), 1e-300);
  sum_g += pi[static_cast<std::size_t>(k)];
 }

 if (!std::isfinite(sum_g) || sum_g <= 0.0) {
  throw std::runtime_error("samplePi_bayesr_ST_csr: invalid gamma sum.");
 }

 for (int k = 0; k < K; ++k) pi[static_cast<std::size_t>(k)] /= sum_g;
}

inline void normalize_pi_bayesr(std::vector<double>& pi) {
 double s = 0.0;
 for (double p : pi) {
  if (!std::isfinite(p) || p < 0.0) {
   throw std::runtime_error("normalize_pi_bayesr: pi must be finite and non-negative.");
  }
  s += p;
 }
 if (!std::isfinite(s) || s <= 0.0) {
  throw std::runtime_error("normalize_pi_bayesr: pi must have positive sum.");
 }
 for (double& p : pi) p /= s;
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_bayesr(
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
  std::vector<double> pi,
  std::vector<double> mixture_var,
  std::vector<double> alpha,
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
  int nchains,
  bool keep_chains,
  std::vector<int> chain_seeds,
  int updateE_start = 0,
  int updateE_every = 1,
  bool updateLDswap = false,
  double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8,
  int ld_swap_max_friends = 50,
  int ld_swap_moves = 1,
  Rcpp::Nullable<Rcpp::NumericVector> selection_s_prior_scale = R_NilValue
) {
 const int nt = static_cast<int>(wy.size());
 if (nt <= 0) throw std::runtime_error("stblr_cpg_omp_csr_bayesr: nt must be positive.");

 const int m = static_cast<int>(wy[0].size());
 if (m <= 0) throw std::runtime_error("stblr_cpg_omp_csr_bayesr: m must be positive.");
 if (nit <= 0 || nburn < 0 || nthin <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: invalid nit/nburn/nthin.");
 }
 if (nchains <= 0) throw std::runtime_error("stblr_cpg_omp_csr_bayesr: nchains must be positive.");
 if (updateE_start < 0) throw std::runtime_error("stblr_cpg_omp_csr_bayesr: updateE_start must be non-negative.");
 if (updateE_every <= 0) throw std::runtime_error("stblr_cpg_omp_csr_bayesr: updateE_every must be positive.");
 if (!std::isfinite(ld_swap_prob) || ld_swap_prob < 0.0 || ld_swap_prob > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: ld_swap_prob must be in [0, 1].");
 }
 if (!std::isfinite(ld_swap_r2) || ld_swap_r2 < 0.0 || ld_swap_r2 > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: ld_swap_r2 must be in [0, 1].");
 }
 if (ld_swap_max_friends <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: ld_swap_max_friends must be positive.");
 }
 if (ld_swap_moves < 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: ld_swap_moves must be non-negative.");
 }
 if (!chain_seeds.empty() && static_cast<int>(chain_seeds.size()) != nchains) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: chain_seeds must have length nchains.");
 }

 bool use_selection_s_prior_scale = selection_s_prior_scale.isNotNull();
 Rcpp::NumericVector selection_s_prior_scale_vec;
 if (use_selection_s_prior_scale) {
  selection_s_prior_scale_vec = Rcpp::NumericVector(selection_s_prior_scale);
  use_selection_s_prior_scale = selection_s_prior_scale_vec.size() > 0;
 }

 if (use_selection_s_prior_scale &&
     static_cast<int>(selection_s_prior_scale_vec.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: selection_s_prior_scale must have length m.");
 }

 const int K = static_cast<int>(mixture_var.size());
 if (K < 2) throw std::runtime_error("stblr_cpg_omp_csr_bayesr: mixture_var must have at least two components.");
 if (static_cast<int>(pi.size()) != K || static_cast<int>(alpha.size()) != K) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: pi, alpha, and mixture_var must have equal length.");
 }
 if (!std::isfinite(mixture_var[0]) || mixture_var[0] != 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: mixture_var[0] must be 0 for the null component.");
 }
 for (int k = 1; k < K; ++k) {
  if (!std::isfinite(mixture_var[static_cast<std::size_t>(k)]) ||
      mixture_var[static_cast<std::size_t>(k)] <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: non-null mixture_var values must be positive.");
  }
 }
 for (double a : alpha) {
  if (!std::isfinite(a) || a <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: alpha must be finite and positive.");
  }
 }
 normalize_pi_bayesr(pi);

 if ((int)ww.size() != nt || (int)b_init.size() != nt ||
     (int)yy.size() != nt || (int)n.size() != nt ||
     (int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: inconsistent trait dimensions.");
 }
 if ((int)B.n_rows != nt || (int)B.n_cols != nt ||
     (int)E.n_rows != nt || (int)E.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: B and E must be nt x nt.");
 }
 if (use_r_init && static_cast<int>(r_init.size()) != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: r_init must have length nt when enabled.");
 }
 if (use_comp_init && static_cast<int>(comp_init.size()) != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: comp_init must have length nt when enabled.");
 }

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_init_mat(nt, m, arma::fill::zeros);
 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
 arma::vec yy_vec(nt, arma::fill::zeros);
 arma::rowvec prior_scale;

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m || (int)ww[t].size() != m ||
      (int)b_init[t].size() != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: inconsistent marker dimensions.");
  }
  if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: priors must be nt x nt.");
  }
  if (use_r_init && static_cast<int>(r_init[t].size()) != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: each r_init[t] must have length m.");
  }
  if (use_comp_init && static_cast<int>(comp_init[t].size()) != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: each comp_init[t] must have length m.");
  }

  yy_vec(static_cast<arma::uword>(t)) = yy[static_cast<std::size_t>(t)];
  for (int i = 0; i < m; ++i) {
   wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
   ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
   b_init_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = b_init[t][i];
  }
  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
   sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
  }
 }

 if (use_selection_s_prior_scale) {
  prior_scale.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double scale_i = selection_s_prior_scale_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(scale_i) || scale_i <= 0.0) {
    throw std::runtime_error(
     "stblr_cpg_omp_csr_bayesr: selection_s_prior_scale must contain positive finite values."
    );
   }
   prior_scale(static_cast<arma::uword>(i)) = scale_i;
  }
 }

 for (int t = 1; t < nt; ++t) {
  if (n[t] != n[0]) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: shared-LD scaling assumes equal n across traits.");
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m));
 for (int i = 0; i < m; ++i) {
  const double w0 = ww_mat(0, static_cast<arma::uword>(i));
  if (!std::isfinite(w0) || w0 <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: ww contains invalid value.");
  }
  xx[static_cast<std::size_t>(i)] = w0;
  for (int t = 1; t < nt; ++t) {
   const double wt = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   const double tol = 1e-8 * std::max(1.0, std::abs(w0));
   if (!std::isfinite(wt) || wt <= 0.0 || std::abs(wt - w0) > tol) {
    throw std::runtime_error("stblr_cpg_omp_csr_bayesr: ww differs across traits.");
   }
  }
 }

 STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
 BayesRLDLDFriends ld_swap_friends;
 if (updateLDswap) {
  ld_swap_friends = build_ld_swap_friends_bayesr_ST_csr(
   m,
   ld,
   xx,
   ld_swap_r2,
   ld_swap_max_friends
  );
 } else {
  ld_swap_friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
 }

 arma::vec mixture_var_vec(K, arma::fill::zeros);
 for (int k = 0; k < K; ++k) mixture_var_vec(static_cast<arma::uword>(k)) = mixture_var[k];

 std::vector<int> order(static_cast<std::size_t>(m));
 std::iota(order.begin(), order.end(), 0);

 const int ntasks = stblr_num_chain_tasks(nt, nchains);
 const int nthreads = stblr_num_threads_for_tasks(ncores, ntasks);
 const int trace_len = nit + nburn;

 arma::mat bm_task(ntasks, m, arma::fill::zeros);
 arma::mat dm_task(ntasks, m, arma::fill::zeros);
 arma::mat component_mean_task(ntasks, m, arma::fill::zeros);
 arma::mat b_task(ntasks, m, arma::fill::zeros);
 arma::mat r_task(ntasks, m, arma::fill::zeros);
 arma::mat component_task(ntasks, m, arma::fill::zeros);
 arma::mat vbs_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat vgs_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat ves_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat vles_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat vlds_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat final_pi_task(ntasks, K, arma::fill::zeros);
 arma::mat mean_pi_task(ntasks, K, arma::fill::zeros);
 arma::vec final_vb_task(ntasks, arma::fill::zeros);
 arma::vec final_vg_task(ntasks, arma::fill::zeros);
 arma::vec final_ve_task(ntasks, arma::fill::zeros);
 arma::vec min_sse_task(ntasks, arma::fill::zeros);
 arma::ivec min_sse_iter_task(ntasks, arma::fill::value(-1));
 arma::vec min_residual_scale_task(ntasks, arma::fill::zeros);
 arma::ivec max_nonzero_components_task(ntasks, arma::fill::zeros);
 arma::vec max_abs_effect_task(ntasks, arma::fill::zeros);
 arma::vec max_fitted_quadratic_task(ntasks, arma::fill::zeros);
 arma::ivec n_updateE_task(ntasks, arma::fill::zeros);
 arma::vec ld_swap_attempted_task(ntasks, arma::fill::zeros);
 arma::vec ld_swap_accepted_task(ntasks, arma::fill::zeros);
 std::vector<arma::mat> comp_prob_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> ncomp_task(static_cast<std::size_t>(ntasks));
 std::vector<int> failed(static_cast<std::size_t>(ntasks), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(ntasks));

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int task = 0; task < ntasks; ++task) {
  const int t = stblr_task_trait(task, nchains);
  const int chain = stblr_task_chain(task, nchains);
  const arma::uword task_u = static_cast<arma::uword>(task);

  try {
   unsigned int task_seed = 0u;
   if (!chain_seeds.empty()) {
    task_seed = stblr_seed_with_chain_base(chain_seeds[static_cast<std::size_t>(chain)], t);
   } else if (nchains == 1) {
    task_seed = stblr_trait_seed(seed, t);
   } else {
    task_seed = stblr_chain_seed(seed, t, chain);
   }

   std::mt19937 gen_t(task_seed);
   std::vector<int> order_t = order;
   std::shuffle(order_t.begin(), order_t.end(), gen_t);

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));
   arma::rowvec b_t = b_init_mat.row(static_cast<arma::uword>(t));
   arma::rowvec r_t(m, arma::fill::zeros);
   arma::Row<int> comp_t(m, arma::fill::zeros);

   if (use_comp_init) {
    for (int i = 0; i < m; ++i) {
     const double ci = comp_init[t][i];
     if (!std::isfinite(ci) || ci != std::floor(ci) || ci < 0.0 || ci >= K) {
      throw std::runtime_error("comp_init contains invalid component index.");
     }
     comp_t(static_cast<arma::uword>(i)) = static_cast<int>(ci);
    }
   } else {
    for (int i = 0; i < m; ++i) {
     comp_t(static_cast<arma::uword>(i)) = (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
    }
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) {
     r_t(static_cast<arma::uword>(i)) = r_init[t][i];
    }
    if (!r_t.is_finite()) throw std::runtime_error("r_init contains NaN/Inf.");
   } else {
    rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   if (!std::isfinite(vb_t) || vb_t <= 0.0 || !std::isfinite(ve_t) || ve_t <= 0.0) {
    throw std::runtime_error("initial B/E diagonal values must be finite and positive.");
   }

   std::vector<double> pi_t = pi;
   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec component_mean_t(m, arma::fill::zeros);
   arma::mat comp_prob_t(m, K, arma::fill::zeros);
   arma::vec pi_mean_t(K, arma::fill::zeros);
   double nsamples_t = 0.0;
   double min_sse_t = std::numeric_limits<double>::infinity();
   int min_sse_iter_t = -1;
   double min_residual_scale_t = std::numeric_limits<double>::infinity();
   int max_nonzero_components_t = 0;
   double max_abs_effect_t = 0.0;
   double max_fitted_quadratic_t = 0.0;
   int n_updateE_t = 0;
   double ld_swap_attempted_t = 0.0;
   double ld_swap_accepted_t = 0.0;

   double vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_bayesr_ST_csr(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   for (int it = 0; it < trace_len; ++it) {
    if (use_selection_s_prior_scale) {
     for (int isort = 0; isort < m; ++isort) {
      sampleBetaR_ST_csr(
       order_t[static_cast<std::size_t>(isort)],
       pi_t,
       mixture_var_vec,
       vb_t,
       prior_scale,
       vei_t,
       ww_t,
       r_t,
       b_t,
       comp_t,
       ld,
       gen_t
      );
     }
    } else {
     for (int isort = 0; isort < m; ++isort) {
      sampleBetaR_ST_csr_unscaled(
       order_t[static_cast<std::size_t>(isort)],
       pi_t,
       mixture_var_vec,
       vb_t,
       vei_t,
       ww_t,
       r_t,
       b_t,
       comp_t,
       ld,
       gen_t
      );
     }
    }

    if (updateLDswap && ld_swap_moves > 0 && ld_swap_prob > 0.0) {
     std::uniform_real_distribution<double> runif(0.0, 1.0);

     if (runif(gen_t) < ld_swap_prob) {
      for (int move = 0; move < ld_swap_moves; ++move) {
       bool attempted = false;
       const bool accepted = use_selection_s_prior_scale ?
        attempt_ld_swap_bayesr_ST_csr(
         m,
         t,
         chain,
         it,
         vei_t,
         vb_t,
         yy_vec(static_cast<arma::uword>(t)),
         ww_t,
         wy_t,
         mixture_var_vec,
         prior_scale,
         r_t,
         b_t,
         comp_t,
         ld,
         ld_swap_friends,
         gen_t,
         attempted
        ) :
        attempt_ld_swap_bayesr_ST_csr_unscaled(
         m,
         t,
         chain,
         it,
         vei_t,
         yy_vec(static_cast<arma::uword>(t)),
         ww_t,
         wy_t,
         r_t,
         b_t,
         comp_t,
         ld,
         ld_swap_friends,
         gen_t,
         attempted
        );
       if (attempted) {
        ld_swap_attempted_t += 1.0;
        if (accepted) ld_swap_accepted_t += 1.0;
       }
      }
     }
    }

    if (updateB) {
     if (use_selection_s_prior_scale) {
      sampleB_bayesr_ST_csr(
       m,
       nub,
       vb_t,
       b_t,
       comp_t,
       mixture_var_vec,
       prior_scale,
       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
       gen_t
      );
     } else {
      sampleB_bayesr_ST_csr_unscaled(
       m,
       nub,
       vb_t,
       b_t,
       comp_t,
       mixture_var_vec,
       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
       gen_t
      );
     }
    }

    const bool do_updateE =
     updateE &&
     it >= updateE_start &&
     ((it - updateE_start) % updateE_every == 0);

    if (do_updateE) {
     ensure_null_effects_bayesr_ST_csr(m, t, chain, it, b_t, comp_t);
     rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
     const BayesRUpdateEDiagnostics diag = residual_diagnostics_bayesr_ST_csr(
      m,
      nue,
      b_t,
      r_t,
      comp_t,
      wy_t,
      sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      yy_vec(static_cast<arma::uword>(t))
     );
     if (diag.sse < min_sse_t) {
      min_sse_t = diag.sse;
      min_sse_iter_t = it;
     }
     min_residual_scale_t = std::min(min_residual_scale_t, diag.residual_scale);
     max_nonzero_components_t = std::max(max_nonzero_components_t, diag.nonzero_components);
     max_abs_effect_t = std::max(max_abs_effect_t, diag.max_abs_b);
     max_fitted_quadratic_t = std::max(max_fitted_quadratic_t, std::abs(diag.bXb));
     ++n_updateE_t;
     check_residual_scale_bayesr_ST_csr(
      m,
      t,
      chain,
      it,
      nue,
      ve_t,
      vb_t,
      mixture_var_vec,
      b_t,
      r_t,
      comp_t,
      wy_t,
      sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      yy_vec(static_cast<arma::uword>(t)),
      n[t],
      adjE
     );
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

    if (updatePi) samplePi_bayesr_ST_csr(comp_t, pi_t, alpha, gen_t);

    vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
    vle_t = computeLE_bayesr_ST_csr(m, b_t, ww_t, n[t]);
    vld_t = vg_t - vle_t;
    vei_t = ve_t + adjE * vg_t;
    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error("adjusted residual variance became invalid.");
    }

    vbs_task(task_u, static_cast<arma::uword>(it)) = vb_t;
    vgs_task(task_u, static_cast<arma::uword>(it)) = vg_t;
    ves_task(task_u, static_cast<arma::uword>(it)) = ve_t;
    vles_task(task_u, static_cast<arma::uword>(it)) = vle_t;
    vlds_task(task_u, static_cast<arma::uword>(it)) = vld_t;

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;
     for (int k = 0; k < K; ++k) {
      pi_mean_t(static_cast<arma::uword>(k)) += pi_t[static_cast<std::size_t>(k)];
     }
     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      const int ci = comp_t(iu);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += ci > 0 ? 1.0 : 0.0;
      component_mean_t(iu) += static_cast<double>(ci);
      comp_prob_t(iu, static_cast<arma::uword>(ci)) += 1.0;
     }
    }
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;
   bm_t /= nsamples_t;
   dm_t /= nsamples_t;
   component_mean_t /= nsamples_t;
   comp_prob_t /= nsamples_t;
   pi_mean_t /= nsamples_t;

   bm_task.row(task_u) = bm_t;
   dm_task.row(task_u) = dm_t;
   component_mean_task.row(task_u) = component_mean_t;
   b_task.row(task_u) = b_t;
   r_task.row(task_u) = r_t;
   for (int i = 0; i < m; ++i) {
    component_task(task_u, static_cast<arma::uword>(i)) =
     static_cast<double>(comp_t(static_cast<arma::uword>(i)));
   }
   for (int k = 0; k < K; ++k) {
    final_pi_task(task_u, static_cast<arma::uword>(k)) = pi_t[static_cast<std::size_t>(k)];
    mean_pi_task(task_u, static_cast<arma::uword>(k)) = pi_mean_t(static_cast<arma::uword>(k));
   }
   final_vb_task(task_u) = vb_t;
   final_vg_task(task_u) = vg_t;
   final_ve_task(task_u) = ve_t;
   min_sse_task(task_u) = std::isfinite(min_sse_t) ? min_sse_t : NA_REAL;
   min_sse_iter_task(task_u) = min_sse_iter_t;
   min_residual_scale_task(task_u) = std::isfinite(min_residual_scale_t) ? min_residual_scale_t : NA_REAL;
   max_nonzero_components_task(task_u) = max_nonzero_components_t;
   max_abs_effect_task(task_u) = max_abs_effect_t;
   max_fitted_quadratic_task(task_u) = max_fitted_quadratic_t;
   n_updateE_task(task_u) = n_updateE_t;
   ld_swap_attempted_task(task_u) = ld_swap_attempted_t;
   ld_swap_accepted_task(task_u) = ld_swap_accepted_t;
   comp_prob_task[static_cast<std::size_t>(task)] = comp_prob_t;
   ncomp_task[static_cast<std::size_t>(task)] = arma::sum(comp_prob_t, 0).t();
  } catch (const std::exception& e) {
   failed[static_cast<std::size_t>(task)] = 1;
   errors[static_cast<std::size_t>(task)] = e.what();
  } catch (...) {
   failed[static_cast<std::size_t>(task)] = 1;
   errors[static_cast<std::size_t>(task)] = "unknown error";
  }
 }

 for (int task = 0; task < ntasks; ++task) {
  if (failed[static_cast<std::size_t>(task)]) {
   throw std::runtime_error(
    "stblr_cpg_omp_csr_bayesr failed for trait " +
    std::to_string(stblr_task_trait(task, nchains)) +
    ", chain " + std::to_string(stblr_task_chain(task, nchains)) +
    ": " + errors[static_cast<std::size_t>(task)]
   );
  }
 }

 const double inv_chains = 1.0 / static_cast<double>(nchains);
 arma::mat bm(nt, m, arma::fill::zeros);
 arma::mat dm(nt, m, arma::fill::zeros);
 arma::mat bm_sd(nt, m, arma::fill::zeros);
 arma::mat dm_sd(nt, m, arma::fill::zeros);
 arma::mat bm_min(nt, m, arma::fill::zeros);
 arma::mat dm_min(nt, m, arma::fill::zeros);
 arma::mat bm_max(nt, m, arma::fill::zeros);
 arma::mat dm_max(nt, m, arma::fill::zeros);
 arma::mat component_mean(nt, m, arma::fill::zeros);
 arma::mat b_out(nt, m, arma::fill::zeros);
 arma::mat r_out(nt, m, arma::fill::zeros);
 arma::mat component_out(nt, m, arma::fill::zeros);
 arma::mat vbs(nt, trace_len, arma::fill::zeros);
 arma::mat vgs(nt, trace_len, arma::fill::zeros);
 arma::mat ves(nt, trace_len, arma::fill::zeros);
 arma::mat vle(nt, trace_len, arma::fill::zeros);
 arma::mat vld(nt, trace_len, arma::fill::zeros);
 arma::mat final_pi(nt, K, arma::fill::zeros);
 arma::mat mean_pi(nt, K, arma::fill::zeros);
 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::mat updateE_diagnostics(ntasks, 9, arma::fill::zeros);
 arma::mat ld_swap_diagnostics(nt, 3, arma::fill::zeros);
 arma::mat ld_swap_chain_diagnostics(ntasks, 5, arma::fill::zeros);
 std::vector<arma::mat> comp_prob(static_cast<std::size_t>(nt));
 arma::mat ncomp(nt, K, arma::fill::zeros);

 for (int task = 0; task < ntasks; ++task) {
  const arma::uword task_u = static_cast<arma::uword>(task);
  const int task_trait = stblr_task_trait(task, nchains);
  const int task_chain = stblr_task_chain(task, nchains);
  updateE_diagnostics(task_u, 0) = stblr_task_trait(task, nchains);
  updateE_diagnostics(task_u, 1) = stblr_task_chain(task, nchains);
  updateE_diagnostics(task_u, 2) = n_updateE_task(task_u);
  updateE_diagnostics(task_u, 3) = min_sse_task(task_u);
  updateE_diagnostics(task_u, 4) = min_sse_iter_task(task_u);
  updateE_diagnostics(task_u, 5) = min_residual_scale_task(task_u);
  updateE_diagnostics(task_u, 6) = max_nonzero_components_task(task_u);
  updateE_diagnostics(task_u, 7) = max_abs_effect_task(task_u);
  updateE_diagnostics(task_u, 8) = max_fitted_quadratic_task(task_u);

  const double attempted = ld_swap_attempted_task(task_u);
  const double accepted = ld_swap_accepted_task(task_u);
  ld_swap_chain_diagnostics(task_u, 0) = task_trait;
  ld_swap_chain_diagnostics(task_u, 1) = task_chain;
  ld_swap_chain_diagnostics(task_u, 2) = attempted;
  ld_swap_chain_diagnostics(task_u, 3) = accepted;
  ld_swap_chain_diagnostics(task_u, 4) =
   attempted > 0.0 ? accepted / attempted : 0.0;
  ld_swap_diagnostics(static_cast<arma::uword>(task_trait), 0) += attempted;
  ld_swap_diagnostics(static_cast<arma::uword>(task_trait), 1) += accepted;
 }

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  ld_swap_diagnostics(tu, 2) =
   ld_swap_diagnostics(tu, 0) > 0.0
   ? ld_swap_diagnostics(tu, 1) / ld_swap_diagnostics(tu, 0)
   : 0.0;
 }

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  bm_min.row(tu).fill(std::numeric_limits<double>::infinity());
  dm_min.row(tu).fill(std::numeric_limits<double>::infinity());
  bm_max.row(tu).fill(-std::numeric_limits<double>::infinity());
  dm_max.row(tu).fill(-std::numeric_limits<double>::infinity());
  comp_prob[static_cast<std::size_t>(t)] = arma::mat(m, K, arma::fill::zeros);

  for (int chain = 0; chain < nchains; ++chain) {
   const int task = t * nchains + chain;
   const arma::uword task_u = static_cast<arma::uword>(task);
   bm.row(tu) += bm_task.row(task_u);
   dm.row(tu) += dm_task.row(task_u);
   component_mean.row(tu) += component_mean_task.row(task_u);
   b_out.row(tu) += b_task.row(task_u);
   r_out.row(tu) += r_task.row(task_u);
   component_out.row(tu) += component_task.row(task_u);
   vbs.row(tu) += vbs_task.row(task_u);
   vgs.row(tu) += vgs_task.row(task_u);
   ves.row(tu) += ves_task.row(task_u);
   vle.row(tu) += vles_task.row(task_u);
   vld.row(tu) += vlds_task.row(task_u);
   final_pi.row(tu) += final_pi_task.row(task_u);
   mean_pi.row(tu) += mean_pi_task.row(task_u);
   final_vb(tu) += final_vb_task(task_u);
   final_vg(tu) += final_vg_task(task_u);
   final_ve(tu) += final_ve_task(task_u);
   comp_prob[static_cast<std::size_t>(t)] += comp_prob_task[static_cast<std::size_t>(task)];
   ncomp.row(tu) += ncomp_task[static_cast<std::size_t>(task)].t();

   for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    bm_min(tu, iu) = std::min(bm_min(tu, iu), bm_task(task_u, iu));
    dm_min(tu, iu) = std::min(dm_min(tu, iu), dm_task(task_u, iu));
    bm_max(tu, iu) = std::max(bm_max(tu, iu), bm_task(task_u, iu));
    dm_max(tu, iu) = std::max(dm_max(tu, iu), dm_task(task_u, iu));
   }
  }

  bm.row(tu) *= inv_chains;
  dm.row(tu) *= inv_chains;
  component_mean.row(tu) *= inv_chains;
  b_out.row(tu) *= inv_chains;
  r_out.row(tu) *= inv_chains;
  component_out.row(tu) *= inv_chains;
  vbs.row(tu) *= inv_chains;
  vgs.row(tu) *= inv_chains;
  ves.row(tu) *= inv_chains;
  vle.row(tu) *= inv_chains;
  vld.row(tu) *= inv_chains;
  final_pi.row(tu) *= inv_chains;
  mean_pi.row(tu) *= inv_chains;
  final_vb(tu) *= inv_chains;
  final_vg(tu) *= inv_chains;
  final_ve(tu) *= inv_chains;
  comp_prob[static_cast<std::size_t>(t)] *= inv_chains;
  ncomp.row(tu) *= inv_chains;

  if (nchains > 1) {
   for (int chain = 0; chain < nchains; ++chain) {
    const int task = t * nchains + chain;
    const arma::uword task_u = static_cast<arma::uword>(task);
    arma::rowvec bm_diff = bm_task.row(task_u) - bm.row(tu);
    arma::rowvec dm_diff = dm_task.row(task_u) - dm.row(tu);
    bm_sd.row(tu) += bm_diff % bm_diff;
    dm_sd.row(tu) += dm_diff % dm_diff;
   }
   bm_sd.row(tu) = arma::sqrt(bm_sd.row(tu) / static_cast<double>(nchains - 1));
   dm_sd.row(tu) = arma::sqrt(dm_sd.row(tu) / static_cast<double>(nchains - 1));
  }
 }

 Rcpp::List comp_prob_out(nt);
 for (int t = 0; t < nt; ++t) comp_prob_out[t] = comp_prob[static_cast<std::size_t>(t)];

 Rcpp::List chains_out;
 if (keep_chains) {
  chains_out = Rcpp::List(nt);
  for (int t = 0; t < nt; ++t) {
   Rcpp::List trait_chains(nchains);
   for (int chain = 0; chain < nchains; ++chain) {
    const int task = t * nchains + chain;
    const arma::uword task_u = static_cast<arma::uword>(task);
    Rcpp::NumericVector updateE_diag(9);
    for (int j = 0; j < 9; ++j) {
     updateE_diag[j] = updateE_diagnostics(task_u, static_cast<arma::uword>(j));
    }
    Rcpp::NumericVector ld_swap_diag(3);
    ld_swap_diag[0] = ld_swap_chain_diagnostics(task_u, 2);
    ld_swap_diag[1] = ld_swap_chain_diagnostics(task_u, 3);
    ld_swap_diag[2] = ld_swap_chain_diagnostics(task_u, 4);
    trait_chains[chain] = Rcpp::List::create(
     Rcpp::Named("dm") = dm_task.row(task_u).t(),
     Rcpp::Named("bm") = bm_task.row(task_u).t(),
     Rcpp::Named("comp_prob") = comp_prob_task[static_cast<std::size_t>(task)],
     Rcpp::Named("dm_component_mean") = component_mean_task.row(task_u).t(),
     Rcpp::Named("final_pi") = final_pi_task.row(task_u).t(),
     Rcpp::Named("mean_pi") = mean_pi_task.row(task_u).t(),
     Rcpp::Named("vbs") = vbs_task.row(task_u).t(),
     Rcpp::Named("vgs") = vgs_task.row(task_u).t(),
     Rcpp::Named("ves") = ves_task.row(task_u).t(),
     Rcpp::Named("vle") = vles_task.row(task_u).t(),
     Rcpp::Named("vld") = vlds_task.row(task_u).t(),
     Rcpp::Named("updateE_diagnostics") = updateE_diag,
     Rcpp::Named("ld_swap") = ld_swap_diag
    );
   }
   chains_out[t] = trait_chains;
  }
 }

 arma::mat covb(nt, nt, arma::fill::zeros);
 arma::mat covg(nt, nt, arma::fill::zeros);
 arma::mat cove(nt, nt, arma::fill::zeros);
 arma::mat vb(nt, nt, arma::fill::zeros);
 arma::mat vg(nt, nt, arma::fill::zeros);
 arma::mat ve(nt, nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  covb(tu, tu) = final_vb(tu);
  covg(tu, tu) = final_vg(tu);
  cove(tu, tu) = final_ve(tu);
  vb(tu, tu) = final_vb(tu);
  vg(tu, tu) = final_vg(tu);
  ve(tu, tu) = final_ve(tu);
 }

 Rcpp::List out = Rcpp::List::create(
  Rcpp::Named("bm") = bm.t(),
  Rcpp::Named("dm") = dm.t(),
  Rcpp::Named("wy") = wy_mat.t(),
  Rcpp::Named("r") = r_out.t(),
  Rcpp::Named("b") = b_out.t(),
  Rcpp::Named("component") = component_out.t(),
  Rcpp::Named("marker_index") = arma::regspace<arma::vec>(0, m - 1),
  Rcpp::Named("vbs") = vbs.t(),
  Rcpp::Named("vgs") = vgs.t(),
  Rcpp::Named("ves") = ves.t(),
  Rcpp::Named("covb") = covb,
  Rcpp::Named("covg") = covg,
  Rcpp::Named("cove") = cove,
  Rcpp::Named("vb") = vb,
  Rcpp::Named("vg") = vg,
  Rcpp::Named("ve") = ve,
  Rcpp::Named("pi") = final_pi,
  Rcpp::Named("pim") = mean_pi,
  Rcpp::Named("vle") = vle.t(),
  Rcpp::Named("vld") = vld.t(),
  Rcpp::Named("bm_sd") = bm_sd.t(),
  Rcpp::Named("bm_min") = bm_min.t(),
  Rcpp::Named("bm_max") = bm_max.t(),
  Rcpp::Named("dm_sd") = dm_sd.t(),
  Rcpp::Named("dm_min") = dm_min.t(),
  Rcpp::Named("dm_max") = dm_max.t(),
  Rcpp::Named("comp_prob") = comp_prob_out,
  Rcpp::Named("dm_component_mean") = component_mean.t(),
  Rcpp::Named("ncomp") = ncomp,
  Rcpp::Named("mixture_var") = mixture_var_vec,
  Rcpp::Named("updateE_diagnostics") = updateE_diagnostics,
  Rcpp::Named("ld_swap") = ld_swap_diagnostics,
  Rcpp::Named("ld_swap_chains") = ld_swap_chain_diagnostics
 );
 if (keep_chains) out["chains"] = chains_out;
 return out;
}
