// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "st_chain_utils.h"
#include "st_block_eigen.h"
#include "st_block_eigen_rcpp.h"
#include "st_csr_common.h"
#include "st_ld_operator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
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

static std::vector<std::string> stblr_bayesr_copy_character_vector(
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

static std::vector<std::vector<int>> stblr_bayesr_copy_int_list(
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

static std::vector<int> stblr_bayesr_copy_rows0_or_empty(
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  int n_bed
) {
 std::vector<int> out;
 if (rows.isNull()) return out;
 Rcpp::IntegerVector r(rows);
 out.resize(static_cast<std::size_t>(r.size()));
 bool identity_rows = (r.size() == n_bed);
 for (int i = 0; i < r.size(); ++i) {
  if (r[i] == NA_INTEGER || r[i] < 1 || r[i] > n_bed) {
   throw std::runtime_error("rows contains an index outside [1, n_bed].");
  }
  if (r[i] != i + 1) identity_rows = false;
  out[static_cast<std::size_t>(i)] = r[i] - 1;
 }
 if (identity_rows) out.clear();
 return out;
}


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

template <class OpT>
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
  const OpT& op,
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

  op.apply_offdiag(i, diff, r);
 }

 b(iu) = b_new;
 comp(iu) = k_new;
}

template <class OpT>
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
  const OpT& op,
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

  op.apply_offdiag(i, diff, r);
 }

 b(iu) = b_new;
 comp(iu) = k_new;
}

inline double bayesr_maf_effect_scale(double s, double log_h) {
 return std::exp((s + 1.0) * log_h);
}

inline void fill_maf_effect_s_prior_scale_bayesr(
  int m,
  double s,
  const arma::rowvec& log_h,
  arma::rowvec& prior_scale
) {
 prior_scale.set_size(static_cast<arma::uword>(m));
 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double scale_i = bayesr_maf_effect_scale(s, log_h(iu));
  if (!std::isfinite(scale_i) || scale_i <= 0.0) {
   throw std::runtime_error("dynamic BayesR maf_effect_s prior scale became invalid.");
  }
  prior_scale(iu) = scale_i;
 }
}

struct BayesRLDLDFriends {
 std::vector<uint64_t> ptr;
 std::vector<int> idx;
 std::vector<double> r2;
};

template <class OpT>
struct BayesROperatorContext {
 OpT op;
 BayesRLDLDFriends ld_swap_friends;
 Rcpp::List diagnostics;

 BayesROperatorContext(
   OpT&& op_,
   BayesRLDLDFriends&& ld_swap_friends_,
   const Rcpp::List& diagnostics_
 ) :
  op(std::move(op_)),
  ld_swap_friends(std::move(ld_swap_friends_)),
  diagnostics(diagnostics_) {}
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

template <class OpT>
inline void set_marker_state_bayesr_ST_csr(
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

template <class OpT>
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
  const OpT& op,
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
 set_marker_state_bayesr_ST_csr(j, 0.0, 0, ww, r, b, comp, op);
 set_marker_state_bayesr_ST_csr(k, b_j_old, comp_j_old, ww, r, b, comp, op);

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
 // cancel the prior density. Fixed maf_effect_s makes the moved component's
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

template <class OpT>
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
  const OpT& op,
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
 set_marker_state_bayesr_ST_csr(j, 0.0, 0, ww, r, b, comp, op);
 set_marker_state_bayesr_ST_csr(k, b_j_old, comp_j_old, ww, r, b, comp, op);

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

inline double logpost_maf_effect_s_bayesr(
  double s,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  double vb,
  const arma::vec& mixture_var,
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
  const int k = comp(j);
  if (k <= 0) continue;
  if (k >= static_cast<int>(mixture_var.n_elem)) {
   return -std::numeric_limits<double>::infinity();
  }
  const double ck = mixture_var(static_cast<arma::uword>(k));
  if (!std::isfinite(ck) || ck <= 0.0) {
   return -std::numeric_limits<double>::infinity();
  }
  const double log_q = (s + 1.0) * log_h(j);
  const double q = std::exp(log_q);
  if (!std::isfinite(q) || q <= 0.0) {
   return -std::numeric_limits<double>::infinity();
  }
  lp += -0.5 * (log_q + b(j) * b(j) / (vb_safe * ck * q));
 }
 return lp;
}

inline bool update_maf_effect_s_bayesr(
  double& maf_effect_s_current,
  const arma::rowvec& b,
  const arma::Row<int>& comp,
  double vb,
  const arma::vec& mixture_var,
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

 const double lp_current = logpost_maf_effect_s_bayesr(
  maf_effect_s_current, b, comp, vb, mixture_var, log_h, prior_lower, prior_upper
 );
 const double lp_prop = logpost_maf_effect_s_bayesr(
  maf_effect_s_prop, b, comp, vb, mixture_var, log_h, prior_lower, prior_upper
 );
 const double log_alpha = lp_prop - lp_current;
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 if (std::log(std::max(runif(gen), 1e-300)) < log_alpha) {
  maf_effect_s_current = maf_effect_s_prop;
  return true;
 }
 return false;
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

#define SBLR_CSR_BAYESR_CORE_IMPL_TRANSLATION_UNIT 1
#include "blr_csr_bayesr_core_impl.h"
#undef SBLR_CSR_BAYESR_CORE_IMPL_TRANSLATION_UNIT

template <class MakeOperator>
Rcpp::List stblr_cpg_omp_csr_bayesr_impl(
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
  int updateE_start,
  int updateE_every,
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
  bool convergence_probability,
  bool convergence_b,
  bool convergence_d,
  bool convergence_component,
  MakeOperator make_operator
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
 for (int marker : convergence_markers) if (marker < 0 || marker >= m) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: convergence marker index is out of range.");
 }

 bool use_maf_effect_s_prior_scale = maf_effect_s_prior_scale.isNotNull();
 Rcpp::NumericVector maf_effect_s_prior_scale_vec;
 if (use_maf_effect_s_prior_scale) {
  maf_effect_s_prior_scale_vec = Rcpp::NumericVector(maf_effect_s_prior_scale);
  use_maf_effect_s_prior_scale = maf_effect_s_prior_scale_vec.size() > 0;
 }

 if (use_maf_effect_s_prior_scale &&
     static_cast<int>(maf_effect_s_prior_scale_vec.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_prior_scale must have length m.");
 }
 if (estimate_maf_effect_s && use_maf_effect_s_prior_scale) {
  throw std::runtime_error(
   "stblr_cpg_omp_csr_bayesr: fixed maf_effect_s_prior_scale and estimate_maf_effect_s cannot both be used."
  );
 }
 if (estimate_maf_effect_s) {
  if (!std::isfinite(maf_effect_s_init)) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_init must be finite.");
  }
  if (maf_effect_s_prior.size() != 2) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_prior must have length 2.");
  }
  if (!std::isfinite(maf_effect_s_prior[0]) || !std::isfinite(maf_effect_s_prior[1]) ||
      maf_effect_s_prior[0] >= maf_effect_s_prior[1]) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_prior must have finite lower < upper.");
  }
  if (maf_effect_s_init < maf_effect_s_prior[0] || maf_effect_s_init > maf_effect_s_prior[1]) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_init must lie within maf_effect_s_prior.");
  }
  if (!std::isfinite(maf_effect_s_proposal_sd) || maf_effect_s_proposal_sd <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_proposal_sd must be positive.");
  }
  if (maf_effect_s_log_h.isNull()) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_log_h is required when estimate_maf_effect_s = TRUE.");
  }
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
 arma::rowvec maf_effect_s_log_h_row;

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

 if (use_maf_effect_s_prior_scale) {
  prior_scale.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double scale_i = maf_effect_s_prior_scale_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(scale_i) || scale_i <= 0.0) {
    throw std::runtime_error(
     "stblr_cpg_omp_csr_bayesr: maf_effect_s_prior_scale must contain positive finite values."
    );
   }
   prior_scale(static_cast<arma::uword>(i)) = scale_i;
  }
 }
 if (estimate_maf_effect_s) {
  Rcpp::NumericVector log_h_vec(maf_effect_s_log_h);
  if (static_cast<int>(log_h_vec.size()) != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_log_h must have length m.");
  }
  maf_effect_s_log_h_row.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
   const double log_h_i = log_h_vec[static_cast<std::size_t>(i)];
   if (!std::isfinite(log_h_i)) {
    throw std::runtime_error("stblr_cpg_omp_csr_bayesr: maf_effect_s_log_h must contain finite values.");
   }
   maf_effect_s_log_h_row(static_cast<arma::uword>(i)) = log_h_i;
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

 arma::rowvec xx_row(static_cast<arma::uword>(m));
 for (int i = 0; i < m; ++i) {
  xx_row(static_cast<arma::uword>(i)) = xx[static_cast<std::size_t>(i)];
 }
 auto operator_context = make_operator(
  m, xx, xx_row, wy_mat, updateLDswap, ld_swap_r2, ld_swap_max_friends
 );
 const auto& op = operator_context.op;
 const BayesRLDLDFriends& ld_swap_friends = operator_context.ld_swap_friends;

 CsrBayesRExecutionContext<typename std::decay<decltype(op)>::type> execution_context;
 execution_context.op=&const_cast<typename std::decay<decltype(op)>::type&>(op);
 execution_context.ld_swap_friends=&ld_swap_friends; execution_context.mixture_var=&mixture_var;
 execution_context.pi=&pi; execution_context.alpha=&alpha; execution_context.n=&n;
 execution_context.chain_seeds=&chain_seeds; execution_context.comp_init=&comp_init; execution_context.r_init=&r_init;
 execution_context.wy_mat=&wy_mat; execution_context.b_init_mat=&b_init_mat; execution_context.B=&B; execution_context.E=&E;
 execution_context.ssb_prior_mat=&ssb_prior_mat; execution_context.sse_prior_mat=&sse_prior_mat; execution_context.yy_vec=&yy_vec;
 execution_context.prior_scale=&prior_scale; execution_context.maf_effect_s_log_h_row=&maf_effect_s_log_h_row;
 execution_context.convergence_markers=&convergence_markers;
 execution_context.K=K; execution_context.m=m; execution_context.nt=nt; execution_context.nchains=nchains;
 execution_context.ncores=ncores; execution_context.nit=nit; execution_context.nburn=nburn; execution_context.nthin=nthin;
 execution_context.seed=seed; execution_context.updateE_start=updateE_start; execution_context.updateE_every=updateE_every;
 execution_context.ld_swap_moves=ld_swap_moves; execution_context.use_comp_init=use_comp_init; execution_context.use_r_init=use_r_init;
 execution_context.estimate_maf_effect_s=estimate_maf_effect_s; execution_context.use_maf_effect_s_prior_scale=use_maf_effect_s_prior_scale;
 execution_context.updateLDswap=updateLDswap; execution_context.updateB=updateB; execution_context.updateE=updateE; execution_context.updatePi=updatePi;
 execution_context.convergence_probability=convergence_probability;
 execution_context.convergence_b=convergence_b; execution_context.convergence_d=convergence_d;
 execution_context.convergence_component=convergence_component;
 execution_context.ld_swap_prob=ld_swap_prob; execution_context.maf_effect_s_init=maf_effect_s_init;
 execution_context.maf_effect_s_prior_lower=maf_effect_s_prior[0]; execution_context.maf_effect_s_prior_upper=maf_effect_s_prior[1];
 execution_context.maf_effect_s_proposal_sd=maf_effect_s_proposal_sd; execution_context.nub=nub; execution_context.nue=nue; execution_context.adjE=adjE;
 CsrBayesRExecutionResult execution_result=run_csr_bayesr(execution_context);
 auto stblr_csr_bayesr_result_to_raw = [&](const CsrBayesRExecutionResult& execution_result) -> Rcpp::List {
 const arma::vec& mixture_var_vec=execution_result.mixture_var_vec;
 const arma::mat &bm_task=execution_result.bm_task,&dm_task=execution_result.dm_task,&component_mean_task=execution_result.component_mean_task;
 const arma::mat &component_task=execution_result.component_task,&vbs_task=execution_result.vbs_task,&vgs_task=execution_result.vgs_task;
 const arma::mat &ves_task=execution_result.ves_task,&vles_task=execution_result.vles_task,&vlds_task=execution_result.vlds_task;
 const arma::mat &pis_task=execution_result.pis_task,&maf_effect_s_task=execution_result.maf_effect_s_task;
 const arma::mat &final_pi_task=execution_result.final_pi_task,&mean_pi_task=execution_result.mean_pi_task;
 const arma::vec &maf_effect_s_attempted_task=execution_result.maf_effect_s_attempted_task,&maf_effect_s_accepted_task=execution_result.maf_effect_s_accepted_task;
 const std::vector<arma::mat>& comp_prob_task=execution_result.comp_prob_task; const std::vector<arma::vec>& ncomp_task=execution_result.ncomp_task;
 const std::vector<arma::mat>& convergence_pi_task=execution_result.convergence_pi_task;
 const std::vector<arma::mat>& convergence_b_task=execution_result.convergence_b_task;
 const std::vector<arma::imat>& convergence_d_task=execution_result.convergence_d_task;
 const std::vector<arma::imat>& convergence_component_task=execution_result.convergence_component_task;
 const arma::mat &bm=execution_result.bm,&dm=execution_result.dm,&bm_sd=execution_result.bm_sd,&dm_sd=execution_result.dm_sd;
 const arma::mat &bm_min=execution_result.bm_min,&dm_min=execution_result.dm_min,&bm_max=execution_result.bm_max,&dm_max=execution_result.dm_max;
 const arma::mat &component_mean=execution_result.component_mean,&b_out=execution_result.b_out,&r_out=execution_result.r_out,&component_out=execution_result.component_out;
 const arma::mat &vbs=execution_result.vbs,&vgs=execution_result.vgs,&ves=execution_result.ves,&vle=execution_result.vle,&vld=execution_result.vld,&pis=execution_result.pis;
 const arma::mat &maf_effect_s=execution_result.maf_effect_s,&final_pi=execution_result.final_pi,&mean_pi=execution_result.mean_pi;
 const arma::mat &updateE_diagnostics=execution_result.updateE_diagnostics,&ld_swap_diagnostics=execution_result.ld_swap_diagnostics;
 const arma::mat& ld_swap_chain_diagnostics=execution_result.ld_swap_chain_diagnostics;
 const arma::vec &final_vb=execution_result.final_vb,&final_vg=execution_result.final_vg,&final_ve=execution_result.final_ve;
 const arma::vec &maf_effect_s_attempted=execution_result.maf_effect_s_attempted,&maf_effect_s_accepted=execution_result.maf_effect_s_accepted,&nsamples=execution_result.nsamples;
 const std::vector<arma::mat>& comp_prob=execution_result.comp_prob; const arma::mat& ncomp=execution_result.ncomp;
 const arma::mat &covb=execution_result.covb,&covg=execution_result.covg,&cove=execution_result.cove,&vb=execution_result.vb,&vg=execution_result.vg,&ve=execution_result.ve;

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

 Rcpp::CharacterVector component_names(K);
 for (int k = 0; k < K; ++k) {
  component_names[k] = "component_" + std::to_string(k);
 }

 Rcpp::List comp_prob_out(nt);
 for (int t = 0; t < nt; ++t) {
  comp_prob_out[t] = comp_prob[static_cast<std::size_t>(t)];
 }

 Rcpp::NumericVector selection_mean(nt);
 Rcpp::NumericVector maf_effect_sd(nt);
 Rcpp::NumericVector selection_min(nt);
 Rcpp::NumericVector selection_max(nt);
 Rcpp::NumericVector selection_acceptance(nt);
 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  selection_acceptance[t] =
   maf_effect_s_attempted(tu) > 0.0
   ? maf_effect_s_accepted(tu) / maf_effect_s_attempted(tu)
   : 0.0;

  if (estimate_maf_effect_s) {
   double mean_s = 0.0;
   double min_s = std::numeric_limits<double>::infinity();
   double max_s = -std::numeric_limits<double>::infinity();
   int ns = 0;
   for (int it = nburn; it < n_trace; ++it) {
    const double val = maf_effect_s(tu, static_cast<arma::uword>(it));
    mean_s += val;
    min_s = std::min(min_s, val);
    max_s = std::max(max_s, val);
    ++ns;
   }
   if (ns > 0) mean_s /= static_cast<double>(ns);
   double ss = 0.0;
   if (ns > 1) {
    for (int it = nburn; it < n_trace; ++it) {
     const double diff = maf_effect_s(tu, static_cast<arma::uword>(it)) - mean_s;
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

 Rcpp::IntegerVector n_used(nt);
 Rcpp::NumericVector nsamples_vec(nt);
 for (int t = 0; t < nt; ++t) {
  n_used[t] = n[static_cast<std::size_t>(t)];
  nsamples_vec[t] = nsamples(static_cast<arma::uword>(t));
 }

 Rcpp::List marker = Rcpp::List::create(
  Rcpp::Named("bm") = marker_matrix(bm),
  Rcpp::Named("dm") = marker_matrix(dm),
  Rcpp::Named("wy") = marker_matrix(wy_mat),
  Rcpp::Named("r") = marker_matrix(r_out),
  Rcpp::Named("b") = marker_matrix(b_out),
  Rcpp::Named("state") = marker_matrix(component_out)
 );
 if (return_chain_summaries) {
  marker["bm_sd"] = marker_matrix(bm_sd);
  marker["bm_min"] = marker_matrix(bm_min);
  marker["bm_max"] = marker_matrix(bm_max);
  marker["dm_sd"] = marker_matrix(dm_sd);
  marker["dm_min"] = marker_matrix(dm_min);
  marker["dm_max"] = marker_matrix(dm_max);
 }

 Rcpp::List trace = Rcpp::List::create(
  Rcpp::Named("vbs") = trace_matrix(vbs),
  Rcpp::Named("vgs") = trace_matrix(vgs),
  Rcpp::Named("ves") = trace_matrix(ves),
  Rcpp::Named("vle") = trace_matrix(vle),
  Rcpp::Named("vld") = trace_matrix(vld),
  Rcpp::Named("pis") = trace_matrix(pis)
 );

 Rcpp::List variance = Rcpp::List::create(
  Rcpp::Named("covb") = covb,
  Rcpp::Named("covg") = covg,
  Rcpp::Named("cove") = cove,
  Rcpp::Named("vb") = vb,
  Rcpp::Named("vg") = vg,
  Rcpp::Named("ve") = ve
 );

 Rcpp::List diagnostics = Rcpp::List::create(
  Rcpp::Named("nsamples") = nsamples_vec,
  Rcpp::Named("n_used") = n_used,
  Rcpp::Named("log_cpo") = Rcpp::NumericVector(nt),
  Rcpp::Named("mean_log_cpo") = Rcpp::NumericVector(nt),
  Rcpp::Named("seconds_mean") = Rcpp::NumericVector(nt),
  Rcpp::Named("seconds_max") = Rcpp::NumericVector(nt),
  Rcpp::Named("ld_swap") = updateLDswap ? Rcpp::wrap(ld_swap_diagnostics) : R_NilValue,
  Rcpp::Named("updateE") = updateE_diagnostics
 );
 if (operator_context.diagnostics.size() > 0) {
  diagnostics["block_eigen"] = operator_context.diagnostics;
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
    Rcpp::NumericVector chain_component_mean(m);
    for (int i = 0; i < m; ++i) {
     const arma::uword iu = static_cast<arma::uword>(i);
     chain_bm[i] = bm_task(task_u, iu);
     chain_dm[i] = dm_task(task_u, iu);
     chain_state[i] = component_task(task_u, iu);
     chain_component_mean[i] = component_mean_task(task_u, iu);
    }
    Rcpp::NumericVector updateE_diag(9);
    for (int j = 0; j < 9; ++j) {
     updateE_diag[j] = updateE_diagnostics(task_u, static_cast<arma::uword>(j));
    }
    Rcpp::NumericMatrix chain_ld(1, 3);
    chain_ld(0, 0) = ld_swap_chain_diagnostics(task_u, 2);
    chain_ld(0, 1) = ld_swap_chain_diagnostics(task_u, 3);
    chain_ld(0, 2) = ld_swap_chain_diagnostics(task_u, 4);

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

    Rcpp::List convergence_trace = Rcpp::List::create(
     Rcpp::Named("component_pi") = convergence_pi_task[static_cast<std::size_t>(task)],
     Rcpp::Named("b") = convergence_b_task[static_cast<std::size_t>(task)],
     Rcpp::Named("d") = convergence_d_task[static_cast<std::size_t>(task)],
     Rcpp::Named("component") = convergence_component_task[static_cast<std::size_t>(task)],
     Rcpp::Named("marker_index") = Rcpp::wrap(convergence_markers)
    );
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
      Rcpp::Named("final") = final_pi_task.row(task_u).t(),
      Rcpp::Named("mean") = mean_pi_task.row(task_u).t()
     ),
     Rcpp::Named("component") = Rcpp::List::create(
      Rcpp::Named("prob") = comp_prob_task[static_cast<std::size_t>(task)],
      Rcpp::Named("ncomp") = ncomp_task[static_cast<std::size_t>(task)].t(),
      Rcpp::Named("dm_component_mean") = chain_component_mean
     ),
     Rcpp::Named("selection") = chain_selection,
     Rcpp::Named("convergence_trace") = convergence_trace,
     Rcpp::Named("diagnostics") = Rcpp::List::create(
      Rcpp::Named("ld_swap") = updateLDswap ? Rcpp::wrap(chain_ld) : R_NilValue,
      Rcpp::Named("updateE") = updateE_diag
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
  Rcpp::Named("trace") = estimate_maf_effect_s ? Rcpp::wrap(trace_matrix(maf_effect_s)) : R_NilValue,
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
   Rcpp::Named("model") = "bayesr",
   Rcpp::Named("backend") = "csr_bayesr",
   Rcpp::Named("data_level") = "summary",
   Rcpp::Named("prior_type") = "component",
   Rcpp::Named("m") = m,
   Rcpp::Named("nt") = nt,
   Rcpp::Named("n_trace") = n_trace,
   Rcpp::Named("nit") = nit,
   Rcpp::Named("nburn") = nburn,
   Rcpp::Named("nthin") = nthin,
   Rcpp::Named("nchains") = nchains,
   Rcpp::Named("keep_chains") = keep_chains,
   Rcpp::Named("n_components") = K,
   Rcpp::Named("n_annotations") = 0,
   Rcpp::Named("n_groups") = 0
  ),
  Rcpp::Named("marker") = marker,
  Rcpp::Named("trace") = trace,
  Rcpp::Named("variance") = variance,
  Rcpp::Named("pi") = Rcpp::List::create(
   Rcpp::Named("final") = final_pi,
   Rcpp::Named("mean") = mean_pi,
   Rcpp::Named("names") = component_names
  ),
  Rcpp::Named("diagnostics") = diagnostics,
  Rcpp::Named("chains") = chains,
  Rcpp::Named("prior") = Rcpp::List::create(),
  Rcpp::Named("group") = Rcpp::List::create(),
  Rcpp::Named("annotation") = Rcpp::List::create(),
  Rcpp::Named("component") = Rcpp::List::create(
   Rcpp::Named("names") = component_names,
   Rcpp::Named("mixture_var") = mixture_var_vec,
   Rcpp::Named("prob") = comp_prob_out,
   Rcpp::Named("ncomp") = ncomp,
   Rcpp::Named("dm_component_mean") = marker_matrix(component_mean)
  ),
  Rcpp::Named("selection") = selection
 );
 raw.attr("class") = Rcpp::CharacterVector::create("stblr_raw_v1", "stblr_raw", "list");
 return raw;
 };
 return stblr_csr_bayesr_result_to_raw(execution_result);
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
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale = R_NilValue,
  bool estimate_maf_effect_s = false,
  double maf_effect_s_init = 0.0,
  Rcpp::NumericVector maf_effect_s_prior = Rcpp::NumericVector::create(-3.0, 2.0),
  double maf_effect_s_proposal_sd = 0.35,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h = R_NilValue,
  Rcpp::IntegerVector convergence_markers = Rcpp::IntegerVector::create(),
  bool convergence_probability = false,
  bool convergence_b = false,
  bool convergence_d = false,
  bool convergence_component = false
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
  BayesRLDLDFriends friends;
  if (update_ld_swap) {
   friends = build_ld_swap_friends_bayesr_ST_csr(
    m, ld, xx, ld_swap_r2_value, ld_swap_max_friends_value
   );
  } else {
   friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  }
  return BayesROperatorContext<CsrOperator>(
   std::move(op), std::move(friends), Rcpp::List::create()
  );
 };

 return stblr_cpg_omp_csr_bayesr_impl(
  wy, ww, yy, b_init, comp_init, use_comp_init, r_init, use_r_init,
  rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, pi,
  mixture_var, alpha, nub, nue, updateB, updateE, updatePi, adjE, n,
  nit, nburn, nthin, ncores, seed, nchains, keep_chains, chain_seeds,
  updateE_start, updateE_every, updateLDswap, ld_swap_prob, ld_swap_r2,
  ld_swap_max_friends, ld_swap_moves, maf_effect_s_prior_scale,
  estimate_maf_effect_s, maf_effect_s_init, maf_effect_s_prior,
  maf_effect_s_proposal_sd, maf_effect_s_log_h,
  Rcpp::as<std::vector<int>>(convergence_markers),
  convergence_probability, convergence_b, convergence_d,
  convergence_component, make_csr_operator
 );
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_bayesr_block_eigen(
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
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale = R_NilValue,
  bool estimate_maf_effect_s = false,
  double maf_effect_s_init = 0.0,
  Rcpp::NumericVector maf_effect_s_prior = Rcpp::NumericVector::create(-3.0, 2.0),
  double maf_effect_s_proposal_sd = 0.35,
  Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h = R_NilValue,
  Rcpp::IntegerVector convergence_markers = Rcpp::IntegerVector::create(),
  bool convergence_probability = false,
  bool convergence_b = false,
  bool convergence_d = false,
  bool convergence_component = false,
  Rcpp::CharacterVector bed_files = Rcpp::CharacterVector::create(),
  int n_bed = 0,
  Rcpp::List cls = R_NilValue,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::NumericVector af = Rcpp::NumericVector::create(),
  Rcpp::IntegerVector block_start = Rcpp::IntegerVector::create(),
  std::string eigen_filter = "hard_truncate",
  double eigen_tau = 0.01,
  double eigen_eta = 0.0
) {
 if (updateLDswap) {
  throw std::runtime_error(
   "LD-swap is not yet supported with the experimental block-eigen operator."
  );
 }
 if (bed_files.size() <= 0) {
  throw std::runtime_error("bed_files must contain at least one BED file.");
 }
 if (n_bed <= 0) throw std::runtime_error("n_bed must be positive.");
 if (cls.size() != bed_files.size()) {
  throw std::runtime_error("cls must have one element per BED file.");
 }

 const std::vector<std::string> bed_files_cpp =
  stblr_bayesr_copy_character_vector(bed_files, "bed_files");
 const std::vector<std::vector<int>> cls_cpp =
  stblr_bayesr_copy_int_list(cls, "cls");
 const std::vector<int> rows0 = stblr_bayesr_copy_rows0_or_empty(rows, n_bed);
 const std::vector<double> af_cpp = Rcpp::as<std::vector<double>>(af);
 const std::vector<int> block_start_cpp =
  Rcpp::as<std::vector<int>>(block_start);
 const EigenFilterMode mode = parse_block_eigen_filter_mode(eigen_filter);

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
   throw std::runtime_error("af length must equal m for block-eigen BayesR.");
  }
  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
   bed_files_cpp, n_bed, rows0.empty() ? nullptr : rows0.data(),
   static_cast<int>(rows0.size()), cls_cpp
  );
  if (G.m != m) throw std::runtime_error("BED marker count does not match m.");

  std::vector<BlockEigenDiag> block_diag;
  BlockEigenOperator op = build_block_eigen(
   G, af_cpp, block_start_cpp, mode, eigen_tau, eigen_eta,
   wy_mat, ncores, &block_diag
  );
  BayesRLDLDFriends friends;
  friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  Rcpp::List diagnostics = Rcpp::List::create(
   Rcpp::Named("blocks") = block_eigen_diagnostics_to_data_frame(block_diag)
  );
  return BayesROperatorContext<BlockEigenOperator>(
   std::move(op), std::move(friends), diagnostics
  );
 };

 return stblr_cpg_omp_csr_bayesr_impl(
  wy, ww, yy, b_init, comp_init, use_comp_init, r_init, use_r_init,
  rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, pi,
  mixture_var, alpha, nub, nue, updateB, updateE, updatePi, adjE, n,
  nit, nburn, nthin, ncores, seed, nchains, keep_chains, chain_seeds,
  updateE_start, updateE_every, updateLDswap, ld_swap_prob, ld_swap_r2,
  ld_swap_max_friends, ld_swap_moves, maf_effect_s_prior_scale,
  estimate_maf_effect_s, maf_effect_s_init, maf_effect_s_prior,
  maf_effect_s_proposal_sd, maf_effect_s_log_h,
  Rcpp::as<std::vector<int>>(convergence_markers),
  convergence_probability, convergence_b, convergence_d,
  convergence_component, make_block_eigen_operator
 );
}
