// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "st_chain_utils.h"
#include "st_csr_common.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
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
// STBLR summary-stat CSR with learned annotation-informed priors
// =============================================================================
//
// Exported function:
//
//   stblr_cpg_omp_csr_annot(...)
//
// Model extension relative to stblr_cpg_omp_csr_prior():
//
//   logit(pi_it) = logit(pi_global_t) + center_i(A_i eta_pi_t)
//
//   vb_it = vb_t * exp(center_i(A_i eta_vb_t)), capped by
//           [vb_multiplier_min, vb_multiplier_max]
//
// The annotation effects eta_pi_t and eta_vb_t are learned inside the sampler
// by random-walk Metropolis-Hastings. This is a practical hierarchical
// annotation-informed BayesC implementation.
//
// Notes:
//   * A is an m x K dense annotation matrix. Rows are markers, columns are
//     annotations. Overlapping annotations are allowed.
//   * For very large K or extremely sparse A, replace dense A by a CSR
//     annotation object later. Dense A is usually fine for K in the tens/hundreds.
//   * Updated to align with the current CSR BLR samplers:
//       - explicit Beta(pi_prior_a, pi_prior_b) prior for global pi
//       - initial VG/VLE/VLD computed from the starting state
//       - VLE and VLD traces returned as slots 20 and 21
//   * Annotation outputs are preserved:
//       result[18][t] = eta_pi posterior mean, length K
//       result[19][t] = eta_vb posterior mean, length K
//       result[20][t] = VLE trace
//       result[21][t] = VLD trace = VG - VLE
//
// =============================================================================

// struct STLDCSR {
//  std::vector<uint64_t> ptr;
//  std::vector<int> idx;
//  std::vector<float> xij;
// };
//
// inline void read_exact_file(
//   const std::string& path,
//   void* data,
//   std::size_t nbytes
// ) {
//  FILE* fs = std::fopen(path.c_str(), "rb");
//
//  if (!fs) {
//   throw std::runtime_error("Could not open file: " + path);
//  }
//
//  const std::size_t got = std::fread(data, 1, nbytes, fs);
//  std::fclose(fs);
//
//  if (got != nbytes) {
//   throw std::runtime_error("Short read from file: " + path);
//  }
// }
//
// inline uint64_t parse_uint64_from_meta(
//   const std::string& value,
//   const std::string& key
// ) {
//  if (value.empty()) {
//   throw std::runtime_error("Empty metadata value for key: " + key);
//  }
//
//  char* endptr = nullptr;
//  const unsigned long long out = std::strtoull(value.c_str(), &endptr, 10);
//
//  if (endptr == value.c_str() || *endptr != '\0') {
//   throw std::runtime_error("Invalid unsigned integer metadata value for key: " + key);
//  }
//
//  return static_cast<uint64_t>(out);
// }
//
// inline STLDCSR read_and_build_st_ld_csr(
//   const std::string& prefix,
//   int m,
//   const std::vector<double>& xx
// ) {
//  const std::string row_file  = prefix + ".row_ptr.u64.bin";
//  const std::string col_file  = prefix + ".col_idx.u32.0based.bin";
//  const std::string val_file  = prefix + ".values.f32.bin";
//  const std::string meta_file = prefix + ".meta.txt";
//
//  if (m <= 0) {
//   throw std::runtime_error("read_and_build_st_ld_csr: m must be positive.");
//  }
//
//  if (static_cast<int>(xx.size()) != m) {
//   throw std::runtime_error("read_and_build_st_ld_csr: xx must have length m.");
//  }
//
//  std::ifstream meta(meta_file.c_str());
//  if (!meta.is_open()) {
//   throw std::runtime_error("Could not open metadata file: " + meta_file);
//  }
//
//  int m_meta = -1;
//  uint64_t nnz_u64 = 0;
//  bool have_nnz = false;
//
//  std::string line;
//  while (std::getline(meta, line)) {
//   const std::string key_m   = "n_variants=";
//   const std::string key_nnz = "nnz=";
//
//   if (line.rfind(key_m, 0) == 0) {
//    m_meta = std::stoi(line.substr(key_m.size()));
//   } else if (line.rfind(key_nnz, 0) == 0) {
//    nnz_u64 = parse_uint64_from_meta(line.substr(key_nnz.size()), "nnz");
//    have_nnz = true;
//   }
//  }
//  meta.close();
//
//  if (m_meta <= 0) {
//   throw std::runtime_error("Could not read n_variants from metadata.");
//  }
//
//  if (m_meta != m) {
//   throw std::runtime_error("LD metadata n_variants does not match marker dimension.");
//  }
//
//  if (!have_nnz) {
//   throw std::runtime_error("Could not read nnz from metadata.");
//  }
//
//  const std::size_t nnz = static_cast<std::size_t>(nnz_u64);
//
//  std::vector<uint64_t> row_ptr(static_cast<std::size_t>(m) + 1);
//  std::vector<uint32_t> col_idx_u32(nnz);
//  std::vector<float> values_r(nnz);
//
//  read_exact_file(row_file, row_ptr.data(), row_ptr.size() * sizeof(uint64_t));
//  read_exact_file(col_file, col_idx_u32.data(), col_idx_u32.size() * sizeof(uint32_t));
//  read_exact_file(val_file, values_r.data(), values_r.size() * sizeof(float));
//
//  if (row_ptr[0] != 0 || row_ptr[static_cast<std::size_t>(m)] != nnz_u64) {
//   throw std::runtime_error("Invalid LD row_ptr: expected 0-based row_ptr ending at nnz.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   if (row_ptr[static_cast<std::size_t>(i + 1)] < row_ptr[static_cast<std::size_t>(i)]) {
//    throw std::runtime_error("Invalid LD row_ptr: row pointers are not nondecreasing.");
//   }
//
//   if (!std::isfinite(xx[static_cast<std::size_t>(i)]) || xx[static_cast<std::size_t>(i)] <= 0.0) {
//    throw std::runtime_error(
//      "read_and_build_st_ld_csr: xx contains invalid value at marker " +
//       std::to_string(i)
//    );
//   }
//  }
//
//  std::vector<uint64_t> degree(static_cast<std::size_t>(m), 0);
//
//  for (int i = 0; i < m; ++i) {
//   const uint64_t start = row_ptr[static_cast<std::size_t>(i)];
//   const uint64_t end   = row_ptr[static_cast<std::size_t>(i + 1)];
//
//   if (end > nnz_u64) {
//    throw std::runtime_error("Invalid LD row_ptr: row end exceeds nnz.");
//   }
//
//   for (uint64_t p = start; p < end; ++p) {
//    const uint32_t j_u32 = col_idx_u32[static_cast<std::size_t>(p)];
//
//    if (j_u32 >= static_cast<uint32_t>(m)) {
//     throw std::runtime_error("LD column index out of range.");
//    }
//
//    const int j = static_cast<int>(j_u32);
//
//    if (j == i) continue;
//
//    ++degree[static_cast<std::size_t>(i)];
//    ++degree[static_cast<std::size_t>(j)];
//   }
//  }
//
//  STLDCSR ld;
//  ld.ptr.resize(static_cast<std::size_t>(m) + 1);
//  ld.ptr[0] = 0;
//
//  for (int i = 0; i < m; ++i) {
//   ld.ptr[static_cast<std::size_t>(i + 1)] =
//    ld.ptr[static_cast<std::size_t>(i)] + degree[static_cast<std::size_t>(i)];
//  }
//
//  const uint64_t nnz_sym = ld.ptr[static_cast<std::size_t>(m)];
//
//  ld.idx.resize(static_cast<std::size_t>(nnz_sym));
//  ld.xij.resize(static_cast<std::size_t>(nnz_sym));
//
//  std::vector<uint64_t> offset = ld.ptr;
//
//  double max_abs_rij = 0.0;
//  double max_abs_xij = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const uint64_t start = row_ptr[static_cast<std::size_t>(i)];
//   const uint64_t end   = row_ptr[static_cast<std::size_t>(i + 1)];
//
//   for (uint64_t p = start; p < end; ++p) {
//    const int j = static_cast<int>(col_idx_u32[static_cast<std::size_t>(p)]);
//
//    if (j == i) continue;
//
//    const double rij = static_cast<double>(values_r[static_cast<std::size_t>(p)]);
//
//    if (!std::isfinite(rij)) {
//     throw std::runtime_error("LD value contains NaN/Inf.");
//    }
//
//    max_abs_rij = std::max(max_abs_rij, std::abs(rij));
//
//    if (std::abs(rij) > 1.0001) {
//     throw std::runtime_error(
//       "LD value is not a correlation. Did you pass X_i'X_j instead of r_ij?"
//     );
//    }
//
//    const double xij =
//     rij * std::sqrt(xx[static_cast<std::size_t>(i)] * xx[static_cast<std::size_t>(j)]);
//
//    if (!std::isfinite(xij)) {
//     throw std::runtime_error("Computed X_i'X_j contains NaN/Inf.");
//    }
//
//    max_abs_xij = std::max(max_abs_xij, std::abs(xij));
//
//    const float xij_f = static_cast<float>(xij);
//
//    const uint64_t pos_i = offset[static_cast<std::size_t>(i)]++;
//    ld.idx[static_cast<std::size_t>(pos_i)] = j;
//    ld.xij[static_cast<std::size_t>(pos_i)] = xij_f;
//
//    const uint64_t pos_j = offset[static_cast<std::size_t>(j)]++;
//    ld.idx[static_cast<std::size_t>(pos_j)] = i;
//    ld.xij[static_cast<std::size_t>(pos_j)] = xij_f;
//   }
//  }
//
//  for (int i = 0; i < m; ++i) {
//   if (offset[static_cast<std::size_t>(i)] != ld.ptr[static_cast<std::size_t>(i + 1)]) {
//    throw std::runtime_error("Internal LD CSR fill-count mismatch.");
//   }
//  }
//
//  Rcpp::Rcout
//  << "ST annotation flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
//  << ", symmetric nnz=" << static_cast<double>(nnz_sym)
//  << ", max_abs_rij=" << max_abs_rij
//  << ", max_abs_xij=" << max_abs_xij
//  << "\n";
//
//  return ld;
// }
//
// inline void rebuild_residual_st_csr(
//   int m,
//   const arma::rowvec& wy,
//   const arma::rowvec& ww,
//   const arma::rowvec& b,
//   arma::rowvec& r,
//   const STLDCSR& ld
// ) {
//  r = wy;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//   const double bi = b(iu);
//
//   if (bi == 0.0) continue;
//
//   r(iu) -= ww(iu) * bi;
//
//   const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
//   const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];
//
//   for (uint64_t p = start; p < end; ++p) {
//    const int j = ld.idx[static_cast<std::size_t>(p)];
//    r(static_cast<arma::uword>(j)) -=
//     static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * bi;
//   }
//  }
// }
//
// inline double clamp_prob(double x) {
//  if (!std::isfinite(x)) {
//   throw std::runtime_error("clamp_prob: probability is NaN/Inf.");
//  }
//  return std::min(std::max(x, 1e-300), 1.0 - 1e-12);
// }
//
// inline double inv_logit_stable(double x) {
//  if (x >= 0.0) {
//   const double z = std::exp(-x);
//   return 1.0 / (1.0 + z);
//  }
//  const double z = std::exp(x);
//  return z / (1.0 + z);
// }

// inline double log1pexp_stable(double x) {
//  if (x > 40.0) return x;
//  if (x < -40.0) return std::exp(x);
//  return std::log1p(std::exp(x));
// }
//
// inline double logit_prob(double p) {
//  p = clamp_prob(p);
//  return std::log(p) - std::log1p(-p);
// }
//
// inline double clamp_double(double x, double lo, double hi) {
//  return std::min(std::max(x, lo), hi);
// }

inline arma::rowvec centered_annotation_linear_predictor(
  const arma::mat& A,
  const arma::vec& eta
) {
 arma::rowvec z = (A * eta).t();
 const double mu = arma::mean(z);
 z -= mu;
 return z;
}

inline arma::rowvec make_pi_from_annotation(
  const arma::mat& A,
  const arma::vec& eta_pi,
  double base_pi,
  double pi_min,
  double pi_max
) {
 arma::rowvec z = centered_annotation_linear_predictor(A, eta_pi);
 const double base = logit_prob(base_pi);
 arma::rowvec out(z.n_elem, arma::fill::zeros);

 for (arma::uword i = 0; i < z.n_elem; ++i) {
  const double p = inv_logit_stable(base + z(i));
  out(i) = clamp_double(p, pi_min, pi_max);
 }

 return out;
}

inline arma::rowvec make_vb_multiplier_from_annotation(
  const arma::mat& A,
  const arma::vec& eta_vb,
  double mult_min,
  double mult_max
) {
 arma::rowvec z = centered_annotation_linear_predictor(A, eta_vb);
 arma::rowvec out(z.n_elem, arma::fill::ones);

 for (arma::uword i = 0; i < z.n_elem; ++i) {
  const double mult = std::exp(clamp_double(z(i), std::log(mult_min), std::log(mult_max)));
  out(i) = clamp_double(mult, mult_min, mult_max);
 }

 return out;
}

inline double logpost_eta_pi(
  const arma::mat& A,
  const arma::vec& eta_pi,
  const arma::Row<int>& d,
  double base_pi,
  double pi_min,
  double pi_max,
  double sigma_eta
) {
 arma::rowvec z = centered_annotation_linear_predictor(A, eta_pi);
 const double base = logit_prob(base_pi);
 double lp = 0.0;

 for (arma::uword i = 0; i < z.n_elem; ++i) {
  const double eta_i = base + z(i);

  // If caps are inactive for most markers, this is the exact Bernoulli-logit
  // likelihood. When caps bind, it is still a stable pseudo-posterior update.
  double p = inv_logit_stable(eta_i);
  p = clamp_double(p, pi_min, pi_max);

  if (d(i) > 0) lp += std::log(std::max(p, 1e-300));
  else          lp += std::log(std::max(1.0 - p, 1e-300));
 }

 if (!std::isfinite(sigma_eta) || sigma_eta <= 0.0) {
  throw std::runtime_error("logpost_eta_pi: sigma_eta must be positive.");
 }

 lp += -0.5 * arma::dot(eta_pi, eta_pi) / (sigma_eta * sigma_eta);
 return lp;
}

inline double logpost_eta_vb(
  const arma::mat& A,
  const arma::vec& eta_vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  double vb_global,
  double mult_min,
  double mult_max,
  double sigma_eta
) {
 if (!std::isfinite(vb_global) || vb_global <= 0.0) {
  throw std::runtime_error("logpost_eta_vb: vb_global must be positive.");
 }

 arma::rowvec z = centered_annotation_linear_predictor(A, eta_vb);
 double lp = 0.0;

 for (arma::uword i = 0; i < z.n_elem; ++i) {
  if (d(i) <= 0) continue;

  const double log_mult = clamp_double(z(i), std::log(mult_min), std::log(mult_max));
  const double mult = std::exp(log_mult);
  const double vbi = std::max(vb_global * mult, 1e-300);
  const double bi = b(i);

  // Active-effect normal likelihood ignoring constants:
  // -0.5 * log(vbi) - 0.5 * b_i^2 / vbi
  lp += -0.5 * std::log(vbi) - 0.5 * bi * bi / vbi;
 }

 if (!std::isfinite(sigma_eta) || sigma_eta <= 0.0) {
  throw std::runtime_error("logpost_eta_vb: sigma_eta must be positive.");
 }

 lp += -0.5 * arma::dot(eta_vb, eta_vb) / (sigma_eta * sigma_eta);
 return lp;
}

inline void rw_update_eta_pi(
  const arma::mat& A,
  arma::vec& eta_pi,
  const arma::Row<int>& d,
  double base_pi,
  double pi_min,
  double pi_max,
  double sigma_eta,
  double rw_sd,
  std::mt19937& gen,
  int& accepted,
  int& proposed
) {
 if (eta_pi.n_elem == 0 || rw_sd <= 0.0) return;

 std::normal_distribution<double> norm01(0.0, 1.0);
 std::uniform_real_distribution<double> runif(0.0, 1.0);

 arma::vec prop = eta_pi;
 for (arma::uword k = 0; k < prop.n_elem; ++k) {
  prop(k) += rw_sd * norm01(gen);
 }

 const double lp_old = logpost_eta_pi(A, eta_pi, d, base_pi, pi_min, pi_max, sigma_eta);
 const double lp_new = logpost_eta_pi(A, prop,   d, base_pi, pi_min, pi_max, sigma_eta);

 ++proposed;

 const double logr = lp_new - lp_old;
 if (std::isfinite(logr) && (logr >= 0.0 || std::log(runif(gen)) < logr)) {
  eta_pi = prop;
  ++accepted;
 }
}

inline void rw_update_eta_vb(
  const arma::mat& A,
  arma::vec& eta_vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  double vb_global,
  double mult_min,
  double mult_max,
  double sigma_eta,
  double rw_sd,
  std::mt19937& gen,
  int& accepted,
  int& proposed
) {
 if (eta_vb.n_elem == 0 || rw_sd <= 0.0) return;

 std::normal_distribution<double> norm01(0.0, 1.0);
 std::uniform_real_distribution<double> runif(0.0, 1.0);

 arma::vec prop = eta_vb;
 for (arma::uword k = 0; k < prop.n_elem; ++k) {
  prop(k) += rw_sd * norm01(gen);
 }

 const double lp_old = logpost_eta_vb(A, eta_vb, b, d, vb_global, mult_min, mult_max, sigma_eta);
 const double lp_new = logpost_eta_vb(A, prop,   b, d, vb_global, mult_min, mult_max, sigma_eta);

 ++proposed;

 const double logr = lp_new - lp_old;
 if (std::isfinite(logr) && (logr >= 0.0 || std::log(runif(gen)) < logr)) {
  eta_vb = prop;
  ++accepted;
 }
}

inline void sampleBetaC_ST_csr_prior(
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
  throw std::runtime_error("sampleBetaC_ST_csr_prior: invalid ww value.");
 }

 if (!std::isfinite(vb_t) || vb_t <= 0.0) {
  throw std::runtime_error("sampleBetaC_ST_csr_prior: invalid global vb.");
 }

 if (!std::isfinite(vb_mult_i) || vb_mult_i <= 0.0) {
  throw std::runtime_error("sampleBetaC_ST_csr_prior: invalid vb multiplier.");
 }

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm01(0.0, 1.0);

 pi1_i = clamp_prob(pi1_i);
 const double pi0_i = std::max(1.0 - pi1_i, 1e-300);
 const double vbi = std::max(vb_t * vb_mult_i, 1e-12);
 const double vei_safe = std::max(vei_i, 1e-300);

 const double score = r(iu) + wi * b(iu);
 const double denom = std::max(vei_safe + wi * vbi, 1e-300);

 const double logBF =
  0.5 * std::log(vei_safe / denom)
  + 0.5 * score * score * vbi / (vei_safe * denom);

 const double logp1 = std::log(pi1_i) + logBF;
 const double logp0 = std::log(pi0_i);
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

inline void sampleB_ST_csr_prior(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  const arma::rowvec& vb_multiplier,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb_scaled = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);

  if (d(iu) > 0) {
   const double mult = vb_multiplier(iu);

   if (!std::isfinite(mult) || mult <= 0.0) {
    throw std::runtime_error("sampleB_ST_csr_prior: invalid vb multiplier.");
   }

   ssb_scaled += b(iu) * b(iu) / mult;
   dfb += 1.0;
  }
 }

 const double scale = ssb_scaled + nub * ssb_prior;

 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleB_ST_csr_prior: invalid scale.");
 }

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);

 vb = std::max(scale / chi2, 1e-12);
}

// inline void sampleE_ST_csr(
//   int m,
//   double nue,
//   double& ve,
//   const arma::rowvec& b,
//   const arma::rowvec& wy,
//   const arma::rowvec& r,
//   double sse_prior,
//   double yy,
//   int n,
//   std::mt19937& gen
// ) {
//  double b_dot_r_plus_wy = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//   b_dot_r_plus_wy += b(iu) * (r(iu) + wy(iu));
//  }
//
//  const double sse = yy - b_dot_r_plus_wy;
//  const double scale = sse + nue * sse_prior;
//
//  if (!std::isfinite(scale) || scale <= 0.0) {
//   throw std::runtime_error("sampleE_ST_csr: invalid residual scale.");
//  }
//
//  std::chi_squared_distribution<double> rchisq(n + nue);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  const double ve_new = scale / chi2;
//
//  if (!std::isfinite(ve_new) || ve_new <= 0.0) {
//   throw std::runtime_error("sampleE_ST_csr: sampled ve is invalid.");
//  }
//
//  ve = std::max(ve_new, 1e-12);
// }
//
// inline double computeG_ST_csr(
//   const arma::rowvec& b,
//   const arma::rowvec& wy,
//   const arma::rowvec& r,
//   int n
// ) {
//  double ssg = 0.0;
//  const arma::uword m = b.n_elem;
//
//  for (arma::uword i = 0; i < m; ++i) {
//   ssg += b(i) * (wy(i) - r(i));
//  }
//
//  return ssg / static_cast<double>(n);
// }

inline double computeLE_ST_csr_annot(
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

struct LDLDFriendsAnnot {
 std::vector<uint64_t> ptr;
 std::vector<int> idx;
 std::vector<double> r2;
};

inline LDLDFriendsAnnot build_ld_swap_friends_st_csr_annot(
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

 LDLDFriendsAnnot friends;
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

inline void set_marker_effect_st_csr_annot(
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

inline double residual_sse_st_csr_annot(
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

inline int count_excluded_ld_friends_annot(
  int i,
  const arma::Row<int>& d,
  const LDLDFriendsAnnot& friends
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

inline int collect_ld_swap_candidates_annot(
  int m,
  const arma::Row<int>& d,
  const LDLDFriendsAnnot& friends,
  std::vector<int>& candidates,
  std::vector<int>& n_excluded
) {
 candidates.clear();
 n_excluded.clear();

 for (int i = 0; i < m; ++i) {
  if (d(static_cast<arma::uword>(i)) <= 0) continue;

  const int nf = count_excluded_ld_friends_annot(i, d, friends);
  if (nf > 0) {
   candidates.push_back(i);
   n_excluded.push_back(nf);
  }
 }

 return static_cast<int>(candidates.size());
}

inline double log_marker_prior_ratio_annot(
  double b_move,
  int j,
  int k,
  double vb_t,
  const arma::rowvec& pi_marker,
  bool use_pi_marker,
  const arma::rowvec& vb_multiplier,
  bool use_vb_multiplier
) {
 double log_ratio = 0.0;

 if (use_pi_marker) {
  const double pi_j = clamp_prob(pi_marker(static_cast<arma::uword>(j)));
  const double pi_k = clamp_prob(pi_marker(static_cast<arma::uword>(k)));
  log_ratio +=
   std::log(pi_k) + std::log(std::max(1.0 - pi_j, 1e-300)) -
   std::log(pi_j) - std::log(std::max(1.0 - pi_k, 1e-300));
 }

 if (use_vb_multiplier) {
  const double mult_j = vb_multiplier(static_cast<arma::uword>(j));
  const double mult_k = vb_multiplier(static_cast<arma::uword>(k));
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
 }

 return log_ratio;
}

inline bool attempt_ld_swap_st_csr_annot(
  int m,
  double vei,
  double yy,
  double vb_t,
  const arma::rowvec& ww,
  const arma::rowvec& wy,
  const arma::rowvec& pi_marker,
  bool use_pi_marker,
  const arma::rowvec& vb_multiplier,
  bool use_vb_multiplier,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const STLDCSR& ld,
  const LDLDFriendsAnnot& friends,
  std::mt19937& gen
) {
 if (!std::isfinite(vei) || vei <= 0.0) return false;

 std::vector<int> candidates;
 std::vector<int> n_excluded;
 const int n_candidates = collect_ld_swap_candidates_annot(m, d, friends, candidates, n_excluded);
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

 const double sse_old = residual_sse_st_csr_annot(m, b, wy, r, yy);
 if (!std::isfinite(sse_old)) return false;

 const arma::rowvec r_old = r;

 set_marker_effect_st_csr_annot(j, 0.0, 0, ww, r, b, d, ld);
 set_marker_effect_st_csr_annot(k, b_j_old, 1, ww, r, b, d, ld);

 const double sse_new = residual_sse_st_csr_annot(m, b, wy, r, yy);
 bool accept = false;

 if (std::isfinite(sse_new)) {
  std::vector<int> reverse_candidates;
  std::vector<int> reverse_n_excluded;
  const int n_reverse_candidates =
   collect_ld_swap_candidates_annot(m, d, friends, reverse_candidates, reverse_n_excluded);

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
   const double log_prior_ratio = log_marker_prior_ratio_annot(
    b_j_old, j, k, vb_t, pi_marker, use_pi_marker,
    vb_multiplier, use_vb_multiplier
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

inline void samplePi_ST_annot(
  const arma::Row<int>& d,
  std::vector<double>& pi,
  double pi_prior_a,
  double pi_prior_b,
  std::mt19937& gen
) {
 // pi[1] is the inclusion probability; pi[0] is the exclusion probability.
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

#define SBLR_CSR_LEARNED_ANNOTATION_BAYESC_CORE_IMPL_TRANSLATION_UNIT 1
#include "blr_csr_learned_annotation_bayesc_core_impl.h"
#undef SBLR_CSR_LEARNED_ANNOTATION_BAYESC_CORE_IMPL_TRANSLATION_UNIT

sblr::core::CsrLearnedAnnotationBayesCExecutionResult stblr_cpg_omp_csr_annot_single(
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
  arma::mat A,
  bool learn_pi_annot,
  bool learn_vb_annot,
  arma::mat eta_pi_init,
  arma::mat eta_vb_init,
  double sigma_eta_pi,
  double sigma_eta_vb,
  double rw_sd_eta_pi,
  double rw_sd_eta_vb,
  int annot_update_every,
  double pi_min,
  double pi_max,
  double vb_multiplier_min,
  double vb_multiplier_max,
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
  int ld_swap_moves = 1,
  const std::vector<int>& convergence_markers=std::vector<int>(),
  bool convergence_annotations=false,
  bool convergence_b=false,
  bool convergence_d=false
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: nt must be positive.");
 }

 const int m = static_cast<int>(wy[0].size());

 if (m <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: m must be positive.");
 }

 if (nit <= 0 || nburn < 0 || nthin <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: invalid nit/nburn/nthin.");
 }

 if (annot_update_every <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: annot_update_every must be positive.");
 }

 if (!std::isfinite(pi_prior_a) || pi_prior_a <= 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: pi_prior_a must be finite and positive.");
 }

 if (!std::isfinite(pi_prior_b) || pi_prior_b <= 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: pi_prior_b must be finite and positive.");
 }

 if ((int)ww.size() != nt || (int)b_init.size() != nt ||
     (int)yy.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: inconsistent trait dimensions.");
 }

 if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: priors must be nt x nt.");
 }

 if (pi.size() != 2) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: pi must have length 2, c(pi0, pi1).");
 }

 if ((int)A.n_rows != m) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: A must have m rows.");
 }

 const int K = static_cast<int>(A.n_cols);

 if (K <= 0 && (learn_pi_annot || learn_vb_annot)) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: A must have at least one column when learning annotations.");
 }

 if (learn_pi_annot) {
  if ((int)eta_pi_init.n_rows != K || (int)eta_pi_init.n_cols != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_annot: eta_pi_init must be K x nt.");
  }
 }

 if (learn_vb_annot) {
  if ((int)eta_vb_init.n_rows != K || (int)eta_vb_init.n_cols != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_annot: eta_vb_init must be K x nt.");
  }
 }

 if (!std::isfinite(pi_min) || !std::isfinite(pi_max) || pi_min <= 0.0 || pi_max >= 1.0 || pi_min >= pi_max) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: invalid pi_min/pi_max.");
 }

 if (!std::isfinite(vb_multiplier_min) || !std::isfinite(vb_multiplier_max) ||
     vb_multiplier_min <= 0.0 || vb_multiplier_min >= vb_multiplier_max) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: invalid vb multiplier bounds.");
 }

 if (!std::isfinite(ld_swap_prob) || ld_swap_prob < 0.0 || ld_swap_prob > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: ld_swap_prob must be in [0, 1].");
 }

 if (!std::isfinite(ld_swap_r2) || ld_swap_r2 < 0.0 || ld_swap_r2 > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: ld_swap_r2 must be in [0, 1].");
 }

 if (ld_swap_max_friends <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: ld_swap_max_friends must be positive.");
 }

 if (ld_swap_moves < 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: ld_swap_moves must be non-negative.");
 }

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m ||
      (int)ww[t].size() != m ||
      (int)b_init[t].size() != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_annot: inconsistent marker dimensions.");
  }
 }

 if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: B must be nt x nt.");
 }

 if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: E must be nt x nt.");
 }

 if (use_r_init) {
  if (static_cast<int>(r_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_annot: r_init must have length nt when use_r_init = true.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(r_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_annot: each r_init[t] must have length m.");
   }
  }
 }

 if (use_d_init) {
  if (static_cast<int>(d_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_annot: d_init must have length nt when use_d_init = true.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(d_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_annot: each d_init[t] must have length m.");
   }
  }
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
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword iu = static_cast<arma::uword>(i);

   wy_mat(tu, iu) = wy[t][i];
   ww_mat(tu, iu) = ww[t][i];
   b_mat(tu, iu)  = b_init[t][i];
  }

  if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_annot: priors must be nt x nt.");
  }

  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
   sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
  }
 }

 for (int t = 1; t < nt; ++t) {
  if (n[t] != n[0]) {
   throw std::runtime_error(
     "stblr_cpg_omp_csr_annot: current shared-LD scaling assumes equal n across traits."
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
      "stblr_cpg_omp_csr_annot: ww contains invalid value before LD pre-scaling."
    );
   }

   if (std::abs(w0 - wt) > tol) {
    throw std::runtime_error(
      "stblr_cpg_omp_csr_annot: ww differs across traits; pre-scaled shared ST LD is invalid."
    );
   }
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);

 for (int i = 0; i < m; ++i) {
  const double wi = ww_mat(0, static_cast<arma::uword>(i));
  if (!std::isfinite(wi) || wi <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_annot: ww contains invalid value in trait 0.");
  }
  xx[static_cast<std::size_t>(i)] = wi;
 }

 STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);

 LDLDFriendsAnnot ld_swap_friends;
 if (updateLDswap) {
  ld_swap_friends = build_ld_swap_friends_st_csr_annot(
   m,
   ld,
   xx,
   ld_swap_r2,
   ld_swap_max_friends
  );
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


 sblr::core::LearnedAnnotationBayesCPolicyView annotation_policy;
 annotation_policy.marker_count=static_cast<std::size_t>(m);
 annotation_policy.annotation_count=static_cast<std::size_t>(K);
 annotation_policy.trait_count=static_cast<std::size_t>(nt);
 annotation_policy.annotation=A.memptr();
 annotation_policy.annotation_value_count=A.n_elem;
 annotation_policy.annotation_order.reserve(static_cast<std::size_t>(K));
 for (int k=0; k<K; ++k)
  annotation_policy.annotation_order.push_back("annotation_"+std::to_string(k));
 annotation_policy.layout="column_major_marker_by_annotation";
 annotation_policy.includes_intercept=false;
 annotation_policy.eta_probability_init=eta_pi_init.memptr();
 annotation_policy.eta_probability_count=eta_pi_init.n_elem;
 annotation_policy.eta_multiplier_init=eta_vb_init.memptr();
 annotation_policy.eta_multiplier_count=eta_vb_init.n_elem;
 annotation_policy.learn_probability=learn_pi_annot;
 annotation_policy.learn_multiplier=learn_vb_annot;
 annotation_policy.probability_prior_sd=sigma_eta_pi;
 annotation_policy.multiplier_prior_sd=sigma_eta_vb;
 annotation_policy.probability_proposal_sd=rw_sd_eta_pi;
 annotation_policy.multiplier_proposal_sd=rw_sd_eta_vb;
 annotation_policy.update_every=annot_update_every;
 annotation_policy.probability_min=pi_min;
 annotation_policy.probability_max=pi_max;
 annotation_policy.multiplier_min=vb_multiplier_min;
 annotation_policy.multiplier_max=vb_multiplier_max;

 sblr::core::CsrLearnedAnnotationBayesCExecutionContext context;
 context.marker_count=m; context.trait_count=nt; context.annotation_count=K;
 context.iterations=nit; context.burnin=nburn; context.thinning=nthin;
 context.cores=ncores; context.seed=seed;
 context.use_initial_inclusion=use_d_init;
 context.use_initial_residual=use_r_init;
 context.rebuild_residual_before_update=rebuild_r_before_updateE;
 context.learn_probability=learn_pi_annot;
 context.learn_multiplier=learn_vb_annot;
 context.update_marker_variance=updateB;
 context.update_residual_variance=updateE;
 context.update_global_probability=updatePi;
 context.update_ld_swap=updateLDswap;
 context.annotation_update_every=annot_update_every;
 context.ld_swap_moves=ld_swap_moves;
 context.probability_prior_sd=sigma_eta_pi;
 context.multiplier_prior_sd=sigma_eta_vb;
 context.probability_proposal_sd=rw_sd_eta_pi;
 context.multiplier_proposal_sd=rw_sd_eta_vb;
 context.probability_min=pi_min; context.probability_max=pi_max;
 context.multiplier_min=vb_multiplier_min;
 context.multiplier_max=vb_multiplier_max;
 context.marker_degrees_freedom=nub; context.residual_degrees_freedom=nue;
 context.residual_adjustment=adjE;
 context.global_probability_prior_a=pi_prior_a;
 context.global_probability_prior_b=pi_prior_b;
 context.ld_swap_probability=ld_swap_prob;
 context.marker_score=&wy_mat; context.marker_diagonal=&ww_mat;
 context.effect=&b_mat; context.residual=&r_mat; context.inclusion=&d_mat;
 context.phenotype_sum_squares=&yy_vec;
 context.marker_scale_prior=&ssb_prior_mat;
 context.residual_scale_prior=&sse_prior_mat;
 context.marker_variance_initial=&B; context.residual_variance_initial=&E;
 context.global_probability=&pi; context.sample_size=&n;
 context.initial_inclusion=&d_init; context.initial_residual=&r_init;
 context.annotation=&A;
 context.initial_probability_coefficient=&eta_pi_init;
 context.initial_multiplier_coefficient=&eta_vb_init;
 context.marker_order=&order; context.ld_storage=&ld;
 context.ld_row_ptr_count=ld.ptr.size();
 context.ld_friends_storage=&ld_swap_friends;
 context.annotation_policy=annotation_policy;
 context.convergence_markers=&convergence_markers;
 context.convergence_annotations=convergence_annotations;
 context.convergence_b=convergence_b;
 context.convergence_d=convergence_d;

 sblr::core::CsrLearnedAnnotationBayesCExecutionResult execution_result=
  sblr::core::run_csr_learned_annotation_bayesc(context);
 return execution_result;
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

static Rcpp::NumericMatrix cpg_raw_annotation_matrix(
 const std::vector<std::vector<std::vector<double>>>& raw,
 std::size_t slot,
 int nanno,
 int nt
) {
 Rcpp::NumericMatrix out(nanno, nt);
 if (slot >= raw.size()) return out;
 for (int t = 0; t < nt; ++t) {
  const std::vector<double>& x = raw[slot][static_cast<std::size_t>(t)];
  for (int a = 0; a < nanno && static_cast<std::size_t>(a) < x.size(); ++a) out(a, t) = x[a];
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

static Rcpp::List cpg_annot_chains_raw_v1(
 const std::vector<std::vector<std::vector<double>>>& raw,
 int m,
 int nt,
 int nanno,
 int nchains,
 const std::vector<arma::mat>& convergence_eta_pi,
 const std::vector<arma::mat>& convergence_eta_vb,
 const std::vector<arma::mat>& convergence_b,
 const std::vector<arma::imat>& convergence_d,
 const std::vector<int>& convergence_markers
) {
 if (raw.size() <= 34) return Rcpp::List::create();
 Rcpp::List traits(nt);
 for (int t = 0; t < nt; ++t) {
  Rcpp::List chains(nchains);
  for (int cc = 0; cc < nchains; ++cc) {
   Rcpp::NumericVector dm(m), bm(m), eta_pi(nanno), eta_vb(nanno);
   for (int i = 0; i < m; ++i) {
    const std::size_t idx = static_cast<std::size_t>(cc * m + i);
    if (idx < raw[30][t].size()) dm[i] = raw[30][t][idx];
    if (idx < raw[31][t].size()) bm[i] = raw[31][t][idx];
   }
   for (int a = 0; a < nanno; ++a) {
    const std::size_t idx = static_cast<std::size_t>(cc * nanno + a);
    if (idx < raw[33][t].size()) eta_pi[a] = raw[33][t][idx];
    if (idx < raw[34][t].size()) eta_vb[a] = raw[34][t][idx];
   }
   Rcpp::NumericMatrix ld(1, 3);
   const std::size_t didx = static_cast<std::size_t>(cc * 4);
   if (didx + 2 < raw[32][t].size()) {
    ld(0, 0) = raw[32][t][didx];
    ld(0, 1) = raw[32][t][didx + 1];
    ld(0, 2) = raw[32][t][didx + 2];
   }
   chains[cc] = Rcpp::List::create(
    Rcpp::_["marker"] = Rcpp::List::create(Rcpp::_["bm"] = bm, Rcpp::_["dm"] = dm),
    Rcpp::_["annotation"] = Rcpp::List::create(
     Rcpp::_["eta_pi_mean"] = eta_pi,
     Rcpp::_["eta_vb_mean"] = eta_vb
    ),
    Rcpp::_["convergence_trace"] = Rcpp::List::create(
     Rcpp::_["marker_index"] = Rcpp::wrap(convergence_markers),
     Rcpp::_["inclusion_coefficient"] = convergence_eta_pi[static_cast<std::size_t>(cc*nt+t)],
     Rcpp::_["variance_coefficient"] = convergence_eta_vb[static_cast<std::size_t>(cc*nt+t)],
     Rcpp::_["b"] = convergence_b[static_cast<std::size_t>(cc*nt+t)],
     Rcpp::_["d"] = convergence_d[static_cast<std::size_t>(cc*nt+t)],
     Rcpp::_["component"] = R_NilValue),
    Rcpp::_["diagnostics"] = Rcpp::List::create(Rcpp::_["ld_swap"] = ld)
   );
  }
  traits[t] = chains;
 }
 return traits;
}

static Rcpp::List stblr_csr_learned_annotation_bayesc_result_to_raw(
 const sblr::core::CsrLearnedAnnotationBayesCExecutionResult& result,
 bool updateLDswap,
 int m,
 int nt,
 int nanno,
 int nit,
 int nburn,
 int nthin,
 int nchains,
 bool keep_chains,
 const std::vector<arma::mat>& convergence_eta_pi,
 const std::vector<arma::mat>& convergence_eta_vb,
 const std::vector<arma::mat>& convergence_b,
 const std::vector<arma::imat>& convergence_d,
 const std::vector<int>& convergence_markers
) {
 const std::vector<std::vector<std::vector<double>>>& raw=result.raw;
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
 if (updateLDswap && raw.size() > 22) {
  Rcpp::NumericMatrix ld(nt, 3);
  for (int t = 0; t < nt; ++t) {
   for (int j = 0; j < 3 && static_cast<std::size_t>(j) < raw[22][t].size(); ++j) ld(t, j) = raw[22][t][j];
  }
  ld_swap = ld;
 }
 Rcpp::RObject chains = keep_chains ? Rcpp::RObject(cpg_annot_chains_raw_v1(
  raw,m,nt,nanno,nchains,convergence_eta_pi,convergence_eta_vb,
  convergence_b,convergence_d,convergence_markers)) : Rcpp::RObject(R_NilValue);
 Rcpp::List raw_out = Rcpp::List::create(
  Rcpp::_["schema"] = Rcpp::List::create(Rcpp::_["class"] = "stblr_raw", Rcpp::_["version"] = 1),
  Rcpp::_["meta"] = Rcpp::List::create(
   Rcpp::_["model"] = "bayesc",
   Rcpp::_["backend"] = "csr_annot_bayesc",
   Rcpp::_["data_level"] = "summary",
   Rcpp::_["prior_type"] = "annotation",
   Rcpp::_["m"] = m,
   Rcpp::_["nt"] = nt,
   Rcpp::_["n_trace"] = n_trace,
   Rcpp::_["nit"] = nit,
   Rcpp::_["nburn"] = nburn,
   Rcpp::_["nthin"] = nthin,
   Rcpp::_["nchains"] = nchains,
   Rcpp::_["keep_chains"] = keep_chains,
   Rcpp::_["n_components"] = 2,
   Rcpp::_["n_annotations"] = nanno,
   Rcpp::_["n_groups"] = 0
  ),
  Rcpp::_["marker"] = Rcpp::List::create(
   Rcpp::_["bm"] = cpg_raw_marker_matrix(raw, 0, m, nt),
   Rcpp::_["dm"] = cpg_raw_marker_matrix(raw, 1, m, nt),
   Rcpp::_["wy"] = cpg_raw_marker_matrix(raw, 2, m, nt),
   Rcpp::_["r"] = cpg_raw_marker_matrix(raw, 3, m, nt),
   Rcpp::_["b"] = cpg_raw_marker_matrix(raw, 4, m, nt),
   Rcpp::_["state"] = cpg_raw_marker_matrix(raw, 5, m, nt),
   Rcpp::_["bm_sd"] = raw.size() > 24 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 24, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["bm_min"] = raw.size() > 25 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 25, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["bm_max"] = raw.size() > 26 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 26, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["dm_sd"] = raw.size() > 27 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 27, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["dm_min"] = raw.size() > 28 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 28, m, nt)) : Rcpp::RObject(R_NilValue),
   Rcpp::_["dm_max"] = raw.size() > 29 ? Rcpp::RObject(cpg_raw_marker_matrix(raw, 29, m, nt)) : Rcpp::RObject(R_NilValue)
  ),
  Rcpp::_["trace"] = Rcpp::List::create(
   Rcpp::_["vbs"] = cpg_raw_trace_matrix(raw, 7, n_trace, nt),
   Rcpp::_["vgs"] = cpg_raw_trace_matrix(raw, 8, n_trace, nt),
   Rcpp::_["ves"] = cpg_raw_trace_matrix(raw, 9, n_trace, nt),
   Rcpp::_["vle"] = cpg_raw_trace_matrix(raw, 20, n_trace, nt),
   Rcpp::_["vld"] = cpg_raw_trace_matrix(raw, 21, n_trace, nt),
   Rcpp::_["pis"] = cpg_raw_trace_matrix(raw, 23, n_trace, nt)
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
  Rcpp::_["group"] = Rcpp::List::create(),
  Rcpp::_["annotation"] = Rcpp::List::create(
   Rcpp::_["annotation_names"] = R_NilValue,
   Rcpp::_["eta_pi_mean"] = cpg_raw_annotation_matrix(raw, 18, nanno, nt),
   Rcpp::_["eta_pi_final"] = R_NilValue,
   Rcpp::_["eta_vb_mean"] = cpg_raw_annotation_matrix(raw, 19, nanno, nt),
   Rcpp::_["eta_vb_final"] = R_NilValue,
   Rcpp::_["marker_pi_mean"] = R_NilValue,
   Rcpp::_["marker_pi_final"] = R_NilValue,
   Rcpp::_["marker_vb_multiplier_mean"] = R_NilValue,
   Rcpp::_["marker_vb_multiplier_final"] = R_NilValue,
   Rcpp::_["eta_pi_acceptance"] = R_NilValue,
   Rcpp::_["eta_vb_acceptance"] = R_NilValue
  ),
  Rcpp::_["component"] = Rcpp::List::create(),
  Rcpp::_["selection"] = cpg_disabled_selection()
 );
 raw_out.attr("class") = Rcpp::CharacterVector::create("stblr_raw_v1", "stblr_raw", "list");
 return raw_out;
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_annot(
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
  arma::mat A,
  bool learn_pi_annot,
  bool learn_vb_annot,
  arma::mat eta_pi_init,
  arma::mat eta_vb_init,
  double sigma_eta_pi,
  double sigma_eta_vb,
  double rw_sd_eta_pi,
  double rw_sd_eta_vb,
  int annot_update_every,
  double pi_min,
  double pi_max,
  double vb_multiplier_min,
  double vb_multiplier_max,
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
  int nchains = 1,
  bool keep_chains = false,
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds = R_NilValue,
  bool updateLDswap = false,
  double ld_swap_prob = 0.05,
  double ld_swap_r2 = 0.8,
  int ld_swap_max_friends = 50,
  int ld_swap_moves = 1,
  Rcpp::IntegerVector convergence_markers=Rcpp::IntegerVector::create(),
  bool convergence_annotations=false,
  bool convergence_b=false,
  bool convergence_d=false
) {
 if (nchains <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: nchains must be positive.");
 }
 std::vector<int> chain_seeds_vec;
 if (chain_seeds.isNotNull()) {
  Rcpp::IntegerVector chain_seeds_r(chain_seeds);
  chain_seeds_vec = Rcpp::as<std::vector<int>>(chain_seeds_r);
 }
 if (!chain_seeds_vec.empty() && static_cast<int>(chain_seeds_vec.size()) != nchains) {
  throw std::runtime_error("stblr_cpg_omp_csr_annot: chain_seeds must have length nchains.");
 }

 std::vector<std::vector<std::vector<double>>> out;
 std::vector<std::vector<std::vector<double>>> sumsq;
 std::vector<std::vector<std::vector<double>>> minv;
 std::vector<std::vector<std::vector<double>>> maxv;
 std::vector<std::vector<double>> ld_swap_attempted_sum;
 std::vector<std::vector<double>> ld_swap_accepted_sum;
 std::vector<std::vector<double>> chain_dm_flat;
 std::vector<std::vector<double>> chain_bm_flat;
 std::vector<std::vector<double>> chain_ld_swap_flat;
 std::vector<std::vector<double>> chain_eta_pi_flat;
 std::vector<std::vector<double>> chain_eta_vb_flat;
 const std::vector<int> convergence_markers_cpp=
  Rcpp::as<std::vector<int>>(convergence_markers);
 std::vector<arma::mat> convergence_eta_pi,convergence_eta_vb,
  convergence_b_trace;
 std::vector<arma::imat> convergence_d_trace;
 const std::size_t task_count=static_cast<std::size_t>(nchains*wy.size());
 convergence_eta_pi.reserve(task_count); convergence_eta_vb.reserve(task_count);
 convergence_b_trace.reserve(task_count); convergence_d_trace.reserve(task_count);

 for (int chain = 0; chain < nchains; ++chain) {
  int chain_seed = seed;
  if (!chain_seeds_vec.empty()) {
   chain_seed = chain_seeds_vec[static_cast<std::size_t>(chain)];
  } else if (nchains > 1) {
   chain_seed = seed + 9176 * (chain + 1);
  }

  sblr::core::CsrLearnedAnnotationBayesCExecutionResult chain_result =
   stblr_cpg_omp_csr_annot_single(
    wy, ww, yy, b_init, d_init, use_d_init, r_init, use_r_init,
    rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, pi,
    A, learn_pi_annot, learn_vb_annot, eta_pi_init, eta_vb_init,
    sigma_eta_pi, sigma_eta_vb, rw_sd_eta_pi, rw_sd_eta_vb,
    annot_update_every, pi_min, pi_max, vb_multiplier_min,
    vb_multiplier_max, nub, nue, updateB, updateE, updatePi, adjE, n,
    nit, nburn, nthin, pi_prior_a, pi_prior_b, ncores, chain_seed,
    updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends, ld_swap_moves,
    convergence_markers_cpp,convergence_annotations,convergence_b,convergence_d
   );
  std::vector<std::vector<std::vector<double>>>& raw=chain_result.raw;
  for (std::size_t t=0;t<chain_result.convergence_eta_pi.size();++t) {
   convergence_eta_pi.push_back(chain_result.convergence_eta_pi[t]);
   convergence_eta_vb.push_back(chain_result.convergence_eta_vb[t]);
   convergence_b_trace.push_back(chain_result.convergence_b[t]);
   convergence_d_trace.push_back(chain_result.convergence_d[t]);
  }

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
   chain_ld_swap_flat.resize(raw[22].size());
   chain_eta_pi_flat.resize(raw[18].size());
   chain_eta_vb_flat.resize(raw[19].size());
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
   if (raw.size() > 22 && raw[22][t].size() >= 2) {
    ld_swap_attempted_sum[t][0] += raw[22][t][0];
    ld_swap_accepted_sum[t][0] += raw[22][t][1];
   }
  }

  if (keep_chains) {
   for (std::size_t t = 0; t < raw[0].size(); ++t) {
    chain_dm_flat[t].insert(chain_dm_flat[t].end(), raw[1][t].begin(), raw[1][t].end());
    chain_bm_flat[t].insert(chain_bm_flat[t].end(), raw[0][t].begin(), raw[0][t].end());
    chain_ld_swap_flat[t].insert(chain_ld_swap_flat[t].end(), raw[22][t].begin(), raw[22][t].end());
    chain_eta_pi_flat[t].insert(chain_eta_pi_flat[t].end(), raw[18][t].begin(), raw[18][t].end());
    chain_eta_vb_flat[t].insert(chain_eta_vb_flat[t].end(), raw[19][t].begin(), raw[19][t].end());
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

 if (out.size() > 22) {
  for (std::size_t t = 0; t < out[22].size(); ++t) {
   if (out[22][t].size() < 4) out[22][t].resize(4);
   const double attempted = ld_swap_attempted_sum[t][0];
   const double accepted = ld_swap_accepted_sum[t][0];
   out[22][t][0] = attempted;
   out[22][t][1] = accepted;
   out[22][t][2] = attempted > 0.0 ? accepted / attempted : 0.0;
   out[22][t][3] = 1.0;
  }
 }

 const bool return_chain_summaries = (nchains > 1) || keep_chains;
 if (return_chain_summaries) {
  std::vector<std::vector<std::vector<double>>> extended(keep_chains ? 35 : 30);
  for (std::size_t k = 0; k < out.size(); ++k) extended[k] = out[k];
  for (int slot = 24; slot <= 29; ++slot) extended[static_cast<std::size_t>(slot)].resize(out[0].size());
  for (std::size_t t = 0; t < out[0].size(); ++t) {
   const std::size_t m_t = out[0][t].size();
   for (int slot = 24; slot <= 29; ++slot) extended[static_cast<std::size_t>(slot)][t].resize(m_t);
   for (std::size_t i = 0; i < m_t; ++i) {
    const double bm_mean = out[0][t][i];
    const double dm_mean = out[1][t][i];
    double bm_sd = 0.0;
    double dm_sd = 0.0;
    if (nchains > 1) {
     bm_sd = std::sqrt(std::max(0.0, (sumsq[0][t][i] - nchains * bm_mean * bm_mean) / static_cast<double>(nchains - 1)));
     dm_sd = std::sqrt(std::max(0.0, (sumsq[1][t][i] - nchains * dm_mean * dm_mean) / static_cast<double>(nchains - 1)));
    }
    extended[24][t][i] = bm_sd;
    extended[25][t][i] = minv[0][t][i];
    extended[26][t][i] = maxv[0][t][i];
    extended[27][t][i] = dm_sd;
    extended[28][t][i] = minv[1][t][i];
    extended[29][t][i] = maxv[1][t][i];
   }
  }
  if (keep_chains) {
   extended[30] = chain_dm_flat;
   extended[31] = chain_bm_flat;
   extended[32] = chain_ld_swap_flat;
   extended[33] = chain_eta_pi_flat;
   extended[34] = chain_eta_vb_flat;
  }
  sblr::core::CsrLearnedAnnotationBayesCExecutionResult execution_result;
  execution_result.raw=std::move(extended);
  return stblr_csr_learned_annotation_bayesc_result_to_raw(
   execution_result, updateLDswap, static_cast<int>(wy[0].size()), static_cast<int>(wy.size()),
   static_cast<int>(A.n_cols), nit, nburn, nthin, nchains, keep_chains,
   convergence_eta_pi,convergence_eta_vb,convergence_b_trace,
   convergence_d_trace,convergence_markers_cpp
  );
 }

 sblr::core::CsrLearnedAnnotationBayesCExecutionResult execution_result;
 execution_result.raw=std::move(out);
 return stblr_csr_learned_annotation_bayesc_result_to_raw(
  execution_result, updateLDswap, static_cast<int>(wy[0].size()), static_cast<int>(wy.size()),
  static_cast<int>(A.n_cols), nit, nburn, nthin, nchains, keep_chains,
  convergence_eta_pi,convergence_eta_vb,convergence_b_trace,
  convergence_d_trace,convergence_markers_cpp
 );
}

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
// // STBLR summary-stat CSR with learned annotation-informed priors
// // =============================================================================
// //
// // Exported function:
// //
// //   stblr_cpg_omp_csr_annot(...)
// //
// // Model extension relative to stblr_cpg_omp_csr_prior():
// //
// //   logit(pi_it) = logit(pi_global_t) + center_i(A_i eta_pi_t)
// //
// //   vb_it = vb_t * exp(center_i(A_i eta_vb_t)), capped by
// //           [vb_multiplier_min, vb_multiplier_max]
// //
// // The annotation effects eta_pi_t and eta_vb_t are learned inside the sampler
// // by random-walk Metropolis-Hastings. This is a practical hierarchical
// // annotation-informed BayesC implementation.
// //
// // Notes:
// //   * A is an m x K dense annotation matrix. Rows are markers, columns are
// //     annotations. Overlapping annotations are allowed.
// //   * For very large K or extremely sparse A, replace dense A by a CSR
// //     annotation object later. Dense A is usually fine for K in the tens/hundreds.
// //   * The 20-slot return structure is preserved. Slots 18 and 19 are used for
// //     compact annotation diagnostics:
// //       result[18][t] = eta_pi posterior mean, length K
// //       result[19][t] = eta_vb posterior mean, length K
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
// //  << "ST annotation flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
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
// //
// // inline double clamp_prob(double x) {
// //  if (!std::isfinite(x)) {
// //   throw std::runtime_error("clamp_prob: probability is NaN/Inf.");
// //  }
// //  return std::min(std::max(x, 1e-300), 1.0 - 1e-12);
// // }
// //
// // inline double inv_logit_stable(double x) {
// //  if (x >= 0.0) {
// //   const double z = std::exp(-x);
// //   return 1.0 / (1.0 + z);
// //  }
// //  const double z = std::exp(x);
// //  return z / (1.0 + z);
// // }
//
// // inline double log1pexp_stable(double x) {
// //  if (x > 40.0) return x;
// //  if (x < -40.0) return std::exp(x);
// //  return std::log1p(std::exp(x));
// // }
// //
// // inline double logit_prob(double p) {
// //  p = clamp_prob(p);
// //  return std::log(p) - std::log1p(-p);
// // }
// //
// // inline double clamp_double(double x, double lo, double hi) {
// //  return std::min(std::max(x, lo), hi);
// // }
//
// inline arma::rowvec centered_annotation_linear_predictor(
//   const arma::mat& A,
//   const arma::vec& eta
// ) {
//  arma::rowvec z = (A * eta).t();
//  const double mu = arma::mean(z);
//  z -= mu;
//  return z;
// }
//
// inline arma::rowvec make_pi_from_annotation(
//   const arma::mat& A,
//   const arma::vec& eta_pi,
//   double base_pi,
//   double pi_min,
//   double pi_max
// ) {
//  arma::rowvec z = centered_annotation_linear_predictor(A, eta_pi);
//  const double base = logit_prob(base_pi);
//  arma::rowvec out(z.n_elem, arma::fill::zeros);
//
//  for (arma::uword i = 0; i < z.n_elem; ++i) {
//   const double p = inv_logit_stable(base + z(i));
//   out(i) = clamp_double(p, pi_min, pi_max);
//  }
//
//  return out;
// }
//
// inline arma::rowvec make_vb_multiplier_from_annotation(
//   const arma::mat& A,
//   const arma::vec& eta_vb,
//   double mult_min,
//   double mult_max
// ) {
//  arma::rowvec z = centered_annotation_linear_predictor(A, eta_vb);
//  arma::rowvec out(z.n_elem, arma::fill::ones);
//
//  for (arma::uword i = 0; i < z.n_elem; ++i) {
//   const double mult = std::exp(clamp_double(z(i), std::log(mult_min), std::log(mult_max)));
//   out(i) = clamp_double(mult, mult_min, mult_max);
//  }
//
//  return out;
// }
//
// inline double logpost_eta_pi(
//   const arma::mat& A,
//   const arma::vec& eta_pi,
//   const arma::Row<int>& d,
//   double base_pi,
//   double pi_min,
//   double pi_max,
//   double sigma_eta
// ) {
//  arma::rowvec z = centered_annotation_linear_predictor(A, eta_pi);
//  const double base = logit_prob(base_pi);
//  double lp = 0.0;
//
//  for (arma::uword i = 0; i < z.n_elem; ++i) {
//   const double eta_i = base + z(i);
//
//   // If caps are inactive for most markers, this is the exact Bernoulli-logit
//   // likelihood. When caps bind, it is still a stable pseudo-posterior update.
//   double p = inv_logit_stable(eta_i);
//   p = clamp_double(p, pi_min, pi_max);
//
//   if (d(i) > 0) lp += std::log(std::max(p, 1e-300));
//   else          lp += std::log(std::max(1.0 - p, 1e-300));
//  }
//
//  if (!std::isfinite(sigma_eta) || sigma_eta <= 0.0) {
//   throw std::runtime_error("logpost_eta_pi: sigma_eta must be positive.");
//  }
//
//  lp += -0.5 * arma::dot(eta_pi, eta_pi) / (sigma_eta * sigma_eta);
//  return lp;
// }
//
// inline double logpost_eta_vb(
//   const arma::mat& A,
//   const arma::vec& eta_vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   double vb_global,
//   double mult_min,
//   double mult_max,
//   double sigma_eta
// ) {
//  if (!std::isfinite(vb_global) || vb_global <= 0.0) {
//   throw std::runtime_error("logpost_eta_vb: vb_global must be positive.");
//  }
//
//  arma::rowvec z = centered_annotation_linear_predictor(A, eta_vb);
//  double lp = 0.0;
//
//  for (arma::uword i = 0; i < z.n_elem; ++i) {
//   if (d(i) <= 0) continue;
//
//   const double log_mult = clamp_double(z(i), std::log(mult_min), std::log(mult_max));
//   const double mult = std::exp(log_mult);
//   const double vbi = std::max(vb_global * mult, 1e-300);
//   const double bi = b(i);
//
//   // Active-effect normal likelihood ignoring constants:
//   // -0.5 * log(vbi) - 0.5 * b_i^2 / vbi
//   lp += -0.5 * std::log(vbi) - 0.5 * bi * bi / vbi;
//  }
//
//  if (!std::isfinite(sigma_eta) || sigma_eta <= 0.0) {
//   throw std::runtime_error("logpost_eta_vb: sigma_eta must be positive.");
//  }
//
//  lp += -0.5 * arma::dot(eta_vb, eta_vb) / (sigma_eta * sigma_eta);
//  return lp;
// }
//
// inline void rw_update_eta_pi(
//   const arma::mat& A,
//   arma::vec& eta_pi,
//   const arma::Row<int>& d,
//   double base_pi,
//   double pi_min,
//   double pi_max,
//   double sigma_eta,
//   double rw_sd,
//   std::mt19937& gen,
//   int& accepted,
//   int& proposed
// ) {
//  if (eta_pi.n_elem == 0 || rw_sd <= 0.0) return;
//
//  std::normal_distribution<double> norm01(0.0, 1.0);
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//
//  arma::vec prop = eta_pi;
//  for (arma::uword k = 0; k < prop.n_elem; ++k) {
//   prop(k) += rw_sd * norm01(gen);
//  }
//
//  const double lp_old = logpost_eta_pi(A, eta_pi, d, base_pi, pi_min, pi_max, sigma_eta);
//  const double lp_new = logpost_eta_pi(A, prop,   d, base_pi, pi_min, pi_max, sigma_eta);
//
//  ++proposed;
//
//  const double logr = lp_new - lp_old;
//  if (std::isfinite(logr) && (logr >= 0.0 || std::log(runif(gen)) < logr)) {
//   eta_pi = prop;
//   ++accepted;
//  }
// }
//
// inline void rw_update_eta_vb(
//   const arma::mat& A,
//   arma::vec& eta_vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   double vb_global,
//   double mult_min,
//   double mult_max,
//   double sigma_eta,
//   double rw_sd,
//   std::mt19937& gen,
//   int& accepted,
//   int& proposed
// ) {
//  if (eta_vb.n_elem == 0 || rw_sd <= 0.0) return;
//
//  std::normal_distribution<double> norm01(0.0, 1.0);
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//
//  arma::vec prop = eta_vb;
//  for (arma::uword k = 0; k < prop.n_elem; ++k) {
//   prop(k) += rw_sd * norm01(gen);
//  }
//
//  const double lp_old = logpost_eta_vb(A, eta_vb, b, d, vb_global, mult_min, mult_max, sigma_eta);
//  const double lp_new = logpost_eta_vb(A, prop,   b, d, vb_global, mult_min, mult_max, sigma_eta);
//
//  ++proposed;
//
//  const double logr = lp_new - lp_old;
//  if (std::isfinite(logr) && (logr >= 0.0 || std::log(runif(gen)) < logr)) {
//   eta_vb = prop;
//   ++accepted;
//  }
// }
//
// inline void sampleBetaC_ST_csr_prior(
//   int i,
//   double pi1_i,
//   double vb_t,
//   double vb_mult_i,
//   double vei_i,
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
//  if (!std::isfinite(wi) || wi <= 0.0) {
//   throw std::runtime_error("sampleBetaC_ST_csr_prior: invalid ww value.");
//  }
//
//  if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//   throw std::runtime_error("sampleBetaC_ST_csr_prior: invalid global vb.");
//  }
//
//  if (!std::isfinite(vb_mult_i) || vb_mult_i <= 0.0) {
//   throw std::runtime_error("sampleBetaC_ST_csr_prior: invalid vb multiplier.");
//  }
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm01(0.0, 1.0);
//
//  pi1_i = clamp_prob(pi1_i);
//  const double pi0_i = std::max(1.0 - pi1_i, 1e-300);
//  const double vbi = std::max(vb_t * vb_mult_i, 1e-12);
//  const double vei_safe = std::max(vei_i, 1e-300);
//
//  const double score = r(iu) + wi * b(iu);
//  const double denom = std::max(vei_safe + wi * vbi, 1e-300);
//
//  const double logBF =
//   0.5 * std::log(vei_safe / denom)
//   + 0.5 * score * score * vbi / (vei_safe * denom);
//
//  const double logp1 = std::log(pi1_i) + logBF;
//  const double logp0 = std::log(pi0_i);
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
//  double b_new = 0.0;
//
//  if (di == 1) {
//   const double lhs = wi + vei_safe / vbi;
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
//  d(iu) = di;
// }
//
// inline void sampleB_ST_csr_prior(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   const arma::rowvec& vb_multiplier,
//   double ssb_prior,
//   std::mt19937& gen
// ) {
//  double ssb_scaled = 0.0;
//  double dfb = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   const arma::uword iu = static_cast<arma::uword>(i);
//
//   if (d(iu) > 0) {
//    const double mult = vb_multiplier(iu);
//
//    if (!std::isfinite(mult) || mult <= 0.0) {
//     throw std::runtime_error("sampleB_ST_csr_prior: invalid vb multiplier.");
//    }
//
//    ssb_scaled += b(iu) * b(iu) / mult;
//    dfb += 1.0;
//   }
//  }
//
//  const double scale = ssb_scaled + nub * ssb_prior;
//
//  if (!std::isfinite(scale) || scale <= 0.0) {
//   throw std::runtime_error("sampleB_ST_csr_prior: invalid scale.");
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
// //
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
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_annot(
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
//   arma::mat A,
//   bool learn_pi_annot,
//   bool learn_vb_annot,
//   arma::mat eta_pi_init,
//   arma::mat eta_vb_init,
//   double sigma_eta_pi,
//   double sigma_eta_vb,
//   double rw_sd_eta_pi,
//   double rw_sd_eta_vb,
//   int annot_update_every,
//   double pi_min,
//   double pi_max,
//   double vb_multiplier_min,
//   double vb_multiplier_max,
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
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: nt must be positive.");
//  }
//
//  const int m = static_cast<int>(wy[0].size());
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: m must be positive.");
//  }
//
//  if (nit <= 0 || nburn < 0 || nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: invalid nit/nburn/nthin.");
//  }
//
//  if (annot_update_every <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: annot_update_every must be positive.");
//  }
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt ||
//      (int)yy.size() != nt || (int)n.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: inconsistent trait dimensions.");
//  }
//
//  if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: priors must be nt x nt.");
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: pi must have length 2, c(pi0, pi1).");
//  }
//
//  if ((int)A.n_rows != m) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: A must have m rows.");
//  }
//
//  const int K = static_cast<int>(A.n_cols);
//
//  if (K <= 0 && (learn_pi_annot || learn_vb_annot)) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: A must have at least one column when learning annotations.");
//  }
//
//  if (learn_pi_annot) {
//   if ((int)eta_pi_init.n_rows != K || (int)eta_pi_init.n_cols != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_annot: eta_pi_init must be K x nt.");
//   }
//  }
//
//  if (learn_vb_annot) {
//   if ((int)eta_vb_init.n_rows != K || (int)eta_vb_init.n_cols != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_annot: eta_vb_init must be K x nt.");
//   }
//  }
//
//  if (!std::isfinite(pi_min) || !std::isfinite(pi_max) || pi_min <= 0.0 || pi_max >= 1.0 || pi_min >= pi_max) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: invalid pi_min/pi_max.");
//  }
//
//  if (!std::isfinite(vb_multiplier_min) || !std::isfinite(vb_multiplier_max) ||
//      vb_multiplier_min <= 0.0 || vb_multiplier_min >= vb_multiplier_max) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: invalid vb multiplier bounds.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m ||
//       (int)ww[t].size() != m ||
//       (int)b_init[t].size() != m) {
//    throw std::runtime_error("stblr_cpg_omp_csr_annot: inconsistent marker dimensions.");
//   }
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_annot: E must be nt x nt.");
//  }
//
//  if (use_r_init) {
//   if (static_cast<int>(r_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_annot: r_init must have length nt when use_r_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(r_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_annot: each r_init[t] must have length m.");
//    }
//   }
//  }
//
//  if (use_d_init) {
//   if (static_cast<int>(d_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_annot: d_init must have length nt when use_d_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(d_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_annot: each d_init[t] must have length m.");
//    }
//   }
//  }
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
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword iu = static_cast<arma::uword>(i);
//
//    wy_mat(tu, iu) = wy[t][i];
//    ww_mat(tu, iu) = ww[t][i];
//    b_mat(tu, iu)  = b_init[t][i];
//   }
//
//   if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_annot: priors must be nt x nt.");
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
//      "stblr_cpg_omp_csr_annot: current shared-LD scaling assumes equal n across traits."
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
//       "stblr_cpg_omp_csr_annot: ww contains invalid value before LD pre-scaling."
//     );
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr_annot: ww differs across traits; pre-scaled shared ST LD is invalid."
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
//    throw std::runtime_error("stblr_cpg_omp_csr_annot: ww contains invalid value in trait 0.");
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
//  arma::vec final_pi(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  arma::mat eta_pi_mean(K, nt, arma::fill::zeros);
//  arma::mat eta_vb_mean(K, nt, arma::fill::zeros);
//  arma::mat eta_pi_final(K, nt, arma::fill::zeros);
//  arma::mat eta_vb_final(K, nt, arma::fill::zeros);
//  arma::mat eta_pi_accept(nt, 2, arma::fill::zeros);
//  arma::mat eta_vb_accept(nt, 2, arma::fill::zeros);
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
//  << "STBLR annotation CSR OpenMP requested threads = "
//  << nthreads
//  << ", omp_get_max_threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
//  Rcpp::Rcout
//  << "STBLR annotation CSR: K=" << K
//  << ", learn_pi_annot=" << learn_pi_annot
//  << ", learn_vb_annot=" << learn_vb_annot
//  << ", annot_update_every=" << annot_update_every
//  << ", pi bounds=[" << pi_min << "," << pi_max << "]"
//  << ", vb multiplier bounds=[" << vb_multiplier_min << "," << vb_multiplier_max << "]"
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
//      throw std::runtime_error("stblr_cpg_omp_csr_annot: r_init contains NaN/Inf.");
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
//    arma::vec eta_pi_t(K, arma::fill::zeros);
//    arma::vec eta_vb_t(K, arma::fill::zeros);
//
//    if (learn_pi_annot) {
//     eta_pi_t = eta_pi_init.col(static_cast<arma::uword>(t));
//    }
//
//    if (learn_vb_annot) {
//     eta_vb_t = eta_vb_init.col(static_cast<arma::uword>(t));
//    }
//
//    arma::rowvec pi_marker_t(m, arma::fill::zeros);
//    arma::rowvec vb_multiplier_t(m, arma::fill::ones);
//
//    if (learn_pi_annot) {
//     pi_marker_t = make_pi_from_annotation(A, eta_pi_t, pi_t[1], pi_min, pi_max);
//    } else {
//     pi_marker_t.fill(clamp_double(pi_t[1], pi_min, pi_max));
//    }
//
//    if (learn_vb_annot) {
//     vb_multiplier_t = make_vb_multiplier_from_annotation(
//      A,
//      eta_vb_t,
//      vb_multiplier_min,
//      vb_multiplier_max
//     );
//    }
//
//    arma::rowvec bm_t(m, arma::fill::zeros);
//    arma::rowvec dm_t(m, arma::fill::zeros);
//    arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
//    arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
//
//    arma::vec eta_pi_accum(K, arma::fill::zeros);
//    arma::vec eta_vb_accum(K, arma::fill::zeros);
//
//    int eta_pi_accepted = 0;
//    int eta_pi_proposed = 0;
//    int eta_vb_accepted = 0;
//    int eta_vb_proposed = 0;
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//     for (int isort = 0; isort < m; ++isort) {
//      const int i = order[static_cast<std::size_t>(isort)];
//      const arma::uword iu = static_cast<arma::uword>(i);
//
//      sampleBetaC_ST_csr_prior(
//       i,
//       pi_marker_t(iu),
//       vb_t,
//       vb_multiplier_t(iu),
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
//     if ((it + 1) % annot_update_every == 0) {
//      if (learn_pi_annot) {
//       rw_update_eta_pi(
//        A,
//        eta_pi_t,
//        d_t,
//        pi_t[1],
//            pi_min,
//            pi_max,
//            sigma_eta_pi,
//            rw_sd_eta_pi,
//            gen_t,
//            eta_pi_accepted,
//            eta_pi_proposed
//       );
//
//       pi_marker_t = make_pi_from_annotation(A, eta_pi_t, pi_t[1], pi_min, pi_max);
//      }
//
//      if (learn_vb_annot) {
//       rw_update_eta_vb(
//        A,
//        eta_vb_t,
//        b_t,
//        d_t,
//        vb_t,
//        vb_multiplier_min,
//        vb_multiplier_max,
//        sigma_eta_vb,
//        rw_sd_eta_vb,
//        gen_t,
//        eta_vb_accepted,
//        eta_vb_proposed
//       );
//
//       vb_multiplier_t = make_vb_multiplier_from_annotation(
//        A,
//        eta_vb_t,
//        vb_multiplier_min,
//        vb_multiplier_max
//       );
//      }
//     }
//
//     if (updateB) {
//      sampleB_ST_csr_prior(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       vb_multiplier_t,
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
//
//      if (learn_pi_annot) {
//       pi_marker_t = make_pi_from_annotation(A, eta_pi_t, pi_t[1], pi_min, pi_max);
//      } else {
//       pi_marker_t.fill(clamp_double(pi_t[1], pi_min, pi_max));
//      }
//     }
//
//     vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
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
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       const arma::uword iu = static_cast<arma::uword>(i);
//       bm_t(iu) += b_t(iu);
//       dm_t(iu) += static_cast<double>(d_t(iu));
//      }
//
//      if (K > 0) {
//       eta_pi_accum += eta_pi_t;
//       eta_vb_accum += eta_vb_t;
//      }
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    bm_t /= nsamples_t;
//    dm_t /= nsamples_t;
//
//    if (K > 0) {
//     eta_pi_accum /= nsamples_t;
//     eta_vb_accum /= nsamples_t;
//    }
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
//    if (K > 0) {
//     eta_pi_mean.col(static_cast<arma::uword>(t)) = eta_pi_accum;
//     eta_vb_mean.col(static_cast<arma::uword>(t)) = eta_vb_accum;
//     eta_pi_final.col(static_cast<arma::uword>(t)) = eta_pi_t;
//     eta_vb_final.col(static_cast<arma::uword>(t)) = eta_vb_t;
//    }
//
//    eta_pi_accept(static_cast<arma::uword>(t), 0) = eta_pi_accepted;
//    eta_pi_accept(static_cast<arma::uword>(t), 1) = eta_pi_proposed;
//    eta_vb_accept(static_cast<arma::uword>(t), 0) = eta_vb_accepted;
//    eta_vb_accept(static_cast<arma::uword>(t), 1) = eta_vb_proposed;
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
//      "stblr_cpg_omp_csr_annot failed for trait " +
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
//   result[18][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(std::max(K, 1)));
//   result[19][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(std::max(K, 1)));
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
//   if (K > 0) {
//    for (int k = 0; k < K; ++k) {
//     const std::size_t ks = static_cast<std::size_t>(k);
//     result[18][ts][ks] = eta_pi_mean(static_cast<arma::uword>(k), static_cast<arma::uword>(t));
//     result[19][ts][ks] = eta_vb_mean(static_cast<arma::uword>(k), static_cast<arma::uword>(t));
//    }
//   } else {
//    result[18][ts][0] = 0.0;
//    result[19][ts][0] = 0.0;
//   }
//  }
//
//  Rcpp::Rcout << "Annotation eta_pi MH accepted/proposed by trait:\n";
//  for (int t = 0; t < nt; ++t) {
//   Rcpp::Rcout
//   << "trait " << t
//   << ": " << eta_pi_accept(static_cast<arma::uword>(t), 0)
//   << "/" << eta_pi_accept(static_cast<arma::uword>(t), 1)
//   << "\n";
//  }
//
//  Rcpp::Rcout << "Annotation eta_vb MH accepted/proposed by trait:\n";
//  for (int t = 0; t < nt; ++t) {
//   Rcpp::Rcout
//   << "trait " << t
//   << ": " << eta_vb_accept(static_cast<arma::uword>(t), 0)
//   << "/" << eta_vb_accept(static_cast<arma::uword>(t), 1)
//   << "\n";
//  }
//
//  return result;
// }
