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
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

// =============================================================================
// STBLR summary-stat CSR with annotation/marker-informed priors
// =============================================================================
//
// This file defines a new exported function:
//
//   stblr_cpg_omp_csr_prior(...)
//
// It follows the current 22-slot CSR BLR return structure and adds:
//
//   use_pi_marker
//   pi_marker[t][i]          absolute prior inclusion probability for marker i,
//                            trait t. If false, global pi_t is used.
//
//   use_vb_multiplier
//   vb_multiplier[t][i]      relative prior variance multiplier. If false, 1.0.
//
// Effective marker prior:
//
//   Global pi_t has explicit Beta(pi_prior_a, pi_prior_b) prior when updatePi=TRUE.
//
//   d_it ~ Bernoulli(pi_it)
//   b_it | d_it = 1 ~ N(0, vb_t * vb_multiplier_it)
//
// If use_vb_multiplier = TRUE, sampleB updates the global scale vb_t using
// sum_i b_i^2 / multiplier_i among active markers, which is the correct scalar
// scale update for b_i ~ N(0, vb_t * multiplier_i).
//
// Return layout:
//
//   0  bm
//   1  dm
//   2  wy
//   3  r
//   4  b
//   5  d
//   6  o
//   7  vbs
//   8  vgs
//   9  ves
//   10 covb
//   11 covg
//   12 cove
//   13 vb
//   14 vg
//   15 ve
//   16 final pi
//   17 posterior mean pi
//   18 diagnostics: currently zeros/reserved, length 4
//   19 diagnostics: nsamples, n, length 2
//   20 vle
//   21 vld
//
// =============================================================================

// -----------------------------------------------------------------------------
// ST-specific LD structure: flat symmetric CSR
// Stores pre-scaled X_i'X_j, not raw LD correlation.
// Disk input is expected to be upper-triangular or otherwise non-symmetric CSR
// with raw correlations r_ij. This builder symmetrizes it.
// -----------------------------------------------------------------------------

// struct STLDCSR {
//  std::vector<uint64_t> ptr;  // length m + 1
//  std::vector<int> idx;       // neighbor marker index
//  std::vector<float> xij;     // pre-scaled X_i'X_j
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
//  // First pass: count symmetric degrees.
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
//  // Second pass: fill symmetric flat CSR.
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
//  // Validate fill counts.
//  for (int i = 0; i < m; ++i) {
//   if (offset[static_cast<std::size_t>(i)] != ld.ptr[static_cast<std::size_t>(i + 1)]) {
//    throw std::runtime_error("Internal LD CSR fill-count mismatch.");
//   }
//  }
//
//  Rcpp::Rcout
//  << "ST prior flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
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

// -----------------------------------------------------------------------------
// Single-trait BayesC marker update with marker-specific pi and vb multiplier
// -----------------------------------------------------------------------------

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

 // score = x_i' residual_without_i
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

// -----------------------------------------------------------------------------
// Single-trait variance and pi updates
// -----------------------------------------------------------------------------

inline void sampleB_ST_csr_prior(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  const arma::rowvec& vb_multiplier,
  bool use_vb_multiplier,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb_scaled = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);

  if (d(iu) > 0) {
   double mult = 1.0;

   if (use_vb_multiplier) {
    mult = vb_multiplier(iu);

    if (!std::isfinite(mult) || mult <= 0.0) {
     throw std::runtime_error("sampleB_ST_csr_prior: invalid vb multiplier.");
    }
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

inline double computeLE_ST_csr_prior(
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

inline void samplePi_ST_prior(
  const arma::Row<int>& d,
  std::vector<double>& pi,
  double pi_prior_a,
  double pi_prior_b,
  std::mt19937& gen
) {
 // pi[1] is the inclusion probability; pi[0] is exclusion probability.
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

inline void samplePi_ST(
  const arma::Row<int>& d,
  std::vector<double>& pi,
  std::mt19937& gen
) {
 double c0 = 1.0;
 double c1 = 1.0;

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

inline void validate_prior_matrix_input(
  const std::vector<std::vector<double>>& x,
  const std::string& name,
  int nt,
  int m,
  bool required
) {
 if (!required) return;

 if (static_cast<int>(x.size()) != nt) {
  throw std::runtime_error(name + " must have length nt when enabled.");
 }

 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(x[static_cast<std::size_t>(t)].size()) != m) {
   throw std::runtime_error(name + "[t] must have length m when enabled.");
  }
 }
}

// -----------------------------------------------------------------------------
// Main exported function: STBLR over traits with marker-specific priors
// -----------------------------------------------------------------------------

std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_prior_single(
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
  bool use_pi_marker,
  std::vector<std::vector<double>> pi_marker,
  bool use_vb_multiplier,
  std::vector<std::vector<double>> vb_multiplier,
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
  int seed
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: nt must be positive.");
 }

 const int m = static_cast<int>(wy[0].size());

 if (m <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: m must be positive.");
 }

 if (nit <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: nit must be positive.");
 }

 if (nburn < 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: nburn must be non-negative.");
 }

 if (nthin <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: nthin must be positive.");
 }

 if (!std::isfinite(pi_prior_a) || pi_prior_a <= 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: pi_prior_a must be finite and positive.");
 }

 if (!std::isfinite(pi_prior_b) || pi_prior_b <= 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: pi_prior_b must be finite and positive.");
 }

 if ((int)ww.size() != nt || (int)b_init.size() != nt ||
     (int)yy.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: inconsistent trait dimensions.");
 }

 if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: priors must be nt x nt.");
 }

 if (pi.size() != 2) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: pi must have length 2, c(pi0, pi1).");
 }

 validate_prior_matrix_input(pi_marker, "pi_marker", nt, m, use_pi_marker);
 validate_prior_matrix_input(vb_multiplier, "vb_multiplier", nt, m, use_vb_multiplier);

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m ||
      (int)ww[t].size() != m ||
      (int)b_init[t].size() != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_prior: inconsistent marker dimensions.");
  }
 }

 if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: B must be nt x nt.");
 }

 if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: E must be nt x nt.");
 }

 if (use_r_init) {
  if (static_cast<int>(r_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_prior_state: r_init must have length nt when use_r_init = true.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(r_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_prior_state: each r_init[t] must have length m.");
   }
  }
 }

 if (use_d_init) {
  if (static_cast<int>(d_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_prior_state: d_init must have length nt when use_d_init = true.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(d_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_prior_state: each d_init[t] must have length m.");
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

 arma::mat pi_marker_mat(nt, m, arma::fill::zeros);
 arma::mat vb_multiplier_mat(nt, m, arma::fill::ones);

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

   if (use_pi_marker) {
    const double pij = pi_marker[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)];

    if (!std::isfinite(pij) || pij <= 0.0 || pij >= 1.0) {
     throw std::runtime_error(
       "stblr_cpg_omp_csr_prior: pi_marker contains invalid value at trait " +
        std::to_string(t) + ", marker " + std::to_string(i)
     );
    }

    pi_marker_mat(tu, iu) = pij;
   }

   if (use_vb_multiplier) {
    const double mult = vb_multiplier[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)];

    if (!std::isfinite(mult) || mult <= 0.0) {
     throw std::runtime_error(
       "stblr_cpg_omp_csr_prior: vb_multiplier contains invalid value at trait " +
        std::to_string(t) + ", marker " + std::to_string(i)
     );
    }

    vb_multiplier_mat(tu, iu) = mult;
   }
  }

  if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_prior: priors must be nt x nt.");
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
     "stblr_cpg_omp_csr_prior: current shared-LD scaling assumes equal n across traits."
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
      "stblr_cpg_omp_csr_prior: ww contains invalid value before LD pre-scaling."
    );
   }

   if (std::abs(w0 - wt) > tol) {
    throw std::runtime_error(
      "stblr_cpg_omp_csr_prior: ww differs across traits; pre-scaled shared ST LD is invalid."
    );
   }
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);

 for (int i = 0; i < m; ++i) {
  const double wi = ww_mat(0, static_cast<arma::uword>(i));
  if (!std::isfinite(wi) || wi <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_prior: ww contains invalid value in trait 0.");
  }
  xx[static_cast<std::size_t>(i)] = wi;
 }

 STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);

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
           [&](int a, int b) {
            return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
           });

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
 omp_set_num_threads(nthreads);

 Rcpp::Rcout
 << "STBLR prior CSR OpenMP requested threads = "
 << nthreads
 << ", omp_get_max_threads = "
 << omp_get_max_threads()
 << ", num procs = "
 << omp_get_num_procs()
 << "\n";
#endif

 Rcpp::Rcout
 << "STBLR prior CSR: use_pi_marker=" << use_pi_marker
 << ", use_vb_multiplier=" << use_vb_multiplier
 << ", updatePi=" << updatePi
 << ", updateB=" << updateB
 << ", pi_prior_a=" << pi_prior_a
 << ", pi_prior_b=" << pi_prior_b
 << "\n";

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
   arma::rowvec pi_marker_t = pi_marker_mat.row(static_cast<arma::uword>(t));
   arma::rowvec vb_multiplier_t = vb_multiplier_mat.row(static_cast<arma::uword>(t));

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
     throw std::runtime_error("stblr_cpg_omp_csr_prior_state: r_init contains NaN/Inf.");
    }
   } else {
    rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_ST_csr_prior(m, b_t, ww_t, n[t]);
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

   for (int it = 0; it < nit + nburn; ++it) {

    // -------------------------------------------------------
    // Marker updates
    // -------------------------------------------------------
    for (int isort = 0; isort < m; ++isort) {
     const int i = order[static_cast<std::size_t>(isort)];
     const arma::uword iu = static_cast<arma::uword>(i);

     const double pi1_i = use_pi_marker ? pi_marker_t(iu) : pi_t[1];
     const double vb_mult_i = use_vb_multiplier ? vb_multiplier_t(iu) : 1.0;

     sampleBetaC_ST_csr_prior(
      i,
      pi1_i,
      vb_t,
      vb_mult_i,
      vei_t,
      ww_t,
      r_t,
      b_t,
      d_t,
      ld,
      gen_t
     );
    }

    // -------------------------------------------------------
    // Variance updates
    // -------------------------------------------------------
    if (updateB) {
     sampleB_ST_csr_prior(
      m,
      nub,
      vb_t,
      b_t,
      d_t,
      vb_multiplier_t,
      use_vb_multiplier,
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
      rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
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
     samplePi_ST_prior(d_t, pi_t, pi_prior_a, pi_prior_b, gen_t);

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

    vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
    vle_t = computeLE_ST_csr_prior(m, b_t, ww_t, n[t]);
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
   final_pi(static_cast<arma::uword>(t)) = pi_t[1];
   final_vle(static_cast<arma::uword>(t)) = vle_t;
   final_vld(static_cast<arma::uword>(t)) = vld_t;
   nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;

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
     "stblr_cpg_omp_csr_prior failed for trait " +
      std::to_string(t) +
      ": " +
      errors[static_cast<std::size_t>(t)]
   );
  }
 }

 // --------------------------------------------------------------------------
 // Build result with same style as MT output
 // --------------------------------------------------------------------------

 std::vector<std::vector<std::vector<double>>> result(22);

 for (int k = 0; k < 22; ++k) {
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

  result[16][static_cast<std::size_t>(t)].resize(2);                                     // final pi trace/reporting
  result[17][static_cast<std::size_t>(t)].resize(2);                                     // posterior mean pi trace/reporting

  result[18][static_cast<std::size_t>(t)].resize(4);
  result[19][static_cast<std::size_t>(t)].resize(2);

  result[20][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn)); // VLE
  result[21][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn)); // VLD
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

  // Diagnostics slot, kept compatible with current CSR formatters.
  // [0] reserved/log_cpo placeholder for CSR prior sampler
  // [1] reserved/mean_log_cpo placeholder for CSR prior sampler
  // [2] runtime seconds
  // [3] runtime seconds, duplicate for compatibility
  result[18][ts][0] = 0.0;
  result[18][ts][1] = 0.0;
  result[18][ts][2] = trait_seconds[static_cast<std::size_t>(t)];
  result[18][ts][3] = trait_seconds[static_cast<std::size_t>(t)];

  result[19][ts][0] = nsamples_vec(static_cast<arma::uword>(t));
  result[19][ts][1] = static_cast<double>(n[t]);
 }

 return result;
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_prior(
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
  bool use_pi_marker,
  std::vector<std::vector<double>> pi_marker,
  bool use_vb_multiplier,
  std::vector<std::vector<double>> vb_multiplier,
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
  Rcpp::Nullable<Rcpp::IntegerVector> chain_seeds = R_NilValue
) {
 if (nchains <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: nchains must be positive.");
 }
 std::vector<int> chain_seeds_vec;
 if (chain_seeds.isNotNull()) {
  Rcpp::IntegerVector chain_seeds_r(chain_seeds);
  chain_seeds_vec = Rcpp::as<std::vector<int>>(chain_seeds_r);
 }
 if (!chain_seeds_vec.empty() && static_cast<int>(chain_seeds_vec.size()) != nchains) {
  throw std::runtime_error("stblr_cpg_omp_csr_prior: chain_seeds must have length nchains.");
 }

 std::vector<std::vector<std::vector<double>>> chain0;
 std::vector<std::vector<std::vector<double>>> out;
 std::vector<std::vector<std::vector<double>>> chain_dm;
 std::vector<std::vector<std::vector<double>>> chain_bm;
 std::vector<std::vector<std::vector<double>>> chain_diag;
 std::vector<std::vector<std::vector<double>>> sumsq;
 std::vector<std::vector<std::vector<double>>> minv;
 std::vector<std::vector<std::vector<double>>> maxv;

 for (int chain = 0; chain < nchains; ++chain) {
  int chain_seed = seed;
  if (!chain_seeds_vec.empty()) {
   chain_seed = chain_seeds_vec[static_cast<std::size_t>(chain)];
  } else if (nchains > 1) {
   chain_seed = seed + 9176 * (chain + 1);
  }

  std::vector<std::vector<std::vector<double>>> raw =
   stblr_cpg_omp_csr_prior_single(
    wy, ww, yy, b_init, d_init, use_d_init, r_init, use_r_init,
    rebuild_r_before_updateE, ld_prefix, B, E, ssb_prior, sse_prior, pi,
    use_pi_marker, pi_marker, use_vb_multiplier, vb_multiplier, nub, nue,
    updateB, updateE, updatePi, adjE, n, nit, nburn, nthin, pi_prior_a,
    pi_prior_b, ncores, chain_seed
   );

  if (chain == 0) {
   chain0 = raw;
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
   chain_dm.resize(raw[1].size());
   chain_bm.resize(raw[0].size());
   chain_diag.resize(raw[0].size());
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

  if (keep_chains) {
   for (std::size_t t = 0; t < raw[0].size(); ++t) {
    chain_dm[t].insert(chain_dm[t].end(), raw[1][t].begin(), raw[1][t].end());
    chain_bm[t].insert(chain_bm[t].end(), raw[0][t].begin(), raw[0][t].end());
    chain_diag[t].insert(chain_diag[t].end(), 4, 0.0);
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

 const bool return_chain_summaries = (nchains > 1) || keep_chains;
 if (return_chain_summaries) {
  std::vector<std::vector<std::vector<double>>> extended(keep_chains ? 32 : 29);
  for (std::size_t k = 0; k < out.size(); ++k) extended[k] = out[k];
  extended[22].resize(out[0].size());
  for (std::size_t t = 0; t < out[0].size(); ++t) {
   extended[22][t].resize(static_cast<std::size_t>(nit + nburn));
   for (std::size_t i = 0; i < extended[22][t].size(); ++i) {
    extended[22][t][i] = (chain0.size() > 16 && chain0[16][t].size() > 1) ? chain0[16][t][1] : 0.0;
   }
  }
  for (int slot = 23; slot <= 28; ++slot) extended[static_cast<std::size_t>(slot)].resize(out[0].size());
  for (std::size_t t = 0; t < out[0].size(); ++t) {
   const std::size_t m_t = out[0][t].size();
   for (int slot = 23; slot <= 28; ++slot) extended[static_cast<std::size_t>(slot)][t].resize(m_t);
   for (std::size_t i = 0; i < m_t; ++i) {
    const double bm_mean = out[0][t][i];
    const double dm_mean = out[1][t][i];
    double bm_sd = 0.0;
    double dm_sd = 0.0;
    if (nchains > 1) {
     bm_sd = std::sqrt(std::max(0.0, (sumsq[0][t][i] - nchains * bm_mean * bm_mean) / static_cast<double>(nchains - 1)));
     dm_sd = std::sqrt(std::max(0.0, (sumsq[1][t][i] - nchains * dm_mean * dm_mean) / static_cast<double>(nchains - 1)));
    }
    extended[23][t][i] = bm_sd;
    extended[24][t][i] = minv[0][t][i];
    extended[25][t][i] = maxv[0][t][i];
    extended[26][t][i] = dm_sd;
    extended[27][t][i] = minv[1][t][i];
    extended[28][t][i] = maxv[1][t][i];
   }
  }
  if (keep_chains) {
   extended[29] = chain_dm;
   extended[30] = chain_bm;
   extended[31] = chain_diag;
  }
  return extended;
 }

 return out;
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
// // STBLR summary-stat CSR with annotation/marker-informed priors
// // =============================================================================
// //
// // This file defines a new exported function:
// //
// //   stblr_cpg_omp_csr_prior(...)
// //
// // It keeps the same 20-slot return structure as stblr_cpg_omp_csr(), but adds:
// //
// //   use_pi_marker
// //   pi_marker[t][i]          absolute prior inclusion probability for marker i,
// //                            trait t. If false, global pi_t is used.
// //
// //   use_vb_multiplier
// //   vb_multiplier[t][i]      relative prior variance multiplier. If false, 1.0.
// //
// // Effective marker prior:
// //
// //   d_it ~ Bernoulli(pi_it)
// //   b_it | d_it = 1 ~ N(0, vb_t * vb_multiplier_it)
// //
// // If use_vb_multiplier = TRUE, sampleB updates the global scale vb_t using
// // sum_i b_i^2 / multiplier_i among active markers, which is the correct scalar
// // scale update for b_i ~ N(0, vb_t * multiplier_i).
// //
// // =============================================================================
//
// // -----------------------------------------------------------------------------
// // ST-specific LD structure: flat symmetric CSR
// // Stores pre-scaled X_i'X_j, not raw LD correlation.
// // Disk input is expected to be upper-triangular or otherwise non-symmetric CSR
// // with raw correlations r_ij. This builder symmetrizes it.
// // -----------------------------------------------------------------------------
//
// // struct STLDCSR {
// //  std::vector<uint64_t> ptr;  // length m + 1
// //  std::vector<int> idx;       // neighbor marker index
// //  std::vector<float> xij;     // pre-scaled X_i'X_j
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
// //  // First pass: count symmetric degrees.
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
// //  // Second pass: fill symmetric flat CSR.
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
// //  // Validate fill counts.
// //  for (int i = 0; i < m; ++i) {
// //   if (offset[static_cast<std::size_t>(i)] != ld.ptr[static_cast<std::size_t>(i + 1)]) {
// //    throw std::runtime_error("Internal LD CSR fill-count mismatch.");
// //   }
// //  }
// //
// //  Rcpp::Rcout
// //  << "ST prior flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
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
//
// // -----------------------------------------------------------------------------
// // Single-trait BayesC marker update with marker-specific pi and vb multiplier
// // -----------------------------------------------------------------------------
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
//
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
//
//  const double vbi = std::max(vb_t * vb_mult_i, 1e-12);
//  const double vei_safe = std::max(vei_i, 1e-300);
//
//  // score = x_i' residual_without_i
//  const double score = r(iu) + wi * b(iu);
//
//  const double denom = std::max(vei_safe + wi * vbi, 1e-300);
//
//  const double logBF =
//   0.5 * std::log(vei_safe / denom)
//   + 0.5 * score * score * vbi / (vei_safe * denom);
//
//  const double logp1 = std::log(pi1_i) + logBF;
//  const double logp0 = std::log(pi0_i);
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
//   const double lhs = wi + vei_safe / vbi;
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
// // Single-trait variance and pi updates
// // -----------------------------------------------------------------------------
//
// inline void sampleB_ST_csr_prior(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   const arma::rowvec& vb_multiplier,
//   bool use_vb_multiplier,
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
//    double mult = 1.0;
//
//    if (use_vb_multiplier) {
//     mult = vb_multiplier(iu);
//
//     if (!std::isfinite(mult) || mult <= 0.0) {
//      throw std::runtime_error("sampleB_ST_csr_prior: invalid vb multiplier.");
//     }
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
// inline void validate_prior_matrix_input(
//   const std::vector<std::vector<double>>& x,
//   const std::string& name,
//   int nt,
//   int m,
//   bool required
// ) {
//  if (!required) return;
//
//  if (static_cast<int>(x.size()) != nt) {
//   throw std::runtime_error(name + " must have length nt when enabled.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(x[static_cast<std::size_t>(t)].size()) != m) {
//    throw std::runtime_error(name + "[t] must have length m when enabled.");
//   }
//  }
// }
//
// // -----------------------------------------------------------------------------
// // Main exported function: STBLR over traits with marker-specific priors
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_prior(
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
//   bool use_pi_marker,
//   std::vector<std::vector<double>> pi_marker,
//   bool use_vb_multiplier,
//   std::vector<std::vector<double>> vb_multiplier,
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
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: nt must be positive.");
//  }
//
//  const int m = static_cast<int>(wy[0].size());
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: m must be positive.");
//  }
//
//  if (nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: nit must be positive.");
//  }
//
//  if (nburn < 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: nburn must be non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: nthin must be positive.");
//  }
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt ||
//      (int)yy.size() != nt || (int)n.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: inconsistent trait dimensions.");
//  }
//
//  if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: priors must be nt x nt.");
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: pi must have length 2, c(pi0, pi1).");
//  }
//
//  validate_prior_matrix_input(pi_marker, "pi_marker", nt, m, use_pi_marker);
//  validate_prior_matrix_input(vb_multiplier, "vb_multiplier", nt, m, use_vb_multiplier);
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m ||
//       (int)ww[t].size() != m ||
//       (int)b_init[t].size() != m) {
//    throw std::runtime_error("stblr_cpg_omp_csr_prior: inconsistent marker dimensions.");
//   }
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_prior: E must be nt x nt.");
//  }
//
//  if (use_r_init) {
//   if (static_cast<int>(r_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_prior_state: r_init must have length nt when use_r_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(r_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_prior_state: each r_init[t] must have length m.");
//    }
//   }
//  }
//
//  if (use_d_init) {
//   if (static_cast<int>(d_init.size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_prior_state: d_init must have length nt when use_d_init = true.");
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(d_init[t].size()) != m) {
//     throw std::runtime_error("stblr_cpg_omp_csr_prior_state: each d_init[t] must have length m.");
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
//  arma::mat pi_marker_mat(nt, m, arma::fill::zeros);
//  arma::mat vb_multiplier_mat(nt, m, arma::fill::ones);
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
//
//    if (use_pi_marker) {
//     const double pij = pi_marker[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)];
//
//     if (!std::isfinite(pij) || pij <= 0.0 || pij >= 1.0) {
//      throw std::runtime_error(
//        "stblr_cpg_omp_csr_prior: pi_marker contains invalid value at trait " +
//         std::to_string(t) + ", marker " + std::to_string(i)
//      );
//     }
//
//     pi_marker_mat(tu, iu) = pij;
//    }
//
//    if (use_vb_multiplier) {
//     const double mult = vb_multiplier[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)];
//
//     if (!std::isfinite(mult) || mult <= 0.0) {
//      throw std::runtime_error(
//        "stblr_cpg_omp_csr_prior: vb_multiplier contains invalid value at trait " +
//         std::to_string(t) + ", marker " + std::to_string(i)
//      );
//     }
//
//     vb_multiplier_mat(tu, iu) = mult;
//    }
//   }
//
//   if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_prior: priors must be nt x nt.");
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
//      "stblr_cpg_omp_csr_prior: current shared-LD scaling assumes equal n across traits."
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
//       "stblr_cpg_omp_csr_prior: ww contains invalid value before LD pre-scaling."
//     );
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error(
//       "stblr_cpg_omp_csr_prior: ww differs across traits; pre-scaled shared ST LD is invalid."
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
//    throw std::runtime_error("stblr_cpg_omp_csr_prior: ww contains invalid value in trait 0.");
//   }
//   xx[static_cast<std::size_t>(i)] = wi;
//  }
//
//  STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
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
//            [&](int a, int b) {
//             return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
//            });
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
//  nthreads = std::max(1, std::min(ncores, nt));
//  omp_set_num_threads(nthreads);
//
//  Rcpp::Rcout
//  << "STBLR prior CSR OpenMP requested threads = "
//  << nthreads
//  << ", omp_get_max_threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
//  Rcpp::Rcout
//  << "STBLR prior CSR: use_pi_marker=" << use_pi_marker
//  << ", use_vb_multiplier=" << use_vb_multiplier
//  << ", updatePi=" << updatePi
//  << ", updateB=" << updateB
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
//    arma::rowvec pi_marker_t = pi_marker_mat.row(static_cast<arma::uword>(t));
//    arma::rowvec vb_multiplier_t = vb_multiplier_mat.row(static_cast<arma::uword>(t));
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
//      throw std::runtime_error("stblr_cpg_omp_csr_prior_state: r_init contains NaN/Inf.");
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
//      const arma::uword iu = static_cast<arma::uword>(i);
//
//      const double pi1_i = use_pi_marker ? pi_marker_t(iu) : pi_t[1];
//      const double vb_mult_i = use_vb_multiplier ? vb_multiplier_t(iu) : 1.0;
//
//      sampleBetaC_ST_csr_prior(
//       i,
//       pi1_i,
//       vb_t,
//       vb_mult_i,
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
//      sampleB_ST_csr_prior(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       vb_multiplier_t,
//       use_vb_multiplier,
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
//      "stblr_cpg_omp_csr_prior failed for trait " +
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
//   result[16][static_cast<std::size_t>(t)].resize(2);                                     // final pi trace/reporting
//   result[17][static_cast<std::size_t>(t)].resize(2);                                     // posterior mean pi trace/reporting
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
