// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
// STBLR summary-stat CSR: SBayesRC-style sampler
// =============================================================================
//
// Exported function:
//
//   stblr_cpg_omp_csr_sbayesrc(...)
//
// Model:
//
//   b_it | class_i = c, component z_it = k ~ N(0, vb_t * gamma_k), k > 0
//   b_it | z_it = 0 = 0
//   P(z_it = k | class_i = c) = pi_class[c, k]
//
// where gamma[0] should normally be 0.0 and remaining components are fixed
// variance multipliers, e.g.
//
//   gamma = c(0, 0.01, 0.1, 1)
//
// This is closer to the SBayesRC structure than the annotation-logit STBLR:
// annotations define mutually exclusive classes, and each class has its own
// mixture proportions over variance components.
//
// Notes:
//   * ann_class is length m. It may be 0-based or 1-based via class_one_based.
//   * pi_class_init is C x Kgamma, where C = number of annotation classes and
//     Kgamma = length(gamma). Rows are class-specific mixture probabilities.
//   * If updatePi = TRUE, each trait updates class-specific mixture weights using
//     a Dirichlet posterior with prior alpha_pi_class.
//   * The CSR residual update is identical to your current STBLR:
//       r <- r - X'X_i * diff
//   * Return structure keeps the usual 20 slots:
//       0 bm, 1 dm, 2 wy, 3 r, 4 b, 5 d, 6 o,
//       7 vbs, 8 vgs, 9 ves,
//       10 covb, 11 covg, 12 cove,
//       13 vb, 14 vg, 15 ve,
//       16 final global active pi, 17 posterior mean global active pi,
//       18 flattened final pi_class by trait, length C*Kgamma,
//       19 gamma vector by trait, length Kgamma.
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
//  << "ST SBayesRC flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
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

inline void normalize_prob_row(arma::rowvec& p) {
 double s = 0.0;
 for (arma::uword k = 0; k < p.n_elem; ++k) {
  if (!std::isfinite(p(k)) || p(k) < 0.0) {
   throw std::runtime_error("normalize_prob_row: invalid probability value.");
  }
  s += p(k);
 }

 if (!std::isfinite(s) || s <= 0.0) {
  throw std::runtime_error("normalize_prob_row: row sum must be positive.");
 }

 p /= s;
}

inline void sampleBeta_SBayesRC_ST_csr(
  int i,
  int cls_i,
  const arma::rowvec& pi_class_row,
  const arma::vec& gamma,
  double vb_t,
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
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: invalid ww value.");
 }

 if (!std::isfinite(vb_t) || vb_t <= 0.0) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: invalid vb.");
 }

 const int Kgamma = static_cast<int>(gamma.n_elem);
 if (static_cast<int>(pi_class_row.n_elem) != Kgamma) {
  throw std::runtime_error("sampleBeta_SBayesRC_ST_csr: pi_class_row/gamma length mismatch.");
 }

 const double vei_safe = std::max(vei_i, 1e-300);
 const double score = r(iu) + wi * b(iu);

 std::vector<double> logp(static_cast<std::size_t>(Kgamma));

 for (int k = 0; k < Kgamma; ++k) {
  const double pik = std::max(static_cast<double>(pi_class_row(static_cast<arma::uword>(k))), 1e-300);
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

  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   r(static_cast<arma::uword>(j)) -=
    static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
  }
 }

 b(iu) = b_new;
 comp(iu) = k_new;
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

inline void update_pi_class_dirichlet(
  arma::mat& pi_class,
  const arma::Row<int>& comp,
  const std::vector<int>& ann_class,
  const arma::mat& alpha_pi_class,
  std::mt19937& gen
) {
 const int C = static_cast<int>(pi_class.n_rows);
 const int Kgamma = static_cast<int>(pi_class.n_cols);
 const int m = static_cast<int>(comp.n_elem);

 arma::mat counts(C, Kgamma, arma::fill::zeros);

 for (int i = 0; i < m; ++i) {
  const int c = ann_class[static_cast<std::size_t>(i)];
  const int k = comp(static_cast<arma::uword>(i));

  if (c < 0 || c >= C) {
   throw std::runtime_error("update_pi_class_dirichlet: class index out of range.");
  }

  if (k < 0 || k >= Kgamma) {
   throw std::runtime_error("update_pi_class_dirichlet: component index out of range.");
  }

  counts(static_cast<arma::uword>(c), static_cast<arma::uword>(k)) += 1.0;
 }

 for (int c = 0; c < C; ++c) {
  double sum_g = 0.0;

  for (int k = 0; k < Kgamma; ++k) {
   const double shape = alpha_pi_class(static_cast<arma::uword>(c), static_cast<arma::uword>(k)) +
    counts(static_cast<arma::uword>(c), static_cast<arma::uword>(k));

   if (!std::isfinite(shape) || shape <= 0.0) {
    throw std::runtime_error("update_pi_class_dirichlet: invalid Dirichlet shape.");
   }

   std::gamma_distribution<double> rg(shape, 1.0);
   const double g = std::max(rg(gen), 1e-300);
   pi_class(static_cast<arma::uword>(c), static_cast<arma::uword>(k)) = g;
   sum_g += g;
  }

  if (!std::isfinite(sum_g) || sum_g <= 0.0) {
   throw std::runtime_error("update_pi_class_dirichlet: invalid gamma sum.");
  }

  for (int k = 0; k < Kgamma; ++k) {
   pi_class(static_cast<arma::uword>(c), static_cast<arma::uword>(k)) /= sum_g;
  }
 }
}

inline std::vector<int> convert_ann_class(
  const std::vector<int>& ann_class_in,
  int m,
  int n_classes,
  bool class_one_based
) {
 if (static_cast<int>(ann_class_in.size()) != m) {
  throw std::runtime_error("ann_class must have length m.");
 }

 std::vector<int> out(static_cast<std::size_t>(m));

 for (int i = 0; i < m; ++i) {
  int c = ann_class_in[static_cast<std::size_t>(i)];
  if (class_one_based) c -= 1;

  if (c < 0 || c >= n_classes) {
   throw std::runtime_error(
     "ann_class contains class outside 0..C-1 or 1..C, marker " + std::to_string(i)
   );
  }

  out[static_cast<std::size_t>(i)] = c;
 }

 return out;
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_sbayesrc_annot1(
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
  std::vector<int> ann_class,
  bool class_one_based,
  arma::vec gamma,
  arma::mat pi_class_init,
  arma::mat alpha_pi_class,
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
  int seed
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

 if ((int)ww.size() != nt || (int)b_init.size() != nt ||
     (int)yy.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: inconsistent trait dimensions.");
 }

 if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: priors must be nt x nt.");
 }

 const int Kgamma = static_cast<int>(gamma.n_elem);

 if (Kgamma < 2) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma must have at least two components, including zero.");
 }

 if (!std::isfinite(gamma(0)) || gamma(0) != 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma[0] must be exactly 0.0.");
 }

 for (int k = 1; k < Kgamma; ++k) {
  if (!std::isfinite(gamma(static_cast<arma::uword>(k))) || gamma(static_cast<arma::uword>(k)) <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma[k] must be positive for k > 0.");
  }
 }

 const int C = static_cast<int>(pi_class_init.n_rows);

 if (C <= 0 || static_cast<int>(pi_class_init.n_cols) != Kgamma) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: pi_class_init must be C x length(gamma).");
 }

 if (static_cast<int>(alpha_pi_class.n_rows) != C || static_cast<int>(alpha_pi_class.n_cols) != Kgamma) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_pi_class must be C x length(gamma).");
 }

 for (int c = 0; c < C; ++c) {
  arma::rowvec p = pi_class_init.row(static_cast<arma::uword>(c));
  normalize_prob_row(p);
  pi_class_init.row(static_cast<arma::uword>(c)) = p;

  for (int k = 0; k < Kgamma; ++k) {
   const double a = alpha_pi_class(static_cast<arma::uword>(c), static_cast<arma::uword>(k));
   if (!std::isfinite(a) || a <= 0.0) {
    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_pi_class must be positive.");
   }
  }
 }

 std::vector<int> cls = convert_ann_class(ann_class, m, C, class_one_based);

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

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);
 arma::Mat<int> comp_mat(nt, m, arma::fill::zeros);

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

 STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);

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

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);

 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);

 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi_active(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);

 std::vector<arma::mat> final_pi_class(static_cast<std::size_t>(nt));
 std::vector<arma::mat> mean_pi_class(static_cast<std::size_t>(nt));
 for (int t = 0; t < nt; ++t) {
  final_pi_class[static_cast<std::size_t>(t)] = arma::mat(C, Kgamma, arma::fill::zeros);
  mean_pi_class[static_cast<std::size_t>(t)] = arma::mat(C, Kgamma, arma::fill::zeros);
 }

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
 << "STBLR SBayesRC CSR OpenMP requested threads = "
 << nthreads
 << ", omp_get_max_threads = "
 << omp_get_max_threads()
 << ", num procs = "
 << omp_get_num_procs()
 << "\n";
#endif

 Rcpp::Rcout
 << "STBLR SBayesRC CSR: m=" << m
 << ", nt=" << nt
 << ", classes=" << C
 << ", components=" << Kgamma
 << ", updatePi=" << updatePi
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

   arma::rowvec b_t(m, arma::fill::zeros);
   for (int i = 0; i < m; ++i) {
    b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   }

   arma::rowvec r_t(m, arma::fill::zeros);
   arma::Row<int> comp_t(m, arma::fill::zeros);

   if (use_comp_init) {
    for (int i = 0; i < m; ++i) {
     const int k = static_cast<int>(std::round(comp_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)]));
     if (k < 0 || k >= Kgamma) {
      throw std::runtime_error("comp_init contains component outside 0..Kgamma-1.");
     }
     comp_t(static_cast<arma::uword>(i)) = k;
    }
   } else {
    for (int i = 0; i < m; ++i) {
     comp_t(static_cast<arma::uword>(i)) =
      (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
    }
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) {
     r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
    }

    if (!r_t.is_finite()) {
     throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: r_init contains NaN/Inf.");
    }
   } else {
    rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = 0.0;
   double vei_t = ve_t + adjE * vg_t;

   arma::mat pi_class_t = pi_class_init;
   arma::mat pi_class_accum(C, Kgamma, arma::fill::zeros);

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);

   double nsamples_t = 0.0;

   for (int it = 0; it < nit + nburn; ++it) {
    for (int isort = 0; isort < m; ++isort) {
     const int i = order[static_cast<std::size_t>(isort)];
     const int c = cls[static_cast<std::size_t>(i)];

     sampleBeta_SBayesRC_ST_csr(
      i,
      c,
      pi_class_t.row(static_cast<arma::uword>(c)),
      gamma,
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

    if (updateB) {
     sampleB_SBayesRC_ST_csr(
      m,
      nub,
      vb_t,
      b_t,
      comp_t,
      gamma,
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
     update_pi_class_dirichlet(
      pi_class_t,
      comp_t,
      cls,
      alpha_pi_class,
      gen_t
     );
    }

    vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);

    if (!std::isfinite(vg_t)) {
     throw std::runtime_error("vg became NaN/Inf after computeG. iter=" + std::to_string(it));
    }

    vei_t = ve_t + adjE * vg_t;

    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error(
       "adjusted residual variance vei became invalid. iter=" +
        std::to_string(it) +
        ", vei=" + std::to_string(vei_t)
     );
    }

    double pi_active = 0.0;
    for (int c = 0; c < C; ++c) {
     double row_active = 0.0;
     for (int k = 1; k < Kgamma; ++k) {
      row_active += pi_class_t(static_cast<arma::uword>(c), static_cast<arma::uword>(k));
     }
     pi_active += row_active;
    }
    pi_active /= static_cast<double>(C);

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    pis_t(static_cast<arma::uword>(it)) = pi_active;

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;

     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += (comp_t(iu) > 0) ? 1.0 : 0.0;
     }

     pi_class_accum += pi_class_t;
    }
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;

   bm_t /= nsamples_t;
   dm_t /= nsamples_t;
   pi_class_accum /= nsamples_t;

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
   comp_mat.row(static_cast<arma::uword>(t)) = comp_t;

   vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
   vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
   ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
   pis_mat.row(static_cast<arma::uword>(t)) = pis_t;

   final_vb(static_cast<arma::uword>(t)) = vb_t;
   final_ve(static_cast<arma::uword>(t)) = ve_t;
   final_vg(static_cast<arma::uword>(t)) = vg_t;
   final_pi_active(static_cast<arma::uword>(t)) = pis_t(static_cast<arma::uword>(nit + nburn - 1));
   nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;

   final_pi_class[static_cast<std::size_t>(t)] = pi_class_t;
   mean_pi_class[static_cast<std::size_t>(t)] = pi_class_accum;

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
     "stblr_cpg_omp_csr_sbayesrc failed for trait " +
      std::to_string(t) +
      ": " +
      errors[static_cast<std::size_t>(t)]
   );
  }
 }

 std::vector<std::vector<std::vector<double>>> result(20);

 for (int k = 0; k < 20; ++k) {
  result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
 }

 for (int t = 0; t < nt; ++t) {
  result[0][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[1][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[2][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[3][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[4][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[5][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[6][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));

  result[7][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
  result[8][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
  result[9][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));

  result[10][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[11][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[12][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[13][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[14][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[15][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));

  result[16][static_cast<std::size_t>(t)].resize(2);
  result[17][static_cast<std::size_t>(t)].resize(2);

  result[18][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(C * Kgamma));
  result[19][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(Kgamma));
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
   result[5][ts][is] = static_cast<double>(comp_mat(tu, iu));
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

  result[16][ts][0] = 1.0 - final_pi_active(static_cast<arma::uword>(t));
  result[16][ts][1] = final_pi_active(static_cast<arma::uword>(t));

  double mean_pi = 0.0;
  int npi = 0;

  for (int it = nburn; it < nit + nburn; ++it) {
   mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
   ++npi;
  }

  if (npi > 0) mean_pi /= static_cast<double>(npi);
  else mean_pi = final_pi_active(static_cast<arma::uword>(t));

  result[17][ts][0] = 1.0 - mean_pi;
  result[17][ts][1] = mean_pi;

  for (int c = 0; c < C; ++c) {
   for (int k = 0; k < Kgamma; ++k) {
    const std::size_t pos = static_cast<std::size_t>(c * Kgamma + k);
    result[18][ts][pos] = mean_pi_class[ts](static_cast<arma::uword>(c), static_cast<arma::uword>(k));
   }
  }

  for (int k = 0; k < Kgamma; ++k) {
   result[19][ts][static_cast<std::size_t>(k)] = gamma(static_cast<arma::uword>(k));
  }
 }

 return result;
}
