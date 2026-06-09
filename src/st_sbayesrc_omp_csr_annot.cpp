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
// STBLR summary-stat CSR: SBayesRC-like comparison sampler
// =============================================================================
//
// Exported function:
//
//   stblr_cpg_omp_csr_sbayesrc(...)
//
// Purpose:
//   A practical SBayesRC-like comparison implementation for the same sparse CSR
//   summary-statistic LD representation used by stblr_cpg_omp_csr().
//
// Main idea:
//   Each marker belongs to one primary annotation class c_i = 0, ..., C - 1.
//   Each marker effect belongs to one mixture component k = 0, ..., Kmix - 1.
//
//   k = 0 is the null component with b_i = 0.
//   k > 0 are non-null normal components:
//
//     b_i | gamma_i = k ~ N(0, vb_t * mixture_var[k])
//
//   Annotation class c has its own mixture probabilities:
//
//     P(gamma_i = k | c_i = c) = pi_class[c, k]
//
//   The pi_class rows are updated by Dirichlet counts, analogous to the
//   annotation-class enrichment idea in SBayesRC.
//
// Important differences from a full published SBayesRC implementation:
//   * This is a compact comparison sampler, not a line-for-line clone.
//   * It supports one primary annotation class per marker. For overlapping
//     annotations, create mutually exclusive classes before calling this function
//     or use your annotation-logit sampler.
//   * It uses the same scalar adjusted residual variance vei_t = ve_t + adjE*vg_t
//     as your current STBLR code.
//   * It estimates annotation-class-specific mixture probabilities, but not
//     annotation-specific LD matrices or full BayesR hyperparameter schemes.
//
// Return structure:
//   Same 20-slot return style as your current STBLR functions.
//
//   result[0]  bm
//   result[1]  dm = posterior non-null probability
//   result[2]  wy
//   result[3]  r
//   result[4]  b
//   result[5]  d = final non-null indicator
//   result[6]  marker index
//   result[7]  vbs trace
//   result[8]  vgs trace
//   result[9]  ves trace
//   result[10:15] diagonal covariance/final variance placeholders
//   result[16] global final non-null pi summary, length 2
//   result[17] global posterior mean non-null pi summary, length 2
//   result[18] flattened posterior mean pi_class, length C * Kmix
//   result[19] flattened final pi_class, length C * Kmix
//
// Suggested model choices:
//   mixture_var = c(0, 0.01, 0.1, 1.0)
//   alpha_class = matrix(1, C, length(mixture_var))
//
// =============================================================================

// -----------------------------------------------------------------------------
// ST-specific LD structure: flat symmetric CSR
// Stores pre-scaled X_i'X_j, not raw LD correlation.
// Disk input is expected to be upper-triangular/non-symmetric CSR with raw r_ij.
// This builder symmetrizes and stores X_i'X_j as float.
// -----------------------------------------------------------------------------

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
//  << "SBayesRC-like ST flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
//  << ", symmetric nnz=" << static_cast<double>(nnz_sym)
//  << ", max_abs_rij=" << max_abs_rij
//  << ", max_abs_xij=" << max_abs_xij
//  << "\n";
//
//  return ld;
// }
//
// // -----------------------------------------------------------------------------
// // Helpers
// // -----------------------------------------------------------------------------
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

inline double safe_positive(double x, const char* name) {
 if (!std::isfinite(x) || x <= 0.0) {
  throw std::runtime_error(std::string(name) + " must be positive and finite.");
 }
 return x;
}

inline arma::mat normalize_rows_to_simplex(
  arma::mat X,
  double eps = 1e-12
) {
 for (arma::uword i = 0; i < X.n_rows; ++i) {
  double s = 0.0;

  for (arma::uword j = 0; j < X.n_cols; ++j) {
   if (!std::isfinite(X(i, j)) || X(i, j) < 0.0) {
    throw std::runtime_error("normalize_rows_to_simplex: invalid entry.");
   }
   X(i, j) = std::max(X(i, j), eps);
   s += X(i, j);
  }

  if (!std::isfinite(s) || s <= 0.0) {
   throw std::runtime_error("normalize_rows_to_simplex: invalid row sum.");
  }

  for (arma::uword j = 0; j < X.n_cols; ++j) {
   X(i, j) /= s;
  }
 }

 return X;
}

inline void validate_class_index(
  const std::vector<int>& annotation_class,
  int m,
  int C
) {
 if (static_cast<int>(annotation_class.size()) != m) {
  throw std::runtime_error("annotation_class must have length m.");
 }

 for (int i = 0; i < m; ++i) {
  if (annotation_class[static_cast<std::size_t>(i)] < 0 ||
      annotation_class[static_cast<std::size_t>(i)] >= C) {
   throw std::runtime_error(
     "annotation_class must be 0-based integers in 0:(C-1). Invalid marker " +
      std::to_string(i)
   );
  }
 }
}

inline int sample_categorical_from_logweights(
  const std::vector<double>& logw,
  std::mt19937& gen
) {
 const int K = static_cast<int>(logw.size());
 if (K <= 0) {
  throw std::runtime_error("sample_categorical_from_logweights: K must be positive.");
 }

 double max_logw = -std::numeric_limits<double>::infinity();
 for (int k = 0; k < K; ++k) {
  max_logw = std::max(max_logw, logw[static_cast<std::size_t>(k)]);
 }

 if (!std::isfinite(max_logw)) {
  throw std::runtime_error("sample_categorical_from_logweights: all log weights invalid.");
 }

 std::vector<double> w(static_cast<std::size_t>(K), 0.0);
 double sw = 0.0;

 for (int k = 0; k < K; ++k) {
  w[static_cast<std::size_t>(k)] = std::exp(logw[static_cast<std::size_t>(k)] - max_logw);
  sw += w[static_cast<std::size_t>(k)];
 }

 if (!std::isfinite(sw) || sw <= 0.0) {
  throw std::runtime_error("sample_categorical_from_logweights: invalid weight sum.");
 }

 std::uniform_real_distribution<double> runif(0.0, sw);
 const double u = runif(gen);
 double cs = 0.0;

 for (int k = 0; k < K; ++k) {
  cs += w[static_cast<std::size_t>(k)];
  if (u <= cs) return k;
 }

 return K - 1;
}

inline arma::rowvec sample_dirichlet_row(
  const arma::rowvec& alpha,
  std::mt19937& gen
) {
 arma::rowvec out(alpha.n_elem, arma::fill::zeros);
 double s = 0.0;

 for (arma::uword k = 0; k < alpha.n_elem; ++k) {
  const double a = safe_positive(alpha(k), "Dirichlet alpha");
  std::gamma_distribution<double> rgamma(a, 1.0);
  out(k) = std::max(rgamma(gen), 1e-300);
  s += out(k);
 }

 if (!std::isfinite(s) || s <= 0.0) {
  throw std::runtime_error("sample_dirichlet_row: invalid gamma sum.");
 }

 out /= s;
 return out;
}

inline void update_pi_class_dirichlet(
  arma::mat& pi_class,
  const arma::Mat<int>& gamma_by_trait_one_trait,
  const std::vector<int>& annotation_class,
  const arma::mat& alpha_class,
  int C,
  int Kmix,
  std::mt19937& gen
) {
 arma::mat counts(C, Kmix, arma::fill::zeros);
 const int m = gamma_by_trait_one_trait.n_cols;

 for (int i = 0; i < m; ++i) {
  const int c = annotation_class[static_cast<std::size_t>(i)];
  const int g = gamma_by_trait_one_trait(0, static_cast<arma::uword>(i));

  if (g < 0 || g >= Kmix) {
   throw std::runtime_error("update_pi_class_dirichlet: invalid mixture state.");
  }

  counts(static_cast<arma::uword>(c), static_cast<arma::uword>(g)) += 1.0;
 }

 for (int c = 0; c < C; ++c) {
  arma::rowvec alpha_post = alpha_class.row(static_cast<arma::uword>(c)) +
   counts.row(static_cast<arma::uword>(c));
  pi_class.row(static_cast<arma::uword>(c)) = sample_dirichlet_row(alpha_post, gen);
 }
}

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
//
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

inline void sampleB_ST_sbayesrc(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& gamma,
  const arma::vec& mixture_var,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb_scaled = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const int g = gamma(iu);

  if (g > 0) {
   const double mult = mixture_var(static_cast<arma::uword>(g));

   if (!std::isfinite(mult) || mult <= 0.0) {
    throw std::runtime_error("sampleB_ST_sbayesrc: invalid non-null mixture variance multiplier.");
   }

   ssb_scaled += b(iu) * b(iu) / mult;
   dfb += 1.0;
  }
 }

 const double scale = ssb_scaled + nub * ssb_prior;

 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleB_ST_sbayesrc: invalid scale.");
 }

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);

 vb = std::max(scale / chi2, 1e-12);
}

// -----------------------------------------------------------------------------
// Marker update: SBayesRC-like mixture update
// -----------------------------------------------------------------------------

inline void sampleBeta_ST_sbayesrc(
  int i,
  int ann_class_i,
  const arma::rowvec& pi_class_row,
  const arma::vec& mixture_var,
  double vb_t,
  double vei_t,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& gamma,
  const STLDCSR& ld,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double wi = ww(iu);

 if (!std::isfinite(wi) || wi <= 0.0) {
  throw std::runtime_error("sampleBeta_ST_sbayesrc: invalid ww value.");
 }

 if (!std::isfinite(vb_t) || vb_t <= 0.0) {
  throw std::runtime_error("sampleBeta_ST_sbayesrc: invalid vb.");
 }

 const int Kmix = static_cast<int>(mixture_var.n_elem);
 if (Kmix < 2) {
  throw std::runtime_error("sampleBeta_ST_sbayesrc: mixture_var must have at least null + one non-null component.");
 }

 if (mixture_var(0) != 0.0) {
  throw std::runtime_error("sampleBeta_ST_sbayesrc: mixture_var[0] must be 0 for null component.");
 }

 const double vei_safe = std::max(vei_t, 1e-300);
 const double score = r(iu) + wi * b(iu);

 std::vector<double> logw(static_cast<std::size_t>(Kmix), 0.0);

 for (int k = 0; k < Kmix; ++k) {
  const double pik = std::max(pi_class_row(static_cast<arma::uword>(k)), 1e-300);

  if (k == 0) {
   logw[static_cast<std::size_t>(k)] = std::log(pik);
  } else {
   const double mult = mixture_var(static_cast<arma::uword>(k));

   if (!std::isfinite(mult) || mult <= 0.0) {
    throw std::runtime_error("sampleBeta_ST_sbayesrc: non-null mixture variance must be positive.");
   }

   const double vbi = std::max(vb_t * mult, 1e-12);
   const double denom = std::max(vei_safe + wi * vbi, 1e-300);

   const double logBF =
    0.5 * std::log(vei_safe / denom) +
    0.5 * score * score * vbi / (vei_safe * denom);

   logw[static_cast<std::size_t>(k)] = std::log(pik) + logBF;
  }
 }

 const int g_new = sample_categorical_from_logweights(logw, gen);
 double b_new = 0.0;

 if (g_new > 0) {
  const double mult = mixture_var(static_cast<arma::uword>(g_new));
  const double vbi = std::max(vb_t * mult, 1e-12);
  const double lhs = wi + vei_safe / vbi;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);

  std::normal_distribution<double> norm01(0.0, 1.0);
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
 gamma(iu) = g_new;
}

// -----------------------------------------------------------------------------
// Main exported function
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_sbayesrc_annot(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> gamma_init,
  bool use_gamma_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  std::string ld_prefix,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<int> annotation_class,
  int n_classes,
  arma::vec mixture_var,
  arma::mat pi_class_init,
  arma::mat alpha_class,
  bool updatePiClass,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
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

 if (n_classes <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: n_classes must be positive.");
 }

 validate_class_index(annotation_class, m, n_classes);

 const int C = n_classes;
 const int Kmix = static_cast<int>(mixture_var.n_elem);

 if (Kmix < 2) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: mixture_var must contain null + at least one non-null component.");
 }

 if (!std::isfinite(mixture_var(0)) || mixture_var(0) != 0.0) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: mixture_var[0] must be exactly 0.");
 }

 for (int k = 1; k < Kmix; ++k) {
  if (!std::isfinite(mixture_var(static_cast<arma::uword>(k))) ||
      mixture_var(static_cast<arma::uword>(k)) <= 0.0) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: non-null mixture_var values must be positive.");
  }
 }

 if ((int)pi_class_init.n_rows != C || (int)pi_class_init.n_cols != Kmix) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: pi_class_init must be n_classes x length(mixture_var).");
 }

 if ((int)alpha_class.n_rows != C || (int)alpha_class.n_cols != Kmix) {
  throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: alpha_class must be n_classes x length(mixture_var).");
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

 if (use_gamma_init) {
  if (static_cast<int>(gamma_init.size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: gamma_init must have length nt when use_gamma_init = true.");
  }

  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(gamma_init[t].size()) != m) {
    throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: each gamma_init[t] must have length m.");
   }
  }
 }

 // --------------------------------------------------------------------------
 // Convert inputs
 // --------------------------------------------------------------------------

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);
 arma::Mat<int> gamma_mat(nt, m, arma::fill::zeros);
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
   throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: priors must be nt x nt.");
  }

  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
   sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
  }
 }

 // --------------------------------------------------------------------------
 // Validate shared scaling and build LD
 // --------------------------------------------------------------------------

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

 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi_nonnull(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);

 std::vector<arma::mat> pi_class_mean(static_cast<std::size_t>(nt));
 std::vector<arma::mat> pi_class_final(static_cast<std::size_t>(nt));

 for (int t = 0; t < nt; ++t) {
  pi_class_mean[static_cast<std::size_t>(t)] = arma::mat(C, Kmix, arma::fill::zeros);
  pi_class_final[static_cast<std::size_t>(t)] = arma::mat(C, Kmix, arma::fill::zeros);
 }

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
 << "STBLR SBayesRC-like CSR OpenMP requested threads = "
 << nthreads
 << ", omp_get_max_threads = "
 << omp_get_max_threads()
 << ", num procs = "
 << omp_get_num_procs()
 << "\n";
#endif

 Rcpp::Rcout
 << "STBLR SBayesRC-like CSR: C=" << C
 << ", Kmix=" << Kmix
 << ", updatePiClass=" << updatePiClass
 << ", mixture_var=";

 for (int k = 0; k < Kmix; ++k) {
  Rcpp::Rcout << mixture_var(static_cast<arma::uword>(k));
  if (k + 1 < Kmix) Rcpp::Rcout << ",";
 }
 Rcpp::Rcout << "\n";

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
   arma::Row<int> gamma_t(m, arma::fill::zeros);
   arma::Row<int> d_t(m, arma::fill::zeros);

   if (use_gamma_init) {
    for (int i = 0; i < m; ++i) {
     const int g = static_cast<int>(gamma_init[t][i]);

     if (g < 0 || g >= Kmix) {
      throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: invalid gamma_init value.");
     }

     gamma_t(static_cast<arma::uword>(i)) = g;
     d_t(static_cast<arma::uword>(i)) = g > 0 ? 1 : 0;
    }
   } else {
    for (int i = 0; i < m; ++i) {
     const arma::uword iu = static_cast<arma::uword>(i);
     gamma_t(iu) = b_t(iu) != 0.0 ? 1 : 0;
     d_t(iu) = gamma_t(iu) > 0 ? 1 : 0;
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

   if (!std::isfinite(vb_t) || vb_t <= 0.0) {
    throw std::runtime_error("initial vb_t must be positive.");
   }

   if (!std::isfinite(ve_t) || ve_t <= 0.0) {
    throw std::runtime_error("initial ve_t must be positive.");
   }

   arma::mat pi_class_t = normalize_rows_to_simplex(pi_class_init);

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);

   arma::mat pi_class_accum(C, Kmix, arma::fill::zeros);

   double nsamples_t = 0.0;

   for (int it = 0; it < nit + nburn; ++it) {
    for (int isort = 0; isort < m; ++isort) {
     const int i = order[static_cast<std::size_t>(isort)];
     const int c = annotation_class[static_cast<std::size_t>(i)];

     sampleBeta_ST_sbayesrc(
      i,
      c,
      pi_class_t.row(static_cast<arma::uword>(c)),
      mixture_var,
      vb_t,
      vei_t,
      ww_t,
      r_t,
      b_t,
      gamma_t,
      ld,
      gen_t
     );
    }

    for (int i = 0; i < m; ++i) {
     d_t(static_cast<arma::uword>(i)) = gamma_t(static_cast<arma::uword>(i)) > 0 ? 1 : 0;
    }

    if (updatePiClass) {
     arma::Mat<int> gamma_one(1, m, arma::fill::zeros);
     gamma_one.row(0) = gamma_t;

     update_pi_class_dirichlet(
      pi_class_t,
      gamma_one,
      annotation_class,
      alpha_class,
      C,
      Kmix,
      gen_t
     );
    }

    if (updateB) {
     sampleB_ST_sbayesrc(
      m,
      nub,
      vb_t,
      b_t,
      gamma_t,
      mixture_var,
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

    vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);

    if (!std::isfinite(vg_t)) {
     throw std::runtime_error(
       "vg became NaN/Inf after computeG. iter=" +
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

    double global_nonnull = 0.0;
    for (int c = 0; c < C; ++c) {
     double row_nonnull = 0.0;
     for (int k = 1; k < Kmix; ++k) {
      row_nonnull += pi_class_t(static_cast<arma::uword>(c), static_cast<arma::uword>(k));
     }
     global_nonnull += row_nonnull;
    }
    global_nonnull /= static_cast<double>(C);

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    pis_t(static_cast<arma::uword>(it)) = global_nonnull;

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;

     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += static_cast<double>(gamma_t(iu) > 0 ? 1 : 0);
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
   gamma_mat.row(static_cast<arma::uword>(t)) = gamma_t;
   d_mat.row(static_cast<arma::uword>(t))  = d_t;

   vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
   vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
   ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
   pis_mat.row(static_cast<arma::uword>(t)) = pis_t;

   final_vb(static_cast<arma::uword>(t)) = vb_t;
   final_ve(static_cast<arma::uword>(t)) = ve_t;
   final_vg(static_cast<arma::uword>(t)) = vg_t;
   final_pi_nonnull(static_cast<arma::uword>(t)) = pis_t(static_cast<arma::uword>(nit + nburn - 1));
   nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;

   pi_class_mean[static_cast<std::size_t>(t)] = pi_class_accum;
   pi_class_final[static_cast<std::size_t>(t)] = pi_class_t;

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

 // --------------------------------------------------------------------------
 // Build result
 // --------------------------------------------------------------------------

 std::vector<std::vector<std::vector<double>>> result(20);

 for (int k = 0; k < 20; ++k) {
  result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
 }

 const int flat_pi_len = C * Kmix;

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

  result[18][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(flat_pi_len));
  result[19][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(flat_pi_len));
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

  result[16][ts][0] = 1.0 - final_pi_nonnull(static_cast<arma::uword>(t));
  result[16][ts][1] = final_pi_nonnull(static_cast<arma::uword>(t));

  double mean_pi = 0.0;
  int npi = 0;

  for (int it = nburn; it < nit + nburn; ++it) {
   mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
   ++npi;
  }

  if (npi > 0) mean_pi /= static_cast<double>(npi);
  else mean_pi = final_pi_nonnull(static_cast<arma::uword>(t));

  result[17][ts][0] = 1.0 - mean_pi;
  result[17][ts][1] = mean_pi;

  int pos = 0;
  for (int c = 0; c < C; ++c) {
   for (int k = 0; k < Kmix; ++k) {
    result[18][ts][static_cast<std::size_t>(pos)] =
     pi_class_mean[ts](static_cast<arma::uword>(c), static_cast<arma::uword>(k));

     result[19][ts][static_cast<std::size_t>(pos)] =
      pi_class_final[ts](static_cast<arma::uword>(c), static_cast<arma::uword>(k));

      ++pos;
   }
  }
 }

 return result;
}

