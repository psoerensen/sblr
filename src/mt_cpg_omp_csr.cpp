// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"

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

// -----------------------------------------------------------------------------
// Shared flat LD structure
// Stores pre-scaled X_i'X_j, not raw LD correlation.
// Disk input is expected to contain raw LD correlations r_ij in 0-based CSR.
// The builder symmetrizes the LD object for direct marker-wise Gibbs updates.
// -----------------------------------------------------------------------------

struct LDCSR {
 std::vector<uint64_t> ptr;  // length m + 1
 std::vector<int> idx;       // neighbor marker index
 std::vector<float> xij;     // pre-scaled X_i'X_j
};

inline void read_exact_file_mt(
  const std::string& path,
  void* data,
  std::size_t nbytes
) {
 FILE* fs = std::fopen(path.c_str(), "rb");

 if (!fs) {
  throw std::runtime_error("Could not open file: " + path);
 }

 const std::size_t got = std::fread(data, 1, nbytes, fs);
 std::fclose(fs);

 if (got != nbytes) {
  throw std::runtime_error("Short read from file: " + path);
 }
}

inline uint64_t parse_uint64_from_meta_mt(
  const std::string& value,
  const std::string& key
) {
 if (value.empty()) {
  throw std::runtime_error("Empty metadata value for key: " + key);
 }

 char* endptr = nullptr;
 const unsigned long long out = std::strtoull(value.c_str(), &endptr, 10);

 if (endptr == value.c_str() || *endptr != '\0') {
  throw std::runtime_error("Invalid unsigned integer metadata value for key: " + key);
 }

 return static_cast<uint64_t>(out);
}

inline LDCSR read_and_build_mt_ld_csr(
  const std::string& prefix,
  int m,
  const std::vector<double>& xx
) {
 const std::string row_file  = prefix + ".row_ptr.u64.bin";
 const std::string col_file  = prefix + ".col_idx.u32.0based.bin";
 const std::string val_file  = prefix + ".values.f32.bin";
 const std::string meta_file = prefix + ".meta.txt";

 if (m <= 0) {
  throw std::runtime_error("read_and_build_mt_ld_csr: m must be positive.");
 }

 if (static_cast<int>(xx.size()) != m) {
  throw std::runtime_error("read_and_build_mt_ld_csr: xx must have length m.");
 }

 std::ifstream meta(meta_file.c_str());
 if (!meta.is_open()) {
  throw std::runtime_error("Could not open metadata file: " + meta_file);
 }

 int m_meta = -1;
 uint64_t nnz_u64 = 0;
 bool have_nnz = false;

 std::string line;
 while (std::getline(meta, line)) {
  const std::string key_m   = "n_variants=";
  const std::string key_nnz = "nnz=";

  if (line.rfind(key_m, 0) == 0) {
   m_meta = std::stoi(line.substr(key_m.size()));
  } else if (line.rfind(key_nnz, 0) == 0) {
   nnz_u64 = parse_uint64_from_meta_mt(line.substr(key_nnz.size()), "nnz");
   have_nnz = true;
  }
 }
 meta.close();

 if (m_meta <= 0) {
  throw std::runtime_error("Could not read n_variants from metadata.");
 }

 if (m_meta != m) {
  throw std::runtime_error("LD metadata n_variants does not match marker dimension.");
 }

 if (!have_nnz) {
  throw std::runtime_error("Could not read nnz from metadata.");
 }

 const std::size_t nnz = static_cast<std::size_t>(nnz_u64);

 std::vector<uint64_t> row_ptr(static_cast<std::size_t>(m) + 1);
 std::vector<uint32_t> col_idx_u32(nnz);
 std::vector<float> values_r(nnz);

 read_exact_file_mt(
  row_file,
  row_ptr.data(),
  row_ptr.size() * sizeof(uint64_t)
 );

 read_exact_file_mt(
  col_file,
  col_idx_u32.data(),
  col_idx_u32.size() * sizeof(uint32_t)
 );

 read_exact_file_mt(
  val_file,
  values_r.data(),
  values_r.size() * sizeof(float)
 );

 if (row_ptr[0] != 0 || row_ptr[static_cast<std::size_t>(m)] != nnz_u64) {
  throw std::runtime_error("Invalid LD row_ptr: expected 0-based row_ptr ending at nnz.");
 }

 for (int i = 0; i < m; ++i) {
  if (row_ptr[static_cast<std::size_t>(i + 1)] < row_ptr[static_cast<std::size_t>(i)]) {
   throw std::runtime_error("Invalid LD row_ptr: row pointers are not nondecreasing.");
  }

  if (!std::isfinite(xx[static_cast<std::size_t>(i)]) || xx[static_cast<std::size_t>(i)] <= 0.0) {
   throw std::runtime_error(
     "read_and_build_mt_ld_csr: xx contains invalid value at marker " +
      std::to_string(i)
   );
  }
 }

 // First pass: count symmetric degrees.
 std::vector<uint64_t> degree(static_cast<std::size_t>(m), 0);

 for (int i = 0; i < m; ++i) {
  const uint64_t start = row_ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = row_ptr[static_cast<std::size_t>(i + 1)];

  if (end > nnz_u64) {
   throw std::runtime_error("Invalid LD row_ptr: row end exceeds nnz.");
  }

  for (uint64_t p = start; p < end; ++p) {
   const uint32_t j_u32 = col_idx_u32[static_cast<std::size_t>(p)];

   if (j_u32 >= static_cast<uint32_t>(m)) {
    throw std::runtime_error("LD column index out of range.");
   }

   const int j = static_cast<int>(j_u32);

   if (j == i) continue;

   ++degree[static_cast<std::size_t>(i)];
   ++degree[static_cast<std::size_t>(j)];
  }
 }

 LDCSR ld;
 ld.ptr.resize(static_cast<std::size_t>(m) + 1);
 ld.ptr[0] = 0;

 for (int i = 0; i < m; ++i) {
  ld.ptr[static_cast<std::size_t>(i + 1)] =
   ld.ptr[static_cast<std::size_t>(i)] + degree[static_cast<std::size_t>(i)];
 }

 const uint64_t nnz_sym = ld.ptr[static_cast<std::size_t>(m)];

 ld.idx.resize(static_cast<std::size_t>(nnz_sym));
 ld.xij.resize(static_cast<std::size_t>(nnz_sym));

 std::vector<uint64_t> offset = ld.ptr;

 double max_abs_rij = 0.0;
 double max_abs_xij = 0.0;

 // Second pass: fill symmetric flat CSR.
 for (int i = 0; i < m; ++i) {
  const uint64_t start = row_ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = row_ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = static_cast<int>(col_idx_u32[static_cast<std::size_t>(p)]);

   if (j == i) continue;

   const double rij = static_cast<double>(values_r[static_cast<std::size_t>(p)]);

   if (!std::isfinite(rij)) {
    throw std::runtime_error("LD value contains NaN/Inf.");
   }

   max_abs_rij = std::max(max_abs_rij, std::abs(rij));

   if (std::abs(rij) > 1.0001) {
    throw std::runtime_error(
      "LD value is not a correlation. Did you pass X_i'X_j instead of r_ij?"
    );
   }

   const double xij =
    rij * std::sqrt(xx[static_cast<std::size_t>(i)] * xx[static_cast<std::size_t>(j)]);

   if (!std::isfinite(xij)) {
    throw std::runtime_error("Computed X_i'X_j contains NaN/Inf.");
   }

   max_abs_xij = std::max(max_abs_xij, std::abs(xij));

   const float xij_f = static_cast<float>(xij);

   const uint64_t pos_i = offset[static_cast<std::size_t>(i)]++;
   ld.idx[static_cast<std::size_t>(pos_i)] = j;
   ld.xij[static_cast<std::size_t>(pos_i)] = xij_f;

   const uint64_t pos_j = offset[static_cast<std::size_t>(j)]++;
   ld.idx[static_cast<std::size_t>(pos_j)] = i;
   ld.xij[static_cast<std::size_t>(pos_j)] = xij_f;
  }
 }

 for (int i = 0; i < m; ++i) {
  if (offset[static_cast<std::size_t>(i)] != ld.ptr[static_cast<std::size_t>(i + 1)]) {
   throw std::runtime_error("Internal LD CSR fill-count mismatch.");
  }
 }

 Rcpp::Rcout
 << "MT flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
 << ", symmetric nnz=" << static_cast<double>(nnz_sym)
 << ", max_abs_rij=" << max_abs_rij
 << ", max_abs_xij=" << max_abs_xij
 << "\n";

 return ld;
}

// -----------------------------------------------------------------------------
// Utility
// -----------------------------------------------------------------------------

inline void make_spd_symmetric(
  arma::mat& A,
  double ridge = 1e-10
) {
 A = 0.5 * (A + A.t());

 arma::vec eigval;
 arma::mat eigvec;

 if (!arma::eig_sym(eigval, eigvec, A)) {
  A = arma::diagmat(arma::clamp(A.diag(), ridge, arma::datum::inf));
  return;
 }

 eigval = arma::clamp(eigval, ridge, arma::datum::inf);
 A = eigvec * arma::diagmat(eigval) * eigvec.t();
 A = 0.5 * (A + A.t());
}

// -----------------------------------------------------------------------------
// Variance component updates
// -----------------------------------------------------------------------------

inline void sampleB_cpg_arma_omp_csr(
  int nt,
  int nub,
  arma::mat& B,
  const arma::mat& beta,
  const arma::Mat<int>& d,
  const arma::mat& ssb_prior,
  std::mt19937& gen
) {
 if (nub <= nt - 1) {
  throw std::runtime_error("sampleB_cpg_arma_omp_csr: nub must be > nt - 1.");
 }

 if (!beta.is_finite()) {
  throw std::runtime_error("sampleB_cpg_arma_omp_csr: beta contains NaN/Inf.");
 }

 arma::mat S_post = static_cast<double>(nub) * ssb_prior;
 int n_active = 0;

 for (arma::uword i = 0; i < beta.n_cols; ++i) {
  bool active = false;

  for (int t = 0; t < nt; ++t) {
   if (d(static_cast<arma::uword>(t), i) > 0) {
    active = true;
    break;
   }
  }

  if (active) {
   const arma::vec bi = beta.col(i);
   S_post += bi * bi.t();
   ++n_active;
  }
 }

 S_post = 0.5 * (S_post + S_post.t());

 if (!S_post.is_finite()) {
  throw std::runtime_error("sampleB_cpg_arma_omp_csr: S_post contains NaN/Inf.");
 }

 arma::mat L;
 bool ok = arma::chol(L, S_post, "lower");

 double jitter = 1e-10;
 int tries = 0;

 while (!ok && tries < 8) {
  S_post.diag() += jitter;
  S_post = 0.5 * (S_post + S_post.t());
  ok = arma::chol(L, S_post, "lower");
  jitter *= 10.0;
  ++tries;
 }

 if (!ok) {
  throw std::runtime_error("sampleB_cpg_arma_omp_csr: S_post is not SPD after jitter.");
 }

 const unsigned int df_post =
  static_cast<unsigned int>(nub + std::max(n_active, 1));

 B = rinvwishart_cpp(df_post, S_post, gen);

 B = 0.5 * (B + B.t());
 B.diag() = arma::clamp(B.diag(), 1e-12, arma::datum::inf);

 if (!B.is_finite()) {
  throw std::runtime_error("sampleB_cpg_arma_omp_csr: sampled B contains NaN/Inf.");
 }
}

inline void sampleB_active_shrunk_csr(
  int nt,
  int m,
  int nub,
  arma::mat& B,
  const arma::Mat<int>& d,
  const arma::mat& b,
  const arma::mat& ssb_prior,
  std::mt19937& gen,
  double corr_df = 20.0,
  double ridge = 1e-8
) {
 arma::vec var_draw(nt, arma::fill::zeros);
 arma::mat R(nt, nt, arma::fill::eye);

 for (int t = 0; t < nt; ++t) {
  double ss = 0.0;
  int n_active = 0;

  for (int i = 0; i < m; ++i) {
   if (d(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) > 0) {
    ss += b(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) *
     b(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
    ++n_active;
   }
  }

  const double scale = ss + static_cast<double>(nub) * ssb_prior(t, t);
  const double df = static_cast<double>(nub + n_active);

  if (!std::isfinite(scale) || scale <= 0.0) {
   throw std::runtime_error("sampleB_active_shrunk_csr: invalid variance scale.");
  }

  std::chi_squared_distribution<double> rchisq(std::max(df, 1.0));
  const double x = std::max(rchisq(gen), 1e-300);

  var_draw(t) = std::max(scale / x, ridge);
 }

 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = t1 + 1; t2 < nt; ++t2) {
   double s12 = 0.0;
   int n_shared = 0;

   for (int i = 0; i < m; ++i) {
    if (d(static_cast<arma::uword>(t1), static_cast<arma::uword>(i)) > 0 &&
        d(static_cast<arma::uword>(t2), static_cast<arma::uword>(i)) > 0) {
     s12 += b(static_cast<arma::uword>(t1), static_cast<arma::uword>(i)) *
      b(static_cast<arma::uword>(t2), static_cast<arma::uword>(i));
     ++n_shared;
    }
   }

   double r = 0.0;

   if (n_shared > 0) {
    const double cov_hat = s12 / static_cast<double>(n_shared);
    const double denom = std::sqrt(var_draw(t1) * var_draw(t2));

    if (denom > 0.0 && std::isfinite(denom)) {
     r = cov_hat / denom;
    }

    const double corr_shrink =
     static_cast<double>(n_shared) /
      (static_cast<double>(n_shared) + corr_df);

    r *= corr_shrink;
   }

   r = std::max(-0.95, std::min(0.95, r));

   R(t1, t2) = r;
   R(t2, t1) = r;
  }
 }

 arma::vec sd = arma::sqrt(var_draw);
 arma::mat D = arma::diagmat(sd);

 B = D * R * D;
 B.diag() += ridge;
 B = 0.5 * (B + B.t());

 arma::vec eigval;
 arma::mat eigvec;

 if (!arma::eig_sym(eigval, eigvec, B)) {
  B = arma::diagmat(arma::clamp(B.diag(), ridge, arma::datum::inf));
  return;
 }

 eigval = arma::clamp(eigval, ridge, arma::datum::inf);
 B = eigvec * arma::diagmat(eigval) * eigvec.t();
 B = 0.5 * (B + B.t());
}

inline void sampleE_cpg_arma_omp_csr(
  int nt,
  int nue,
  arma::mat& E,
  const arma::mat& b,
  const arma::mat& wy,
  const arma::mat& r,
  const arma::mat& sse_prior,
  const arma::mat& YY,
  const std::vector<int>& n,
  std::mt19937& gen,
  bool update_full_E = true,
  double ridge = 1e-10
) {
 if ((int)n.size() != nt) {
  throw std::runtime_error("sampleE_cpg_arma_omp_csr: n must have length nt.");
 }

 if ((int)YY.n_rows != nt || (int)YY.n_cols != nt) {
  throw std::runtime_error("sampleE_cpg_arma_omp_csr: YY must be nt x nt.");
 }

 if ((int)b.n_rows != nt || (int)wy.n_rows != nt || (int)r.n_rows != nt) {
  throw std::runtime_error("sampleE_cpg_arma_omp_csr: trait dimension mismatch.");
 }

 if (b.n_cols != wy.n_cols || b.n_cols != r.n_cols) {
  throw std::runtime_error("sampleE_cpg_arma_omp_csr: marker dimension mismatch.");
 }

 if (!update_full_E) {
  E.zeros(nt, nt);

  for (int t = 0; t < nt; ++t) {
   const double fitted_quad = arma::dot(b.row(t), r.row(t) + wy.row(t));
   const double sse = YY(t, t) - fitted_quad;
   const double scale = sse + static_cast<double>(nue) * sse_prior(t, t);

   if (!std::isfinite(scale) || scale <= 0.0) {
    throw std::runtime_error("sampleE_cpg_arma_omp_csr: invalid diagonal residual scale.");
   }

   std::chi_squared_distribution<double> rchisq(n[t] + nue);
   const double chi2 = std::max(rchisq(gen), 1e-300);

   E(t, t) = std::max(scale / chi2, ridge);
  }

  return;
 }

 arma::mat Se(nt, nt, arma::fill::zeros);

 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = t1; t2 < nt; ++t2) {

   const double term12 =
    YY(t1, t2)
   - arma::dot(b.row(t1), wy.row(t2))
   - arma::dot(b.row(t2), wy.row(t1))
   + arma::dot(b.row(t1), wy.row(t2) - r.row(t2));

   double se = term12;

   if (t1 != t2) {
    const double term21 =
     YY(t2, t1)
    - arma::dot(b.row(t2), wy.row(t1))
    - arma::dot(b.row(t1), wy.row(t2))
    + arma::dot(b.row(t2), wy.row(t1) - r.row(t1));

    se = 0.5 * (term12 + term21);
   }

   Se(t1, t2) = se;
   Se(t2, t1) = se;
  }
 }

 Se = 0.5 * (Se + Se.t());

 for (int t = 0; t < nt; ++t) {
  if (!std::isfinite(Se(t, t)) || Se(t, t) <= 0.0) {
   throw std::runtime_error(
     "sampleE_cpg_arma_omp_csr: full residual cross-product has invalid diagonal."
   );
  }
 }

 arma::mat S_post = Se + static_cast<double>(nue) * sse_prior;
 S_post = 0.5 * (S_post + S_post.t());

 if (!S_post.is_finite()) {
  throw std::runtime_error("sampleE_cpg_arma_omp_csr: S_post contains NaN/Inf.");
 }

 arma::vec eigval;
 arma::mat eigvec;

 if (!arma::eig_sym(eigval, eigvec, S_post)) {
  throw std::runtime_error("sampleE_cpg_arma_omp_csr: eig_sym failed for S_post.");
 }

 const double min_eig = arma::min(eigval);

 if (!std::isfinite(min_eig) || min_eig <= -1e-6) {
  throw std::runtime_error(
    "sampleE_cpg_arma_omp_csr: full S_post is not positive semidefinite."
  );
 }

 eigval = arma::clamp(eigval, ridge, arma::datum::inf);
 S_post = eigvec * arma::diagmat(eigval) * eigvec.t();
 S_post = 0.5 * (S_post + S_post.t());

 const unsigned int df_post = static_cast<unsigned int>(n[0] + nue);

 E = rinvwishart_cpp(df_post, S_post, gen);

 E = 0.5 * (E + E.t());
 E.diag() = arma::clamp(E.diag(), ridge, arma::datum::inf);

 make_spd_symmetric(E, ridge);

 if (!E.is_finite()) {
  throw std::runtime_error("sampleE_cpg_arma_omp_csr: sampled E contains NaN/Inf.");
 }
}

inline void computeG_cpg_arma_omp_csr(
  int nt,
  const arma::mat& b,
  const arma::mat& wy,
  const arma::mat& r,
  const std::vector<int>& n,
  arma::mat& G
) {
 G.set_size(nt, nt);

 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = t1; t2 < nt; ++t2) {
   const double denom = std::sqrt(static_cast<double>(n[t1]) * static_cast<double>(n[t2]));

   double s12 = arma::dot(b.row(t1), wy.row(t2) - r.row(t2));

   if (t1 != t2) {
    const double s21 = arma::dot(b.row(t2), wy.row(t1) - r.row(t1));
    s12 = 0.5 * (s12 + s21);
   }

   const double gij = s12 / denom;
   G(t1, t2) = gij;
   G(t2, t1) = gij;
  }
 }
}

// -----------------------------------------------------------------------------
// Residual update using shared flat LD
// -----------------------------------------------------------------------------

inline void rebuild_residual_from_b_ld(
  int nt,
  int m,
  const arma::mat& wy,
  const arma::mat& ww,
  const arma::mat& b,
  arma::mat& r,
  const LDCSR& ld
) {
 r = wy;

 for (int i = 0; i < m; ++i) {
  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (int t = 0; t < nt; ++t) {
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword iu = static_cast<arma::uword>(i);

   const double bi = b(tu, iu);
   if (bi == 0.0) continue;

   r(tu, iu) -= ww(tu, iu) * bi;

   for (uint64_t p = start; p < end; ++p) {
    const int j = ld.idx[static_cast<std::size_t>(p)];
    r(tu, static_cast<arma::uword>(j)) -=
     static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * bi;
   }
  }
 }
}

// -----------------------------------------------------------------------------
// Marker-update workspace and model precomputation
// -----------------------------------------------------------------------------

struct MTModelInfo {
 std::vector<int> active;
};

struct MTMarkerWorkspace {
 arma::vec score;
 arma::vec z_full;
 arma::vec beta_new;
 arma::vec b_new;
 arma::vec sqrtw;
 arma::vec z;
 arma::vec rhs;
 arma::vec y;
 arma::vec mean;
 arma::vec eps;
 arma::mat C;
 arma::mat P;
 arma::mat L;
 std::vector<double> loglik;
 std::vector<double> w;

 MTMarkerWorkspace(int nt, int nmodels)
  : score(nt, arma::fill::zeros),
    z_full(nt, arma::fill::zeros),
    beta_new(nt, arma::fill::zeros),
    b_new(nt, arma::fill::zeros),
    sqrtw(nt, arma::fill::zeros),
    z(nt, arma::fill::zeros),
    rhs(nt, arma::fill::zeros),
    y(nt, arma::fill::zeros),
    mean(nt, arma::fill::zeros),
    eps(nt, arma::fill::zeros),
    C(nt, nt, arma::fill::zeros),
    P(nt, nt, arma::fill::zeros),
    L(nt, nt, arma::fill::zeros),
    loglik(static_cast<std::size_t>(nmodels), -INFINITY),
    w(static_cast<std::size_t>(nmodels), 0.0) {}
};

inline void fill_model_matrices_and_rhs(
  int i,
  const std::vector<int>& active,
  const arma::mat& Ei,
  const arma::mat& Bi,
  const arma::mat& ww,
  const arma::vec& z_full,
  MTMarkerWorkspace& ws,
  int q
) {
 for (int a = 0; a < q; ++a) {
  const int ta = active[static_cast<std::size_t>(a)];
  ws.sqrtw(static_cast<arma::uword>(a)) = std::sqrt(ww(ta, i));
  ws.z(static_cast<arma::uword>(a)) = z_full(ta);
 }

 for (int a = 0; a < q; ++a) {
  const int ta = active[static_cast<std::size_t>(a)];
  const double sqrtw_a = ws.sqrtw(static_cast<arma::uword>(a));

  ws.rhs(static_cast<arma::uword>(a)) = 0.0;

  for (int b2 = 0; b2 < q; ++b2) {
   const int tb = active[static_cast<std::size_t>(b2)];
   const double p_ab =
    sqrtw_a * Ei(ta, tb) * ws.sqrtw(static_cast<arma::uword>(b2));

   ws.P(static_cast<arma::uword>(a), static_cast<arma::uword>(b2)) = p_ab;
   ws.C(static_cast<arma::uword>(a), static_cast<arma::uword>(b2)) = Bi(ta, tb) + p_ab;
   ws.rhs(static_cast<arma::uword>(a)) += p_ab * ws.z(static_cast<arma::uword>(b2));
  }
 }

 // Symmetrize only the active q x q block.
 for (int a = 0; a < q; ++a) {
  for (int b2 = a + 1; b2 < q; ++b2) {
   const double c = 0.5 * (
    ws.C(static_cast<arma::uword>(a), static_cast<arma::uword>(b2)) +
     ws.C(static_cast<arma::uword>(b2), static_cast<arma::uword>(a))
   );
   ws.C(static_cast<arma::uword>(a), static_cast<arma::uword>(b2)) = c;
   ws.C(static_cast<arma::uword>(b2), static_cast<arma::uword>(a)) = c;
  }
 }
}

inline void sampleBetaCPG_Mt_arma_ld(
  int i,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  const std::vector<MTModelInfo>& model_info,
  MTMarkerWorkspace& ws,
  std::vector<double>& cmodel_local,
  const std::vector<double>& pi,
  const arma::mat& Ei,
  const arma::mat& Bi,
  const arma::mat& ww,
  arma::mat& r,
  arma::mat& beta,
  arma::mat& b,
  arma::Mat<int>& d,
  const LDCSR& ld,
  std::mt19937& gen
) {
 std::fill(ws.loglik.begin(), ws.loglik.end(), -INFINITY);
 std::fill(ws.w.begin(), ws.w.end(), 0.0);
 ws.score.zeros();
 ws.z_full.zeros();
 ws.beta_new.zeros();
 ws.b_new.zeros();

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm(0.0, 1.0);

 if ((int)Ei.n_rows != nt || (int)Ei.n_cols != nt) {
  throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: Ei must be nt x nt.");
 }

 if ((int)Bi.n_rows != nt || (int)Bi.n_cols != nt) {
  throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: Bi must be nt x nt.");
 }

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  const arma::uword iu = static_cast<arma::uword>(i);

  if (!std::isfinite(ww(tu, iu)) || ww(tu, iu) <= 0.0) {
   throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: ww contains invalid value.");
  }

  ws.score(tu) = r(tu, iu) + ww(tu, iu) * b(tu, iu);
  ws.z_full(tu) = ws.score(tu) / ww(tu, iu);
 }

 // --------------------------------------------------------------------------
 // 1. Compute model log posterior weights in active subspace.
 // Uses precomputed active indices and reusable workspace.
 // --------------------------------------------------------------------------

 for (int k = 0; k < nmodels; ++k) {
  const std::size_t ks = static_cast<std::size_t>(k);

  if (pi[ks] <= 0.0 || !std::isfinite(pi[ks])) {
   ws.loglik[ks] = -INFINITY;
   continue;
  }

  const std::vector<int>& active = model_info[ks].active;
  const int q = static_cast<int>(active.size());

  if (q == 0) {
   ws.loglik[ks] = std::log(pi[ks]);
   continue;
  }

  fill_model_matrices_and_rhs(
   i,
   active,
   Ei,
   Bi,
   ww,
   ws.z_full,
   ws,
   q
  );

  arma::mat Cq = ws.C.submat(0, 0, q - 1, q - 1);
  arma::vec rhsq = ws.rhs.subvec(0, q - 1);

  arma::mat Lq;
  bool ok = arma::chol(Lq, Cq, "lower");

  if (!ok) {
   double jitter = 1e-10;
   for (int tries = 0; tries < 8 && !ok; ++tries) {
    Cq.diag() += jitter;
    Cq = 0.5 * (Cq + Cq.t());
    ok = arma::chol(Lq, Cq, "lower");
    jitter *= 10.0;
   }
  }

  if (!ok) {
   ws.loglik[ks] = -INFINITY;
   continue;
  }

  const double logdet = 2.0 * arma::sum(arma::log(Lq.diag()));

  arma::vec yq = arma::solve(arma::trimatl(Lq), rhsq);
  arma::vec meanq = arma::solve(arma::trimatu(Lq.t()), yq);

  const double quad = arma::dot(rhsq, meanq);

  ws.loglik[ks] = std::log(pi[ks]) - 0.5 * logdet + 0.5 * quad;
 }

 // --------------------------------------------------------------------------
 // 2. Sample model indicator.
 // --------------------------------------------------------------------------

 const double max_log = *std::max_element(ws.loglik.begin(), ws.loglik.end());

 double sumw = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  const std::size_t ks = static_cast<std::size_t>(k);
  ws.w[ks] = std::isfinite(ws.loglik[ks]) ? std::exp(ws.loglik[ks] - max_log) : 0.0;
  sumw += ws.w[ks];
 }

 if (!std::isfinite(sumw) || sumw <= 0.0) {
  throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: invalid model weights.");
 }

 const double u = runif(gen) * sumw;

 int mselect = nmodels - 1;
 double cum = 0.0;

 for (int k = 0; k < nmodels; ++k) {
  cum += ws.w[static_cast<std::size_t>(k)];
  if (u <= cum) {
   mselect = k;
   break;
  }
 }

 cmodel_local[static_cast<std::size_t>(mselect)] += 1.0;

 const auto& msel = models[static_cast<std::size_t>(mselect)];

 for (int t = 0; t < nt; ++t) {
  d(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) =
   msel[static_cast<std::size_t>(t)];
 }

 // --------------------------------------------------------------------------
 // 3. Sample selected model in active subspace.
 // --------------------------------------------------------------------------

 const std::vector<int>& active = model_info[static_cast<std::size_t>(mselect)].active;
 const int q = static_cast<int>(active.size());

 ws.beta_new.zeros();
 ws.b_new.zeros();

 if (q > 0) {
  fill_model_matrices_and_rhs(
   i,
   active,
   Ei,
   Bi,
   ww,
   ws.z_full,
   ws,
   q
  );

  arma::mat Cq = ws.C.submat(0, 0, q - 1, q - 1);
  arma::vec rhsq = ws.rhs.subvec(0, q - 1);

  arma::mat Lq;
  bool ok = arma::chol(Lq, Cq, "lower");

  if (!ok) {
   double jitter = 1e-10;
   for (int tries = 0; tries < 8 && !ok; ++tries) {
    Cq.diag() += jitter;
    Cq = 0.5 * (Cq + Cq.t());
    ok = arma::chol(Lq, Cq, "lower");
    jitter *= 10.0;
   }
  }

  if (!ok) {
   throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: selected model chol failed.");
  }

  arma::vec yq = arma::solve(arma::trimatl(Lq), rhsq);
  arma::vec meanq = arma::solve(arma::trimatu(Lq.t()), yq);

  arma::vec epsq(q);
  for (int a = 0; a < q; ++a) {
   epsq(static_cast<arma::uword>(a)) = norm(gen);
  }

  arma::vec draw = meanq + arma::solve(arma::trimatu(Lq.t()), epsq);

  for (int a = 0; a < q; ++a) {
   const int t = active[static_cast<std::size_t>(a)];
   ws.beta_new(static_cast<arma::uword>(t)) = draw(static_cast<arma::uword>(a));
   ws.b_new(static_cast<arma::uword>(t)) = draw(static_cast<arma::uword>(a));
  }
 }

 // --------------------------------------------------------------------------
 // 4. Residual update using effective effects.
 // --------------------------------------------------------------------------

 const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
 const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  const arma::uword iu = static_cast<arma::uword>(i);

  const double diff = ws.b_new(tu) - b(tu, iu);

  if (diff != 0.0) {
   r(tu, iu) -= ww(tu, iu) * diff;

   for (uint64_t p = start; p < end; ++p) {
    const int j = ld.idx[static_cast<std::size_t>(p)];
    r(tu, static_cast<arma::uword>(j)) -=
     static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
   }
  }
 }

 beta.col(static_cast<arma::uword>(i)) = ws.beta_new;
 b.col(static_cast<arma::uword>(i)) = ws.b_new;
}

// -----------------------------------------------------------------------------
// Main exported function
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> mtblr_cpg_omp_csr(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<std::vector<double>> yy,
  std::vector<std::vector<double>> b_init,
  std::string ld_prefix,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<std::vector<int>> models,
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
  int seed,
  int method
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) {
  throw std::runtime_error("mtblr_cpg_omp_csr: nt must be positive.");
 }

 const int m = static_cast<int>(wy[0].size());
 const int nmodels = static_cast<int>(models.size());

 if (m <= 0 || nmodels <= 0) {
  throw std::runtime_error("mtblr_cpg_omp_csr: invalid dimensions.");
 }

 if (nit <= 0) {
  throw std::runtime_error("mtblr_cpg_omp_csr: nit must be positive.");
 }

 if (nburn < 0) {
  throw std::runtime_error("mtblr_cpg_omp_csr: nburn must be non-negative.");
 }

 if (nthin <= 0) {
  throw std::runtime_error("mtblr_cpg_omp_csr: nthin must be > 0.");
 }

 if ((int)ww.size() != nt || (int)b_init.size() != nt || (int)yy.size() != nt ||
     (int)ssb_prior.size() != nt || (int)sse_prior.size() != nt ||
     (int)n.size() != nt) {
  throw std::runtime_error("mtblr_cpg_omp_csr: inconsistent trait dimensions.");
 }

 if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
  throw std::runtime_error("mtblr_cpg_omp_csr: B must be nt x nt.");
 }

 if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
  throw std::runtime_error("mtblr_cpg_omp_csr: E must be nt x nt.");
 }

 if ((int)pi.size() != nmodels) {
  throw std::runtime_error("mtblr_cpg_omp_csr: pi must have length nmodels.");
 }

 std::vector<MTModelInfo> model_info(static_cast<std::size_t>(nmodels));
 for (int k = 0; k < nmodels; ++k) {
  model_info[static_cast<std::size_t>(k)].active.reserve(static_cast<std::size_t>(nt));
  for (int t = 0; t < nt; ++t) {
   if (models[static_cast<std::size_t>(k)][static_cast<std::size_t>(t)] != 0) {
    model_info[static_cast<std::size_t>(k)].active.push_back(t);
   }
  }
 }

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m || (int)ww[t].size() != m || (int)b_init[t].size() != m) {
   throw std::runtime_error("mtblr_cpg_omp_csr: inconsistent marker dimensions.");
  }

  if ((int)yy[t].size() != nt || (int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("mtblr_cpg_omp_csr: yy and priors must be nt x nt.");
  }
 }

 for (int k = 0; k < nmodels; ++k) {
  if ((int)models[static_cast<std::size_t>(k)].size() != nt) {
   throw std::runtime_error("mtblr_cpg_omp_csr: each model must have length nt.");
  }
 }

 if (n[0] <= 1) {
  throw std::runtime_error("mtblr_cpg_omp_csr: n[0] must be > 1.");
 }

 for (int t = 1; t < nt; ++t) {
  if (n[t] != n[0]) {
   throw std::runtime_error(
     "mtblr_cpg_omp_csr: current shared-LD scaling assumes equal n across traits."
   );
  }
 }

 double nsamples = 0.0;

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat beta_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);

 arma::mat YY_mat(nt, nt, arma::fill::zeros);
 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  for (int t2 = 0; t2 < nt; ++t2) {
   YY_mat(t, t2) = yy[t][t2];
   ssb_prior_mat(t, t2) = ssb_prior[t][t2];
   sse_prior_mat(t, t2) = sse_prior[t][t2];
  }

  for (int i = 0; i < m; ++i) {
   wy_mat(t, i) = wy[t][i];
   ww_mat(t, i) = ww[t][i];
   b_mat(t, i)  = b_init[t][i];
  }
 }

 // --------------------------------------------------------------------------
 // Marker update order
 // --------------------------------------------------------------------------

 std::vector<double> x2(static_cast<std::size_t>(m), 0.0);
 std::vector<int> order(static_cast<std::size_t>(m));

 for (int i = 0; i < m; ++i) {
  double best = 0.0;

  for (int t = 0; t < nt; ++t) {
   if (ww_mat(t, i) > 0.0) {
    const double bhat = wy_mat(t, i) / ww_mat(t, i);
    best = std::max(best, bhat * bhat);
   }
  }

  x2[static_cast<std::size_t>(i)] = best;
  order[static_cast<std::size_t>(i)] = i;
 }

 std::sort(order.begin(), order.end(),
           [&](int a, int b) { return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)]; });

 // --------------------------------------------------------------------------
 // Shared flat LD object
 // --------------------------------------------------------------------------

 for (int t = 1; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   const double w0 = ww_mat(0, i);
   const double wt = ww_mat(t, i);
   const double tol = 1e-8 * std::max(1.0, std::abs(w0));

   if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
    throw std::runtime_error("mtblr_cpg_omp_csr: ww contains invalid value before LD pre-scaling.");
   }

   if (std::abs(w0 - wt) > tol) {
    throw std::runtime_error("mtblr_cpg_omp_csr: ww differs across traits; pre-scaled shared MT LD is invalid.");
   }
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
 for (int i = 0; i < m; ++i) {
  const double wi = ww_mat(0, i);
  if (!std::isfinite(wi) || wi <= 0.0) {
   throw std::runtime_error("mtblr_cpg_omp_csr: ww contains invalid value in trait 0.");
  }
  xx[static_cast<std::size_t>(i)] = wi;
 }

 LDCSR ld = read_and_build_mt_ld_csr(ld_prefix, m, xx);

 // Initial residual: r = wy - X'X b
 rebuild_residual_from_b_ld(
  nt,
  m,
  wy_mat,
  ww_mat,
  b_mat,
  r_mat,
  ld
 );

 arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
 arma::mat Bi, Ei, G(nt, nt, arma::fill::zeros);

 if (!arma::inv_sympd(Bi, B)) {
  throw std::runtime_error("Initial Bi inversion failed.");
 }

 if (!arma::inv_sympd(Ei, E)) {
  throw std::runtime_error("Initial Ei inversion failed.");
 }

 std::vector<double> cmodel(static_cast<std::size_t>(nmodels), 0.0);
 std::vector<double> pis(static_cast<std::size_t>(nmodels), 0.0);

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);

 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);

 arma::mat cvbm(nt, nt, arma::fill::zeros);
 arma::mat cvem(nt, nt, arma::fill::zeros);
 arma::mat cvgm(nt, nt, arma::fill::zeros);

 std::mt19937 gen(seed);
 MTMarkerWorkspace marker_ws(nt, nmodels);

 for (int it = 0; it < nit + nburn; ++it) {

  std::fill(cmodel.begin(), cmodel.end(), 1.0);

  arma::mat E_marker = E;

  for (int t = 0; t < nt; ++t) {
   E_marker(t, t) += adjE * std::max(G(t, t), 0.0);
  }

  arma::mat Ei_marker;
  if (!arma::inv_sympd(Ei_marker, E_marker)) {
   throw std::runtime_error("Adjusted marker E inversion failed.");
  }

  // ------------------------------------------------------------------------
  // Marker updates. For correctness, these are serial because each marker
  // update modifies the shared residual r_mat.
  // ------------------------------------------------------------------------
  if (method == 4) {
   for (int isort = 0; isort < m; ++isort) {
    const int i = order[static_cast<std::size_t>(isort)];

    sampleBetaCPG_Mt_arma_ld(
     i,
     nt,
     nmodels,
     models,
     model_info,
     marker_ws,
     cmodel,
     pi,
     Ei_marker,
     Bi,
     ww_mat,
     r_mat,
     beta_mat,
     b_mat,
     d_mat,
     ld,
     gen
    );
   }
  }

  if (updateB && method == 4) {
   sampleB_cpg_arma_omp_csr(
    nt,
    static_cast<int>(nub),
    B,
    beta_mat,
    d_mat,
    ssb_prior_mat,
    gen
   );

   if (!arma::inv_sympd(Bi, B)) {
    throw std::runtime_error("Bi inversion failed.");
   }
  }

  if (updateE) {
   sampleE_cpg_arma_omp_csr(
    nt,
    static_cast<int>(nue),
    E,
    b_mat,
    wy_mat,
    r_mat,
    sse_prior_mat,
    YY_mat,
    n,
    gen,
    true
   );

   if (!arma::inv_sympd(Ei, E)) {
    throw std::runtime_error("Ei inversion failed.");
   }
  }

  if (updatePi && method == 4) {
   samplePi_cpg(cmodel, pi, gen);

   if (it >= nburn) {
    for (int k = 0; k < nmodels; ++k) {
     pis[static_cast<std::size_t>(k)] += pi[static_cast<std::size_t>(k)];
    }
   }
  }

  computeG_cpg_arma_omp_csr(nt, b_mat, wy_mat, r_mat, n, G);

  for (int t = 0; t < nt; ++t) {
   vbs_mat(t, it) = B(t, t);
   ves_mat(t, it) = E(t, t);
   vgs_mat(t, it) = G(t, t);
  }

  if (it >= nburn) {
   cvbm += B;
   cvem += E;
   cvgm += G;
  }

  if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
   nsamples += 1.0;

   for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
     if (d_mat(t, i) > 0) {
      dm_mat(t, i) += 1.0;
      bm_mat(t, i) += b_mat(t, i);
     }
    }
   }
  }
 }

 // --------------------------------------------------------------------------
 // Build result with same 20-slot structure
 // --------------------------------------------------------------------------

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

  result[16][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nmodels));
  result[17][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nmodels));

  result[18][static_cast<std::size_t>(t)].resize(4);
  result[19][static_cast<std::size_t>(t)].resize(2);
 }

 const double denom = std::max(nsamples, 1.0);
 const double post_iter_denom = static_cast<double>(std::max(nit, 1));

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int i = 0; i < m; ++i) {
   const std::size_t is = static_cast<std::size_t>(i);

   result[0][ts][is] = bm_mat(t, i) / denom;
   result[1][ts][is] = dm_mat(t, i) / denom;
   result[2][ts][is] = wy_mat(t, i);
   result[3][ts][is] = r_mat(t, i);
   result[4][ts][is] = b_mat(t, i);
   result[5][ts][is] = static_cast<double>(d_mat(t, i));
   result[6][ts][is] = static_cast<double>(i);
  }
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int it = 0; it < nit + nburn; ++it) {
   const std::size_t its = static_cast<std::size_t>(it);

   result[7][ts][its] = vbs_mat(t, it);
   result[8][ts][its] = vgs_mat(t, it);
   result[9][ts][its] = ves_mat(t, it);
  }
 }

 for (int t1 = 0; t1 < nt; ++t1) {
  const std::size_t t1s = static_cast<std::size_t>(t1);

  for (int t2 = 0; t2 < nt; ++t2) {
   const std::size_t t2s = static_cast<std::size_t>(t2);

   result[10][t1s][t2s] = cvbm(t1, t2) / post_iter_denom;
   result[11][t1s][t2s] = cvgm(t1, t2) / post_iter_denom;
   result[12][t1s][t2s] = cvem(t1, t2) / post_iter_denom;
   result[13][t1s][t2s] = B(t1, t2);
   result[14][t1s][t2s] = G(t1, t2);
   result[15][t1s][t2s] = E(t1, t2);
  }
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int k = 0; k < nmodels; ++k) {
   const std::size_t ks = static_cast<std::size_t>(k);
   result[16][ts][ks] = pi[ks];
   result[17][ts][ks] = pis[ks] / post_iter_denom;
  }
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int i = 0; i < 4; ++i) {
   result[18][ts][static_cast<std::size_t>(i)] = 0.0;
  }

  for (int i = 0; i < 2; ++i) {
   result[19][ts][static_cast<std::size_t>(i)] = 0.0;
  }
 }

 return result;
}


// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
//
// #include "cpg_samplers.h"
// #include "distributions.h"
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
// // -----------------------------------------------------------------------------
// // Shared flat LD structure
// // Stores pre-scaled X_i'X_j, not raw LD correlation.
// // Disk input is expected to contain raw LD correlations r_ij in 0-based CSR.
// // The builder symmetrizes the LD object for direct marker-wise Gibbs updates.
// // -----------------------------------------------------------------------------
//
// struct LDCSR {
//  std::vector<uint64_t> ptr;  // length m + 1
//  std::vector<int> idx;       // neighbor marker index
//  std::vector<float> xij;     // pre-scaled X_i'X_j
// };
//
// inline void read_exact_file_mt(
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
// inline uint64_t parse_uint64_from_meta_mt(
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
// inline LDCSR read_and_build_mt_ld_csr(
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
//   throw std::runtime_error("read_and_build_mt_ld_csr: m must be positive.");
//  }
//
//  if (static_cast<int>(xx.size()) != m) {
//   throw std::runtime_error("read_and_build_mt_ld_csr: xx must have length m.");
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
//    nnz_u64 = parse_uint64_from_meta_mt(line.substr(key_nnz.size()), "nnz");
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
//  read_exact_file_mt(
//   row_file,
//   row_ptr.data(),
//   row_ptr.size() * sizeof(uint64_t)
//  );
//
//  read_exact_file_mt(
//   col_file,
//   col_idx_u32.data(),
//   col_idx_u32.size() * sizeof(uint32_t)
//  );
//
//  read_exact_file_mt(
//   val_file,
//   values_r.data(),
//   values_r.size() * sizeof(float)
//  );
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
//      "read_and_build_mt_ld_csr: xx contains invalid value at marker " +
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
//  LDCSR ld;
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
//  for (int i = 0; i < m; ++i) {
//   if (offset[static_cast<std::size_t>(i)] != ld.ptr[static_cast<std::size_t>(i + 1)]) {
//    throw std::runtime_error("Internal LD CSR fill-count mismatch.");
//   }
//  }
//
//  Rcpp::Rcout
//  << "MT flat LD builder: input nnz=" << static_cast<double>(nnz_u64)
//  << ", symmetric nnz=" << static_cast<double>(nnz_sym)
//  << ", max_abs_rij=" << max_abs_rij
//  << ", max_abs_xij=" << max_abs_xij
//  << "\n";
//
//  return ld;
// }
//
// // -----------------------------------------------------------------------------
// // Utility
// // -----------------------------------------------------------------------------
//
// inline void make_spd_symmetric(
//   arma::mat& A,
//   double ridge = 1e-10
// ) {
//  A = 0.5 * (A + A.t());
//
//  arma::vec eigval;
//  arma::mat eigvec;
//
//  if (!arma::eig_sym(eigval, eigvec, A)) {
//   A = arma::diagmat(arma::clamp(A.diag(), ridge, arma::datum::inf));
//   return;
//  }
//
//  eigval = arma::clamp(eigval, ridge, arma::datum::inf);
//  A = eigvec * arma::diagmat(eigval) * eigvec.t();
//  A = 0.5 * (A + A.t());
// }
//
// // -----------------------------------------------------------------------------
// // Variance component updates
// // -----------------------------------------------------------------------------
//
// inline void sampleB_cpg_arma_omp_csr(
//   int nt,
//   int nub,
//   arma::mat& B,
//   const arma::mat& beta,
//   const arma::Mat<int>& d,
//   const arma::mat& ssb_prior,
//   std::mt19937& gen
// ) {
//  if (nub <= nt - 1) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: nub must be > nt - 1.");
//  }
//
//  if (!beta.is_finite()) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: beta contains NaN/Inf.");
//  }
//
//  arma::mat S_post = static_cast<double>(nub) * ssb_prior;
//  int n_active = 0;
//
//  for (arma::uword i = 0; i < beta.n_cols; ++i) {
//   bool active = false;
//
//   for (int t = 0; t < nt; ++t) {
//    if (d(static_cast<arma::uword>(t), i) > 0) {
//     active = true;
//     break;
//    }
//   }
//
//   if (active) {
//    const arma::vec bi = beta.col(i);
//    S_post += bi * bi.t();
//    ++n_active;
//   }
//  }
//
//  S_post = 0.5 * (S_post + S_post.t());
//
//  if (!S_post.is_finite()) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: S_post contains NaN/Inf.");
//  }
//
//  arma::mat L;
//  bool ok = arma::chol(L, S_post, "lower");
//
//  double jitter = 1e-10;
//  int tries = 0;
//
//  while (!ok && tries < 8) {
//   S_post.diag() += jitter;
//   S_post = 0.5 * (S_post + S_post.t());
//   ok = arma::chol(L, S_post, "lower");
//   jitter *= 10.0;
//   ++tries;
//  }
//
//  if (!ok) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: S_post is not SPD after jitter.");
//  }
//
//  const unsigned int df_post =
//   static_cast<unsigned int>(nub + std::max(n_active, 1));
//
//  B = rinvwishart_cpp(df_post, S_post, gen);
//
//  B = 0.5 * (B + B.t());
//  B.diag() = arma::clamp(B.diag(), 1e-12, arma::datum::inf);
//
//  if (!B.is_finite()) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: sampled B contains NaN/Inf.");
//  }
// }
//
// inline void sampleB_active_shrunk_csr(
//   int nt,
//   int m,
//   int nub,
//   arma::mat& B,
//   const arma::Mat<int>& d,
//   const arma::mat& b,
//   const arma::mat& ssb_prior,
//   std::mt19937& gen,
//   double corr_df = 20.0,
//   double ridge = 1e-8
// ) {
//  arma::vec var_draw(nt, arma::fill::zeros);
//  arma::mat R(nt, nt, arma::fill::eye);
//
//  for (int t = 0; t < nt; ++t) {
//   double ss = 0.0;
//   int n_active = 0;
//
//   for (int i = 0; i < m; ++i) {
//    if (d(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) > 0) {
//     ss += b(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) *
//      b(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
//     ++n_active;
//    }
//   }
//
//   const double scale = ss + static_cast<double>(nub) * ssb_prior(t, t);
//   const double df = static_cast<double>(nub + n_active);
//
//   if (!std::isfinite(scale) || scale <= 0.0) {
//    throw std::runtime_error("sampleB_active_shrunk_csr: invalid variance scale.");
//   }
//
//   std::chi_squared_distribution<double> rchisq(std::max(df, 1.0));
//   const double x = std::max(rchisq(gen), 1e-300);
//
//   var_draw(t) = std::max(scale / x, ridge);
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = t1 + 1; t2 < nt; ++t2) {
//    double s12 = 0.0;
//    int n_shared = 0;
//
//    for (int i = 0; i < m; ++i) {
//     if (d(static_cast<arma::uword>(t1), static_cast<arma::uword>(i)) > 0 &&
//         d(static_cast<arma::uword>(t2), static_cast<arma::uword>(i)) > 0) {
//      s12 += b(static_cast<arma::uword>(t1), static_cast<arma::uword>(i)) *
//       b(static_cast<arma::uword>(t2), static_cast<arma::uword>(i));
//      ++n_shared;
//     }
//    }
//
//    double r = 0.0;
//
//    if (n_shared > 0) {
//     const double cov_hat = s12 / static_cast<double>(n_shared);
//     const double denom = std::sqrt(var_draw(t1) * var_draw(t2));
//
//     if (denom > 0.0 && std::isfinite(denom)) {
//      r = cov_hat / denom;
//     }
//
//     const double corr_shrink =
//      static_cast<double>(n_shared) /
//       (static_cast<double>(n_shared) + corr_df);
//
//     r *= corr_shrink;
//    }
//
//    r = std::max(-0.95, std::min(0.95, r));
//
//    R(t1, t2) = r;
//    R(t2, t1) = r;
//   }
//  }
//
//  arma::vec sd = arma::sqrt(var_draw);
//  arma::mat D = arma::diagmat(sd);
//
//  B = D * R * D;
//  B.diag() += ridge;
//  B = 0.5 * (B + B.t());
//
//  arma::vec eigval;
//  arma::mat eigvec;
//
//  if (!arma::eig_sym(eigval, eigvec, B)) {
//   B = arma::diagmat(arma::clamp(B.diag(), ridge, arma::datum::inf));
//   return;
//  }
//
//  eigval = arma::clamp(eigval, ridge, arma::datum::inf);
//  B = eigvec * arma::diagmat(eigval) * eigvec.t();
//  B = 0.5 * (B + B.t());
// }
//
// inline void sampleE_cpg_arma_omp_csr(
//   int nt,
//   int nue,
//   arma::mat& E,
//   const arma::mat& b,
//   const arma::mat& wy,
//   const arma::mat& r,
//   const arma::mat& sse_prior,
//   const arma::mat& YY,
//   const std::vector<int>& n,
//   std::mt19937& gen,
//   bool update_full_E = true,
//   double ridge = 1e-10
// ) {
//  if ((int)n.size() != nt) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: n must have length nt.");
//  }
//
//  if ((int)YY.n_rows != nt || (int)YY.n_cols != nt) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: YY must be nt x nt.");
//  }
//
//  if ((int)b.n_rows != nt || (int)wy.n_rows != nt || (int)r.n_rows != nt) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: trait dimension mismatch.");
//  }
//
//  if (b.n_cols != wy.n_cols || b.n_cols != r.n_cols) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: marker dimension mismatch.");
//  }
//
//  if (!update_full_E) {
//   E.zeros(nt, nt);
//
//   for (int t = 0; t < nt; ++t) {
//    const double fitted_quad = arma::dot(b.row(t), r.row(t) + wy.row(t));
//    const double sse = YY(t, t) - fitted_quad;
//    const double scale = sse + static_cast<double>(nue) * sse_prior(t, t);
//
//    if (!std::isfinite(scale) || scale <= 0.0) {
//     throw std::runtime_error("sampleE_cpg_arma_omp_csr: invalid diagonal residual scale.");
//    }
//
//    std::chi_squared_distribution<double> rchisq(n[t] + nue);
//    const double chi2 = std::max(rchisq(gen), 1e-300);
//
//    E(t, t) = std::max(scale / chi2, ridge);
//   }
//
//   return;
//  }
//
//  arma::mat Se(nt, nt, arma::fill::zeros);
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = t1; t2 < nt; ++t2) {
//
//    const double term12 =
//     YY(t1, t2)
//    - arma::dot(b.row(t1), wy.row(t2))
//    - arma::dot(b.row(t2), wy.row(t1))
//    + arma::dot(b.row(t1), wy.row(t2) - r.row(t2));
//
//    double se = term12;
//
//    if (t1 != t2) {
//     const double term21 =
//      YY(t2, t1)
//     - arma::dot(b.row(t2), wy.row(t1))
//     - arma::dot(b.row(t1), wy.row(t2))
//     + arma::dot(b.row(t2), wy.row(t1) - r.row(t1));
//
//     se = 0.5 * (term12 + term21);
//    }
//
//    Se(t1, t2) = se;
//    Se(t2, t1) = se;
//   }
//  }
//
//  Se = 0.5 * (Se + Se.t());
//
//  for (int t = 0; t < nt; ++t) {
//   if (!std::isfinite(Se(t, t)) || Se(t, t) <= 0.0) {
//    throw std::runtime_error(
//      "sampleE_cpg_arma_omp_csr: full residual cross-product has invalid diagonal."
//    );
//   }
//  }
//
//  arma::mat S_post = Se + static_cast<double>(nue) * sse_prior;
//  S_post = 0.5 * (S_post + S_post.t());
//
//  if (!S_post.is_finite()) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: S_post contains NaN/Inf.");
//  }
//
//  arma::vec eigval;
//  arma::mat eigvec;
//
//  if (!arma::eig_sym(eigval, eigvec, S_post)) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: eig_sym failed for S_post.");
//  }
//
//  const double min_eig = arma::min(eigval);
//
//  if (!std::isfinite(min_eig) || min_eig <= -1e-6) {
//   throw std::runtime_error(
//     "sampleE_cpg_arma_omp_csr: full S_post is not positive semidefinite."
//   );
//  }
//
//  eigval = arma::clamp(eigval, ridge, arma::datum::inf);
//  S_post = eigvec * arma::diagmat(eigval) * eigvec.t();
//  S_post = 0.5 * (S_post + S_post.t());
//
//  const unsigned int df_post = static_cast<unsigned int>(n[0] + nue);
//
//  E = rinvwishart_cpp(df_post, S_post, gen);
//
//  E = 0.5 * (E + E.t());
//  E.diag() = arma::clamp(E.diag(), ridge, arma::datum::inf);
//
//  make_spd_symmetric(E, ridge);
//
//  if (!E.is_finite()) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: sampled E contains NaN/Inf.");
//  }
// }
//
// inline void computeG_cpg_arma_omp_csr(
//   int nt,
//   const arma::mat& b,
//   const arma::mat& wy,
//   const arma::mat& r,
//   const std::vector<int>& n,
//   arma::mat& G
// ) {
//  G.set_size(nt, nt);
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = t1; t2 < nt; ++t2) {
//    const double denom = std::sqrt(static_cast<double>(n[t1]) * static_cast<double>(n[t2]));
//
//    double s12 = arma::dot(b.row(t1), wy.row(t2) - r.row(t2));
//
//    if (t1 != t2) {
//     const double s21 = arma::dot(b.row(t2), wy.row(t1) - r.row(t1));
//     s12 = 0.5 * (s12 + s21);
//    }
//
//    const double gij = s12 / denom;
//    G(t1, t2) = gij;
//    G(t2, t1) = gij;
//   }
//  }
// }
//
// // -----------------------------------------------------------------------------
// // Residual update using shared flat LD
// // -----------------------------------------------------------------------------
//
// inline void rebuild_residual_from_b_ld(
//   int nt,
//   int m,
//   const arma::mat& wy,
//   const arma::mat& ww,
//   const arma::mat& b,
//   arma::mat& r,
//   const LDCSR& ld
// ) {
//  r = wy;
//
//  for (int i = 0; i < m; ++i) {
//   const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
//   const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];
//
//   for (int t = 0; t < nt; ++t) {
//    const arma::uword tu = static_cast<arma::uword>(t);
//    const arma::uword iu = static_cast<arma::uword>(i);
//
//    const double bi = b(tu, iu);
//    if (bi == 0.0) continue;
//
//    r(tu, iu) -= ww(tu, iu) * bi;
//
//    for (uint64_t p = start; p < end; ++p) {
//     const int j = ld.idx[static_cast<std::size_t>(p)];
//     r(tu, static_cast<arma::uword>(j)) -=
//      static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * bi;
//    }
//   }
//  }
// }
//
// // -----------------------------------------------------------------------------
// // Marker update using shared flat LD
// // -----------------------------------------------------------------------------
//
// inline void sampleBetaCPG_Mt_arma_ld(
//   int i,
//   int nt,
//   int nmodels,
//   const std::vector<std::vector<int>>& models,
//   std::vector<double>& cmodel_local,
//   const std::vector<double>& pi,
//   const arma::mat& Ei,
//   const arma::mat& Bi,
//   const arma::mat& ww,
//   arma::mat& r,
//   arma::mat& beta,
//   arma::mat& b,
//   arma::Mat<int>& d,
//   const LDCSR& ld,
//   std::mt19937& gen
// ) {
//  arma::vec score(nt, arma::fill::zeros);
//  arma::vec z_full(nt, arma::fill::zeros);
//  arma::vec beta_new(nt, arma::fill::zeros);
//  arma::vec b_new(nt, arma::fill::zeros);
//
//  std::vector<double> loglik(static_cast<std::size_t>(nmodels), -INFINITY);
//  std::vector<double> w(static_cast<std::size_t>(nmodels), 0.0);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm(0.0, 1.0);
//
//  if ((int)Ei.n_rows != nt || (int)Ei.n_cols != nt) {
//   throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: Ei must be nt x nt.");
//  }
//
//  if ((int)Bi.n_rows != nt || (int)Bi.n_cols != nt) {
//   throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: Bi must be nt x nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const arma::uword tu = static_cast<arma::uword>(t);
//   const arma::uword iu = static_cast<arma::uword>(i);
//
//   if (!std::isfinite(ww(tu, iu)) || ww(tu, iu) <= 0.0) {
//    throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: ww contains invalid value.");
//   }
//
//   score(t) = r(tu, iu) + ww(tu, iu) * b(tu, iu);
//   z_full(t) = score(t) / ww(tu, iu);
//  }
//
//  // --------------------------------------------------------------------------
//  // 1. Compute model log posterior weights in active subspace
//  // --------------------------------------------------------------------------
//
//  for (int k = 0; k < nmodels; ++k) {
//   if (pi[static_cast<std::size_t>(k)] <= 0.0 || !std::isfinite(pi[static_cast<std::size_t>(k)])) {
//    loglik[static_cast<std::size_t>(k)] = -INFINITY;
//    continue;
//   }
//
//   const auto& mk = models[static_cast<std::size_t>(k)];
//
//   std::vector<int> active;
//   active.reserve(static_cast<std::size_t>(nt));
//
//   for (int t = 0; t < nt; ++t) {
//    if (mk[static_cast<std::size_t>(t)]) active.push_back(t);
//   }
//
//   const int q = static_cast<int>(active.size());
//
//   if (q == 0) {
//    loglik[static_cast<std::size_t>(k)] = std::log(pi[static_cast<std::size_t>(k)]);
//    continue;
//   }
//
//   arma::uvec idx(q);
//   arma::vec sqrtw(q, arma::fill::zeros);
//
//   for (int a = 0; a < q; ++a) {
//    const int t = active[static_cast<std::size_t>(a)];
//    idx(static_cast<arma::uword>(a)) = static_cast<arma::uword>(t);
//    sqrtw(static_cast<arma::uword>(a)) = std::sqrt(ww(t, i));
//   }
//
//   arma::mat Ei_sub = Ei.submat(idx, idx);
//   arma::mat C = Bi.submat(idx, idx);
//   arma::vec z = z_full.elem(idx);
//
//   arma::mat D = arma::diagmat(sqrtw);
//   arma::mat P = D * Ei_sub * D;
//
//   C += P;
//   arma::vec rhs = P * z;
//
//   C = 0.5 * (C + C.t());
//
//   arma::mat L;
//   bool ok = arma::chol(L, C, "lower");
//
//   if (!ok) {
//    double jitter = 1e-10;
//    for (int tries = 0; tries < 8 && !ok; ++tries) {
//     C.diag() += jitter;
//     C = 0.5 * (C + C.t());
//     ok = arma::chol(L, C, "lower");
//     jitter *= 10.0;
//    }
//   }
//
//   if (!ok) {
//    loglik[static_cast<std::size_t>(k)] = -INFINITY;
//    continue;
//   }
//
//   const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
//
//   arma::vec y = arma::solve(arma::trimatl(L), rhs);
//   arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
//
//   const double quad = arma::dot(rhs, mean);
//
//   loglik[static_cast<std::size_t>(k)] =
//    std::log(pi[static_cast<std::size_t>(k)]) - 0.5 * logdet + 0.5 * quad;
//  }
//
//  // --------------------------------------------------------------------------
//  // 2. Sample model indicator
//  // --------------------------------------------------------------------------
//
//  const double max_log = *std::max_element(loglik.begin(), loglik.end());
//
//  double sumw = 0.0;
//  for (int k = 0; k < nmodels; ++k) {
//   w[static_cast<std::size_t>(k)] =
//    std::isfinite(loglik[static_cast<std::size_t>(k)]) ?
//    std::exp(loglik[static_cast<std::size_t>(k)] - max_log) : 0.0;
//   sumw += w[static_cast<std::size_t>(k)];
//  }
//
//  if (!std::isfinite(sumw) || sumw <= 0.0) {
//   throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: invalid model weights.");
//  }
//
//  const double u = runif(gen) * sumw;
//
//  int mselect = nmodels - 1;
//  double cum = 0.0;
//
//  for (int k = 0; k < nmodels; ++k) {
//   cum += w[static_cast<std::size_t>(k)];
//   if (u <= cum) {
//    mselect = k;
//    break;
//   }
//  }
//
//  cmodel_local[static_cast<std::size_t>(mselect)] += 1.0;
//
//  const auto& msel = models[static_cast<std::size_t>(mselect)];
//
//  for (int t = 0; t < nt; ++t) {
//   d(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) =
//    msel[static_cast<std::size_t>(t)];
//  }
//
//  // --------------------------------------------------------------------------
//  // 3. Sample selected model in active subspace
//  // --------------------------------------------------------------------------
//
//  std::vector<int> active;
//  active.reserve(static_cast<std::size_t>(nt));
//
//  for (int t = 0; t < nt; ++t) {
//   if (msel[static_cast<std::size_t>(t)]) active.push_back(t);
//  }
//
//  const int q = static_cast<int>(active.size());
//
//  beta_new.zeros();
//  b_new.zeros();
//
//  if (q > 0) {
//   arma::uvec idx(q);
//   arma::vec sqrtw(q, arma::fill::zeros);
//
//   for (int a = 0; a < q; ++a) {
//    const int t = active[static_cast<std::size_t>(a)];
//    idx(static_cast<arma::uword>(a)) = static_cast<arma::uword>(t);
//    sqrtw(static_cast<arma::uword>(a)) = std::sqrt(ww(t, i));
//   }
//
//   arma::mat Ei_sub = Ei.submat(idx, idx);
//   arma::mat C = Bi.submat(idx, idx);
//   arma::vec z = z_full.elem(idx);
//
//   arma::mat D = arma::diagmat(sqrtw);
//   arma::mat P = D * Ei_sub * D;
//
//   C += P;
//   arma::vec rhs = P * z;
//
//   C = 0.5 * (C + C.t());
//
//   arma::mat L;
//   bool ok = arma::chol(L, C, "lower");
//
//   if (!ok) {
//    double jitter = 1e-10;
//    for (int tries = 0; tries < 8 && !ok; ++tries) {
//     C.diag() += jitter;
//     C = 0.5 * (C + C.t());
//     ok = arma::chol(L, C, "lower");
//     jitter *= 10.0;
//    }
//   }
//
//   if (!ok) {
//    throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: selected model chol failed.");
//   }
//
//   arma::vec y = arma::solve(arma::trimatl(L), rhs);
//   arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
//
//   arma::vec eps(q);
//   for (int a = 0; a < q; ++a) {
//    eps(static_cast<arma::uword>(a)) = norm(gen);
//   }
//
//   arma::vec draw = mean + arma::solve(arma::trimatu(L.t()), eps);
//
//   for (int a = 0; a < q; ++a) {
//    const int t = active[static_cast<std::size_t>(a)];
//    beta_new(t) = draw(static_cast<arma::uword>(a));
//    b_new(t) = draw(static_cast<arma::uword>(a));
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // 4. Residual update using effective effects
//  // --------------------------------------------------------------------------
//
//  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
//  const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];
//
//  for (int t = 0; t < nt; ++t) {
//   const arma::uword tu = static_cast<arma::uword>(t);
//   const arma::uword iu = static_cast<arma::uword>(i);
//
//   const double diff = b_new(t) - b(tu, iu);
//
//   if (diff != 0.0) {
//    r(tu, iu) -= ww(tu, iu) * diff;
//
//    for (uint64_t p = start; p < end; ++p) {
//     const int j = ld.idx[static_cast<std::size_t>(p)];
//     r(tu, static_cast<arma::uword>(j)) -=
//      static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
//    }
//   }
//  }
//
//  beta.col(static_cast<arma::uword>(i)) = beta_new;
//  b.col(static_cast<arma::uword>(i)) = b_new;
// }
//
// // -----------------------------------------------------------------------------
// // Main exported function
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> mtblr_cpg_omp_csr(
//   std::vector<std::vector<double>> wy,
//   std::vector<std::vector<double>> ww,
//   std::vector<std::vector<double>> yy,
//   std::vector<std::vector<double>> b_init,
//   std::string ld_prefix,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<std::vector<int>> models,
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
//   int seed,
//   int method
// ) {
//  const int nt = static_cast<int>(wy.size());
//
//  if (nt <= 0) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: nt must be positive.");
//  }
//
//  const int m = static_cast<int>(wy[0].size());
//  const int nmodels = static_cast<int>(models.size());
//
//  if (m <= 0 || nmodels <= 0) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: invalid dimensions.");
//  }
//
//  if (nit <= 0) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: nit must be positive.");
//  }
//
//  if (nburn < 0) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: nburn must be non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: nthin must be > 0.");
//  }
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt || (int)yy.size() != nt ||
//      (int)ssb_prior.size() != nt || (int)sse_prior.size() != nt ||
//      (int)n.size() != nt) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: inconsistent trait dimensions.");
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: E must be nt x nt.");
//  }
//
//  if ((int)pi.size() != nmodels) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: pi must have length nmodels.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m || (int)ww[t].size() != m || (int)b_init[t].size() != m) {
//    throw std::runtime_error("mtblr_cpg_omp_csr: inconsistent marker dimensions.");
//   }
//
//   if ((int)yy[t].size() != nt || (int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("mtblr_cpg_omp_csr: yy and priors must be nt x nt.");
//   }
//  }
//
//  for (int k = 0; k < nmodels; ++k) {
//   if ((int)models[static_cast<std::size_t>(k)].size() != nt) {
//    throw std::runtime_error("mtblr_cpg_omp_csr: each model must have length nt.");
//   }
//  }
//
//  if (n[0] <= 1) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: n[0] must be > 1.");
//  }
//
//  for (int t = 1; t < nt; ++t) {
//   if (n[t] != n[0]) {
//    throw std::runtime_error(
//      "mtblr_cpg_omp_csr: current shared-LD scaling assumes equal n across traits."
//    );
//   }
//  }
//
//  double nsamples = 0.0;
//
//  arma::mat wy_mat(nt, m, arma::fill::zeros);
//  arma::mat ww_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::mat beta_mat(nt, m, arma::fill::zeros);
//  arma::mat r_mat(nt, m, arma::fill::zeros);
//
//  arma::mat YY_mat(nt, nt, arma::fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    YY_mat(t, t2) = yy[t][t2];
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//
//   for (int i = 0; i < m; ++i) {
//    wy_mat(t, i) = wy[t][i];
//    ww_mat(t, i) = ww[t][i];
//    b_mat(t, i)  = b_init[t][i];
//   }
//  }
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
//    if (ww_mat(t, i) > 0.0) {
//     const double bhat = wy_mat(t, i) / ww_mat(t, i);
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
//  // Shared flat LD object
//  // --------------------------------------------------------------------------
//
//  for (int t = 1; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    const double w0 = ww_mat(0, i);
//    const double wt = ww_mat(t, i);
//    const double tol = 1e-8 * std::max(1.0, std::abs(w0));
//
//    if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
//     throw std::runtime_error("mtblr_cpg_omp_csr: ww contains invalid value before LD pre-scaling.");
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error("mtblr_cpg_omp_csr: ww differs across traits; pre-scaled shared MT LD is invalid.");
//    }
//   }
//  }
//
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//  for (int i = 0; i < m; ++i) {
//   const double wi = ww_mat(0, i);
//   if (!std::isfinite(wi) || wi <= 0.0) {
//    throw std::runtime_error("mtblr_cpg_omp_csr: ww contains invalid value in trait 0.");
//   }
//   xx[static_cast<std::size_t>(i)] = wi;
//  }
//
//  LDCSR ld = read_and_build_mt_ld_csr(ld_prefix, m, xx);
//
//  // Initial residual: r = wy - X'X b
//  rebuild_residual_from_b_ld(
//   nt,
//   m,
//   wy_mat,
//   ww_mat,
//   b_mat,
//   r_mat,
//   ld
//  );
//
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
//  arma::mat Bi, Ei, G(nt, nt, arma::fill::zeros);
//
//  if (!arma::inv_sympd(Bi, B)) {
//   throw std::runtime_error("Initial Bi inversion failed.");
//  }
//
//  if (!arma::inv_sympd(Ei, E)) {
//   throw std::runtime_error("Initial Ei inversion failed.");
//  }
//
//  std::vector<double> cmodel(static_cast<std::size_t>(nmodels), 0.0);
//  std::vector<double> pis(static_cast<std::size_t>(nmodels), 0.0);
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//
//  arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
//
//  arma::mat cvbm(nt, nt, arma::fill::zeros);
//  arma::mat cvem(nt, nt, arma::fill::zeros);
//  arma::mat cvgm(nt, nt, arma::fill::zeros);
//
//  std::mt19937 gen(seed);
//
//  for (int it = 0; it < nit + nburn; ++it) {
//
//   std::fill(cmodel.begin(), cmodel.end(), 1.0);
//
//   arma::mat E_marker = E;
//
//   for (int t = 0; t < nt; ++t) {
//    E_marker(t, t) += adjE * std::max(G(t, t), 0.0);
//   }
//
//   arma::mat Ei_marker;
//   if (!arma::inv_sympd(Ei_marker, E_marker)) {
//    throw std::runtime_error("Adjusted marker E inversion failed.");
//   }
//
//   // ------------------------------------------------------------------------
//   // Marker updates. For correctness, these are serial because each marker
//   // update modifies the shared residual r_mat.
//   // ------------------------------------------------------------------------
//   if (method == 4) {
//    for (int isort = 0; isort < m; ++isort) {
//     const int i = order[static_cast<std::size_t>(isort)];
//
//     sampleBetaCPG_Mt_arma_ld(
//      i,
//      nt,
//      nmodels,
//      models,
//      cmodel,
//      pi,
//      Ei_marker,
//      Bi,
//      ww_mat,
//      r_mat,
//      beta_mat,
//      b_mat,
//      d_mat,
//      ld,
//      gen
//     );
//    }
//   }
//
//   if (updateB && method == 4) {
//    sampleB_cpg_arma_omp_csr(
//     nt,
//     static_cast<int>(nub),
//     B,
//     beta_mat,
//     d_mat,
//     ssb_prior_mat,
//     gen
//    );
//
//    if (!arma::inv_sympd(Bi, B)) {
//     throw std::runtime_error("Bi inversion failed.");
//    }
//   }
//
//   if (updateE) {
//    sampleE_cpg_arma_omp_csr(
//     nt,
//     static_cast<int>(nue),
//     E,
//     b_mat,
//     wy_mat,
//     r_mat,
//     sse_prior_mat,
//     YY_mat,
//     n,
//     gen,
//     true
//    );
//
//    if (!arma::inv_sympd(Ei, E)) {
//     throw std::runtime_error("Ei inversion failed.");
//    }
//   }
//
//   if (updatePi && method == 4) {
//    samplePi_cpg(cmodel, pi, gen);
//
//    if (it >= nburn) {
//     for (int k = 0; k < nmodels; ++k) {
//      pis[static_cast<std::size_t>(k)] += pi[static_cast<std::size_t>(k)];
//     }
//    }
//   }
//
//   computeG_cpg_arma_omp_csr(nt, b_mat, wy_mat, r_mat, n, G);
//
//   for (int t = 0; t < nt; ++t) {
//    vbs_mat(t, it) = B(t, t);
//    ves_mat(t, it) = E(t, t);
//    vgs_mat(t, it) = G(t, t);
//   }
//
//   if (it >= nburn) {
//    cvbm += B;
//    cvem += E;
//    cvgm += G;
//   }
//
//   if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//    nsamples += 1.0;
//
//    for (int t = 0; t < nt; ++t) {
//     for (int i = 0; i < m; ++i) {
//      if (d_mat(t, i) > 0) {
//       dm_mat(t, i) += 1.0;
//       bm_mat(t, i) += b_mat(t, i);
//      }
//     }
//    }
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Build result with same 20-slot structure
//  // --------------------------------------------------------------------------
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
//   result[16][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nmodels));
//   result[17][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nmodels));
//
//   result[18][static_cast<std::size_t>(t)].resize(4);
//   result[19][static_cast<std::size_t>(t)].resize(2);
//  }
//
//  const double denom = std::max(nsamples, 1.0);
//  const double post_iter_denom = static_cast<double>(std::max(nit, 1));
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int i = 0; i < m; ++i) {
//    const std::size_t is = static_cast<std::size_t>(i);
//
//    result[0][ts][is] = bm_mat(t, i) / denom;
//    result[1][ts][is] = dm_mat(t, i) / denom;
//    result[2][ts][is] = wy_mat(t, i);
//    result[3][ts][is] = r_mat(t, i);
//    result[4][ts][is] = b_mat(t, i);
//    result[5][ts][is] = static_cast<double>(d_mat(t, i));
//    result[6][ts][is] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int it = 0; it < nit + nburn; ++it) {
//    const std::size_t its = static_cast<std::size_t>(it);
//
//    result[7][ts][its] = vbs_mat(t, it);
//    result[8][ts][its] = vgs_mat(t, it);
//    result[9][ts][its] = ves_mat(t, it);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   const std::size_t t1s = static_cast<std::size_t>(t1);
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    const std::size_t t2s = static_cast<std::size_t>(t2);
//
//    result[10][t1s][t2s] = cvbm(t1, t2) / post_iter_denom;
//    result[11][t1s][t2s] = cvgm(t1, t2) / post_iter_denom;
//    result[12][t1s][t2s] = cvem(t1, t2) / post_iter_denom;
//    result[13][t1s][t2s] = B(t1, t2);
//    result[14][t1s][t2s] = G(t1, t2);
//    result[15][t1s][t2s] = E(t1, t2);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int k = 0; k < nmodels; ++k) {
//    const std::size_t ks = static_cast<std::size_t>(k);
//    result[16][ts][ks] = pi[ks];
//    result[17][ts][ks] = pis[ks] / post_iter_denom;
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
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



// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
//
// #include "cpg_samplers.h"
// #include "distributions.h"
//
// #include <algorithm>
// #include <cmath>
// #include <limits>
// #include <numeric>
// #include <random>
// #include <stdexcept>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // -----------------------------------------------------------------------------
// // Shared LD structure
// // -----------------------------------------------------------------------------
// struct LDNeighbor {
//  int j;
//  //double xij;  // X_i'X_j
//  float xij;
// };
//
//
// using LDNeighbors = std::vector<std::vector<LDNeighbor>>;
//
// inline LDNeighbors build_symmetric_ld_neighbors(
//   int m,
//   const std::vector<double>& row_ptr,
//   const std::vector<int>& col_idx,
//   const std::vector<double>& values,
//   bool col_idx_one_based,
//   const std::vector<double>& xx
// ) {
//  if ((int)row_ptr.size() != m + 1) {
//   throw std::runtime_error("build_symmetric_ld_neighbors: row_ptr must have length m + 1.");
//  }
//
//  if (col_idx.size() != values.size()) {
//   throw std::runtime_error("build_symmetric_ld_neighbors: col_idx and values must have same length.");
//  }
//
//  if ((int)xx.size() != m) {
//   throw std::runtime_error("build_symmetric_ld_neighbors: xx must have length m.");
//  }
//
//  LDNeighbors neigh(m);
//
//  for (int i = 0; i < m; ++i) {
//   const std::size_t start = static_cast<std::size_t>(row_ptr[i]);
//   const std::size_t end   = static_cast<std::size_t>(row_ptr[i + 1]);
//
//   if (end < start || end > values.size()) {
//    throw std::runtime_error("build_symmetric_ld_neighbors: invalid row_ptr values.");
//   }
//
//   if (!std::isfinite(xx[i]) || xx[i] < 0.0) {
//    throw std::runtime_error("build_symmetric_ld_neighbors: xx contains invalid values.");
//   }
//
//   for (std::size_t k = start; k < end; ++k) {
//    int j = col_idx[k];
//    if (col_idx_one_based) --j;
//
//    if (j < 0 || j >= m) {
//     throw std::runtime_error("build_symmetric_ld_neighbors: LD column index out of range.");
//    }
//
//    if (j == i) continue;
//
//    const double rij = values[k];
//
//    if (!std::isfinite(rij)) {
//     throw std::runtime_error("build_symmetric_ld_neighbors: LD value contains NaN/Inf.");
//    }
//
//    if (!std::isfinite(xx[j]) || xx[j] < 0.0) {
//     throw std::runtime_error("build_symmetric_ld_neighbors: xx contains invalid values.");
//    }
//
//    const double xij = rij * std::sqrt(xx[i] * xx[j]);
//
//    if (!std::isfinite(xij)) {
//     throw std::runtime_error("build_symmetric_ld_neighbors: computed X_i'X_j contains NaN/Inf.");
//    }
//
//    neigh[i].push_back({j, xij});
//    neigh[j].push_back({i, xij});
//   }
//  }
//
//  return neigh;
// }
//
// // -----------------------------------------------------------------------------
// // Variance component updates
// // -----------------------------------------------------------------------------
//
// inline void sampleB_cpg_arma_omp_csr(
//   int nt,
//   int nub,
//   arma::mat& B,
//   const arma::mat& beta,
//   const arma::Mat<int>& d,
//   const arma::mat& ssb_prior,
//   std::mt19937& gen
// ) {
//  if (nub <= nt - 1) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: nub must be > nt - 1.");
//  }
//
//  if (!beta.is_finite()) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: beta contains NaN/Inf.");
//  }
//
//  //arma::mat S_post = ssb_prior;
//  arma::mat S_post = static_cast<double>(nub) * ssb_prior;
//  int n_active = 0;
//
//  for (arma::uword i = 0; i < beta.n_cols; ++i) {
//   bool active = false;
//
//   for (int t = 0; t < nt; ++t) {
//    if (d(t, i) > 0) {
//     active = true;
//     break;
//    }
//   }
//
//   if (active) {
//    const arma::vec bi = beta.col(i);
//    S_post += bi * bi.t();
//    ++n_active;
//   }
//  }
//
//  S_post = 0.5 * (S_post + S_post.t());
//
//  if (!S_post.is_finite()) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: S_post contains NaN/Inf.");
//  }
//
//  arma::mat L;
//  bool ok = arma::chol(L, S_post, "lower");
//
//  double jitter = 1e-10;
//  int tries = 0;
//
//  while (!ok && tries < 8) {
//   S_post.diag() += jitter;
//   S_post = 0.5 * (S_post + S_post.t());
//   ok = arma::chol(L, S_post, "lower");
//   jitter *= 10.0;
//   ++tries;
//  }
//
//  if (!ok) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: S_post is not SPD after jitter.");
//  }
//
//  const unsigned int df_post =
//   static_cast<unsigned int>(nub + std::max(n_active, 1));
//
//  B = rinvwishart_cpp(df_post, S_post, gen);
//
//  B = 0.5 * (B + B.t());
//  B.diag() = arma::clamp(B.diag(), 1e-12, arma::datum::inf);
//
//  if (!B.is_finite()) {
//   throw std::runtime_error("sampleB_cpg_arma_omp_csr: sampled B contains NaN/Inf.");
//  }
// }
//
// inline void sampleB_active_shrunk_csr(
//   int nt,
//   int m,
//   int nub,
//   arma::mat& B,
//   const arma::Mat<int>& d,
//   const arma::mat& b,
//   const arma::mat& ssb_prior,
//   std::mt19937& gen,
//   double corr_df = 20.0,
//   double ridge = 1e-8
// ) {
//  arma::vec var_draw(nt, arma::fill::zeros);
//  arma::mat R(nt, nt, arma::fill::eye);
//
//  // ----------------------------
//  // 1. Trait-specific variances
//  // ----------------------------
//  for (int t = 0; t < nt; ++t) {
//   double ss = 0.0;
//   int n_active = 0;
//
//   for (int i = 0; i < m; ++i) {
//    if (d(t, i) > 0) {
//     ss += b(t, i) * b(t, i);
//     ++n_active;
//    }
//   }
//
//   const double scale = ss + nub * ssb_prior(t, t);
//   const double df = static_cast<double>(nub + n_active);
//
//   if (!std::isfinite(scale) || scale <= 0.0) {
//    throw std::runtime_error("sampleB_active_shrunk_csr: invalid variance scale.");
//   }
//
//   std::chi_squared_distribution<double> rchisq(std::max(df, 1.0));
//   const double x = std::max(rchisq(gen), 1e-300);
//
//   var_draw(t) = std::max(scale / x, ridge);
//  }
//
//  // ----------------------------
//  // 2. Pairwise correlations
//  // ----------------------------
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = t1 + 1; t2 < nt; ++t2) {
//    double s12 = 0.0;
//    int n_shared = 0;
//
//    for (int i = 0; i < m; ++i) {
//     if (d(t1, i) > 0 && d(t2, i) > 0) {
//      s12 += b(t1, i) * b(t2, i);
//      ++n_shared;
//     }
//    }
//
//    double r = 0.0;
//
//    if (n_shared > 0) {
//     const double cov_hat = s12 / static_cast<double>(n_shared);
//     const double denom = std::sqrt(var_draw(t1) * var_draw(t2));
//
//     if (denom > 0.0 && std::isfinite(denom)) {
//      r = cov_hat / denom;
//     }
//
//     const double corr_shrink =
//      static_cast<double>(n_shared) /
//       (static_cast<double>(n_shared) + corr_df);
//
//     r *= corr_shrink;
//    }
//
//    r = std::max(-0.95, std::min(0.95, r));
//
//    R(t1, t2) = r;
//    R(t2, t1) = r;
//   }
//  }
//
//  // ----------------------------
//  // 3. Reconstruct B = D R D
//  // ----------------------------
//  arma::vec sd = arma::sqrt(var_draw);
//  arma::mat D = arma::diagmat(sd);
//
//  B = D * R * D;
//  B.diag() += ridge;
//  B = 0.5 * (B + B.t());
//
//  // SPD safeguard
//  arma::vec eigval;
//  arma::mat eigvec;
//
//  if (!arma::eig_sym(eigval, eigvec, B)) {
//   B = arma::diagmat(arma::clamp(B.diag(), ridge, arma::datum::inf));
//   return;
//  }
//
//  eigval = arma::clamp(eigval, ridge, arma::datum::inf);
//  B = eigvec * arma::diagmat(eigval) * eigvec.t();
//  B = 0.5 * (B + B.t());
// }
//
//
// // inline void sampleE_cpg_arma_omp_csr(
// //   int nt,
// //   int nue,
// //   arma::mat& E,
// //   const arma::mat& b,
// //   const arma::mat& wy,
// //   const arma::mat& r,
// //   const arma::mat& sse_prior,
// //   const arma::vec& yy,
// //   const std::vector<int>& n,
// //   std::mt19937& gen
// // ) {
// //  E.zeros();
// //
// //  for (int t = 0; t < nt; ++t) {
// //   const double fitted_quad = arma::dot(b.row(t), r.row(t) + wy.row(t));
// //   const double sse = yy(t) - fitted_quad;
// //   const double prior_scale = nue * sse_prior(t, t);
// //   const double scale = sse + prior_scale;
// //
// //   if (!std::isfinite(scale) || scale <= 0.0) {
// //    int active = 0;
// //    for (arma::uword i = 0; i < b.n_cols; ++i) {
// //     if (b(t, i) != 0.0) ++active;
// //    }
// //
// //    Rcpp::Rcerr
// //    << "\nDIAGNOSTIC sampleE failure"
// //    << "\ntrait        = " << t
// //    << "\nactive       = " << active
// //    << "\nyy           = " << yy(t)
// //    << "\nfitted_quad  = " << fitted_quad
// //    << "\nsse          = " << sse
// //    << "\nprior_scale  = " << prior_scale
// //    << "\nscale        = " << scale
// //    << "\nE_old        = " << E(t, t)
// //    << "\n";
// //
// //    throw std::runtime_error(
// //      "DIAGNOSTIC sampleE_cpg_arma_omp: invalid residual scale. trait=" +
// //       std::to_string(t) +
// //       ", active=" + std::to_string(active) +
// //       ", yy=" + std::to_string(yy(t)) +
// //       ", fitted_quad=" + std::to_string(fitted_quad) +
// //       ", sse=" + std::to_string(sse) +
// //       ", prior_scale=" + std::to_string(prior_scale) +
// //       ", scale=" + std::to_string(scale)
// //    );
// //   }
// //   std::chi_squared_distribution<double> rchisq(n[t] + nue);
// //   const double chi2 = rchisq(gen);
// //
// //   E(t, t) = scale / chi2;
// //  }
// // }
//
// inline void make_spd_symmetric(
//   arma::mat& A,
//   double ridge = 1e-10
// ) {
//  A = 0.5 * (A + A.t());
//
//  arma::vec eigval;
//  arma::mat eigvec;
//
//  if (!arma::eig_sym(eigval, eigvec, A)) {
//   A = arma::diagmat(arma::clamp(A.diag(), ridge, arma::datum::inf));
//   return;
//  }
//
//  eigval = arma::clamp(eigval, ridge, arma::datum::inf);
//  A = eigvec * arma::diagmat(eigval) * eigvec.t();
//  A = 0.5 * (A + A.t());
// }
//
//
// inline void sampleE_cpg_arma_omp_csr(
//   int nt,
//   int nue,
//   arma::mat& E,
//   const arma::mat& b,
//   const arma::mat& wy,
//   const arma::mat& r,
//   const arma::mat& sse_prior,
//   const arma::mat& YY,
//   const std::vector<int>& n,
//   std::mt19937& gen,
//   bool update_full_E = true,
//   double ridge = 1e-10
// ) {
//  if ((int)n.size() != nt) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: n must have length nt.");
//  }
//
//  if ((int)YY.n_rows != nt || (int)YY.n_cols != nt) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: YY must be nt x nt.");
//  }
//
//  if ((int)b.n_rows != nt || (int)wy.n_rows != nt || (int)r.n_rows != nt) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: trait dimension mismatch.");
//  }
//
//  if (b.n_cols != wy.n_cols || b.n_cols != r.n_cols) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: marker dimension mismatch.");
//  }
//
//  // --------------------------------------------------------------------------
//  // Diagonal-only update.
//  // This preserves the old behavior, but uses YY(t,t).
//  // --------------------------------------------------------------------------
//  if (!update_full_E) {
//   E.zeros(nt, nt);
//
//   for (int t = 0; t < nt; ++t) {
//    const double fitted_quad = arma::dot(b.row(t), r.row(t) + wy.row(t));
//    const double sse = YY(t, t) - fitted_quad;
//    const double scale = sse + static_cast<double>(nue) * sse_prior(t, t);
//
//    if (!std::isfinite(scale) || scale <= 0.0) {
//     int active = 0;
//     for (arma::uword i = 0; i < b.n_cols; ++i) {
//      if (b(t, i) != 0.0) ++active;
//     }
//
//     Rcpp::Rcerr
//     << "\nDIAGNOSTIC sampleE diagonal failure"
//     << "\ntrait        = " << t
//     << "\nactive       = " << active
//     << "\nYY_diag      = " << YY(t, t)
//     << "\nfitted_quad  = " << fitted_quad
//     << "\nsse          = " << sse
//     << "\nprior_scale  = " << static_cast<double>(nue) * sse_prior(t, t)
//     << "\nscale        = " << scale
//     << "\n";
//
//     throw std::runtime_error("sampleE_cpg_arma_omp_csr: invalid diagonal residual scale.");
//    }
//
//    std::chi_squared_distribution<double> rchisq(n[t] + nue);
//    const double chi2 = std::max(rchisq(gen), 1e-300);
//
//    E(t, t) = std::max(scale / chi2, ridge);
//   }
//
//   return;
//  }
//
//  // --------------------------------------------------------------------------
//  // Full residual covariance update.
//  //
//  // Residual cross-product:
//  //
//  //   Se = Y'Y - B'X'Y - Y'X B + B'X'X B
//  //
//  // Here:
//  //   wy[t] = X'y_t
//  //   r[t]  = X'y_t - X'X b_t
//  //   therefore X'X b_t = wy[t] - r[t]
//  //
//  // For traits t1,t2:
//  //
//  //   Se(t1,t2)
//  //    = YY(t1,t2)
//  //      - b1' wy2
//  //      - b2' wy1
//  //      + b1' X'X b2
//  //
//  // with b1' X'X b2 = b1' (wy2 - r2).
//  // --------------------------------------------------------------------------
//
//  arma::mat Se(nt, nt, arma::fill::zeros);
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = t1; t2 < nt; ++t2) {
//
//    const double term12 =
//     YY(t1, t2)
//    - arma::dot(b.row(t1), wy.row(t2))
//    - arma::dot(b.row(t2), wy.row(t1))
//    + arma::dot(b.row(t1), wy.row(t2) - r.row(t2));
//
//    double se = term12;
//
//    if (t1 != t2) {
//     const double term21 =
//      YY(t2, t1)
//     - arma::dot(b.row(t2), wy.row(t1))
//     - arma::dot(b.row(t1), wy.row(t2))
//     + arma::dot(b.row(t2), wy.row(t1) - r.row(t1));
//
//     se = 0.5 * (term12 + term21);
//    }
//
//    Se(t1, t2) = se;
//    Se(t2, t1) = se;
//   }
//  }
//
//  Se = 0.5 * (Se + Se.t());
//
//  for (int t = 0; t < nt; ++t) {
//   if (!std::isfinite(Se(t, t)) || Se(t, t) <= 0.0) {
//    int active = 0;
//    for (arma::uword i = 0; i < b.n_cols; ++i) {
//     if (b(t, i) != 0.0) ++active;
//    }
//
//    Rcpp::Rcerr
//    << "\nDIAGNOSTIC full sampleE failure"
//    << "\ntrait        = " << t
//    << "\nactive       = " << active
//    << "\nYY_diag      = " << YY(t, t)
//    << "\nSe_diag      = " << Se(t, t)
//    << "\nprior_diag   = " << sse_prior(t, t)
//    << "\n";
//
//    throw std::runtime_error(
//      "sampleE_cpg_arma_omp_csr: full residual cross-product has invalid diagonal."
//    );
//   }
//  }
//
//  arma::mat S_post = Se + static_cast<double>(nue) * sse_prior;
//  S_post = 0.5 * (S_post + S_post.t());
//
//  if (!S_post.is_finite()) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: S_post contains NaN/Inf.");
//  }
//
//  arma::vec eigval;
//  arma::mat eigvec;
//
//  if (!arma::eig_sym(eigval, eigvec, S_post)) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: eig_sym failed for S_post.");
//  }
//
//  const double min_eig = arma::min(eigval);
//
//  if (!std::isfinite(min_eig) || min_eig <= -1e-6) {
//   Rcpp::Rcerr
//   << "\nDIAGNOSTIC full sampleE S_post not PSD"
//   << "\nmin_eigenvalue = " << min_eig
//   << "\nS_post diag    = " << S_post.diag().t()
//   << "\n";
//
//   throw std::runtime_error(
//     "sampleE_cpg_arma_omp_csr: full S_post is not positive semidefinite."
//   );
//  }
//
//  eigval = arma::clamp(eigval, ridge, arma::datum::inf);
//  S_post = eigvec * arma::diagmat(eigval) * eigvec.t();
//  S_post = 0.5 * (S_post + S_post.t());
//
//  const unsigned int df_post = static_cast<unsigned int>(n[0] + nue);
//
//  E = rinvwishart_cpp(df_post, S_post, gen);
//
//  E = 0.5 * (E + E.t());
//  E.diag() = arma::clamp(E.diag(), ridge, arma::datum::inf);
//
//  make_spd_symmetric(E, ridge);
//
//  if (!E.is_finite()) {
//   throw std::runtime_error("sampleE_cpg_arma_omp_csr: sampled E contains NaN/Inf.");
//  }
// }
//
// inline void computeG_cpg_arma_omp_csr(
//   int nt,
//   const arma::mat& b,
//   const arma::mat& wy,
//   const arma::mat& r,
//   const std::vector<int>& n,
//   arma::mat& G
// ) {
//  const arma::mat Q = wy - r;
//
//  G.set_size(nt, nt);
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = t1; t2 < nt; ++t2) {
//    const double denom = std::sqrt((double)n[t1] * (double)n[t2]);
//
//    double s12 = arma::dot(b.row(t1), Q.row(t2));
//
//    if (t1 != t2) {
//     const double s21 = arma::dot(b.row(t2), Q.row(t1));
//     s12 = 0.5 * (s12 + s21);
//    }
//
//    const double gij = s12 / denom;
//    G(t1, t2) = gij;
//    G(t2, t1) = gij;
//   }
//  }
// }
//
// // -----------------------------------------------------------------------------
// // Residual update using shared LD neighbors
// // -----------------------------------------------------------------------------
//
// inline void rebuild_residual_from_b_ld(
//   int nt,
//   int m,
//   const arma::mat& wy,
//   const arma::mat& ww,
//   const arma::mat& b,
//   arma::mat& r,
//   const LDNeighbors& ld_neighbors
// ) {
//  r = wy;
//
//  for (int i = 0; i < m; ++i) {
//   for (int t = 0; t < nt; ++t) {
//    const double bi = b(t, i);
//    if (bi == 0.0) continue;
//
//    // Diagonal contribution: X_i'X_i b_i
//    r(t, i) -= ww(t, i) * bi;
//
//    // Off-diagonal contribution: X_j'X_i b_i
//    for (const LDNeighbor& nb : ld_neighbors[i]) {
//     r(t, nb.j) -= nb.xij * bi;
//    }
//   }
//  }
// }
//
// // -----------------------------------------------------------------------------
// // Marker update using shared LD neighbors
// // -----------------------------------------------------------------------------
//
//
// inline void sampleBetaCPG_Mt_latent_arma_ld(
//   int i,
//   int nt,
//   int nmodels,
//   const std::vector<std::vector<int>>& models,
//   std::vector<double>& cmodel_local,
//   const std::vector<double>& pi,
//   const arma::mat& Ei,
//   const arma::mat& Bi,
//   const arma::mat& ww,
//   arma::mat& r,
//   arma::mat& beta,
//   arma::mat& b,
//   arma::Mat<int>& d,
//   const LDNeighbors& ld_neighbors,
//   std::mt19937& gen
// ) {
//  arma::vec rhs(nt), rhs_base(nt), z(nt), beta_new(nt), b_new(nt);
//  arma::mat C(nt, nt), L(nt, nt);
//
//  std::vector<double> loglik(nmodels), w(nmodels);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm(0.0, 1.0);
//
//  for (int t = 0; t < nt; ++t) {
//   rhs_base(t) = Ei(t, t) * (r(t, i) + ww(t, i) * b(t, i));
//  }
//
//  for (int k = 0; k < nmodels; ++k) {
//   if (pi[k] <= 0.0) {
//    loglik[k] = -INFINITY;
//    continue;
//   }
//
//   rhs.zeros();
//   C = Bi;
//
//   const auto& mk = models[k];
//
//   for (int t = 0; t < nt; ++t) {
//    if (mk[t]) {
//     rhs(t) = rhs_base(t);
//     C(t, t) += ww(t, i) * Ei(t, t);
//    }
//   }
//
//   if (!arma::chol(L, C, "lower")) {
//    loglik[k] = -INFINITY;
//    continue;
//   }
//
//   const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
//
//   arma::vec y    = arma::solve(arma::trimatl(L), rhs);
//   arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
//
//   const double quad = arma::dot(rhs, mean);
//
//   loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
//  }
//
//  const double max_log = *std::max_element(loglik.begin(), loglik.end());
//
//  double sumw = 0.0;
//  for (int k = 0; k < nmodels; ++k) {
//   w[k] = std::isfinite(loglik[k]) ? std::exp(loglik[k] - max_log) : 0.0;
//   sumw += w[k];
//  }
//
//  if (!std::isfinite(sumw) || sumw <= 0.0) {
//   throw std::runtime_error("sampleBetaCPG_Mt_latent_arma_ld: invalid model weights.");
//  }
//
//  const double u = runif(gen) * sumw;
//
//  int mselect = nmodels - 1;
//  double cum = 0.0;
//
//  for (int k = 0; k < nmodels; ++k) {
//   cum += w[k];
//   if (u <= cum) {
//    mselect = k;
//    break;
//   }
//  }
//
//  cmodel_local[mselect] += 1.0;
//
//  const auto& msel = models[mselect];
//
//  rhs.zeros();
//  C = Bi;
//
//  for (int t = 0; t < nt; ++t) {
//   d(t, i) = msel[t];
//
//   if (msel[t]) {
//    rhs(t) = rhs_base(t);
//    C(t, t) += ww(t, i) * Ei(t, t);
//   }
//  }
//
//  if (!arma::chol(L, C, "lower")) {
//   throw std::runtime_error("sampleBetaCPG_Mt_latent_arma_ld: selected model chol failed.");
//  }
//
//  arma::vec y    = arma::solve(arma::trimatl(L), rhs);
//  arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
//
//  for (int t = 0; t < nt; ++t) {
//   z(t) = norm(gen);
//  }
//
//  beta_new = mean + arma::solve(arma::trimatu(L.t()), z);
//
//  b_new.zeros();
//  for (int t = 0; t < nt; ++t) {
//   if (msel[t]) b_new(t) = beta_new(t);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const double diff = b_new(t) - b(t, i);
//
//   if (diff != 0.0) {
//    // Diagonal contribution
//    r(t, i) -= ww(t, i) * diff;
//
//    // Off-diagonal contribution
//    for (const LDNeighbor& nb : ld_neighbors[i]) {
//     r(t, nb.j) -= nb.xij * diff;
//    }
//   }
//  }
//
//  beta.col(i) = beta_new;
//  b.col(i)    = b_new;
// }
//
// // inline void sampleBetaCPG_Mt_arma_ld(
// //   int i,
// //   int nt,
// //   int nmodels,
// //   const std::vector<std::vector<int>>& models,
// //   std::vector<double>& cmodel_local,
// //   const std::vector<double>& pi,
// //   const arma::mat& Ei,
// //   const arma::mat& Bi,
// //   const arma::mat& ww,
// //   arma::mat& r,
// //   arma::mat& beta,
// //   arma::mat& b,
// //   arma::Mat<int>& d,
// //   const LDNeighbors& ld_neighbors,
// //   std::mt19937& gen
// // ) {
// //  arma::vec rhs_base(nt, arma::fill::zeros);
// //  arma::vec beta_new(nt, arma::fill::zeros);
// //  arma::vec b_new(nt, arma::fill::zeros);
// //
// //  std::vector<double> loglik(nmodels, -INFINITY);
// //  std::vector<double> w(nmodels, 0.0);
// //
// //  std::uniform_real_distribution<double> runif(0.0, 1.0);
// //  std::normal_distribution<double> norm(0.0, 1.0);
// //
// //  // ------------------------------------------------------------
// //  // Base RHS using current EFFECTIVE effect
// //  //
// //  // r(t,i) is current residual score:
// //  //   r = X'y - X'X b
// //  //
// //  // r(t,i) + ww(t,i) * b(t,i) is the partial residual
// //  // with marker i added back.
// //  // ------------------------------------------------------------
// //  for (int t = 0; t < nt; ++t) {
// //   rhs_base(t) = Ei(t, t) * (r(t, i) + ww(t, i) * b(t, i));
// //  }
// //
// //  // ------------------------------------------------------------
// //  // 1. Compute model log posterior weights in active subspace
// //  // ------------------------------------------------------------
// //  for (int k = 0; k < nmodels; ++k) {
// //   if (pi[k] <= 0.0 || !std::isfinite(pi[k])) {
// //    loglik[k] = -INFINITY;
// //    continue;
// //   }
// //
// //   const auto& mk = models[k];
// //
// //   std::vector<int> active;
// //   active.reserve(nt);
// //
// //   for (int t = 0; t < nt; ++t) {
// //    if (mk[t]) active.push_back(t);
// //   }
// //
// //   const int q = static_cast<int>(active.size());
// //
// //   // Null model
// //   if (q == 0) {
// //    loglik[k] = std::log(pi[k]);
// //    continue;
// //   }
// //
// //   arma::uvec idx(q);
// //   for (int a = 0; a < q; ++a) {
// //    idx(a) = static_cast<arma::uword>(active[a]);
// //   }
// //
// //   arma::vec rhs = rhs_base.elem(idx);
// //   arma::mat C = Bi.submat(idx, idx);
// //
// //   for (int a = 0; a < q; ++a) {
// //    const int t = active[a];
// //    C(a, a) += ww(t, i) * Ei(t, t);
// //   }
// //
// //   arma::mat L;
// //   bool ok = arma::chol(L, C, "lower");
// //
// //   if (!ok) {
// //    double jitter = 1e-10;
// //    for (int tries = 0; tries < 8 && !ok; ++tries) {
// //     C.diag() += jitter;
// //     ok = arma::chol(L, C, "lower");
// //     jitter *= 10.0;
// //    }
// //   }
// //
// //   if (!ok) {
// //    loglik[k] = -INFINITY;
// //    continue;
// //   }
// //
// //   const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
// //
// //   arma::vec y = arma::solve(arma::trimatl(L), rhs);
// //   arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
// //
// //   const double quad = arma::dot(rhs, mean);
// //
// //   loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
// //  }
// //
// //  // ------------------------------------------------------------
// //  // 2. Stabilized sampling of model indicator
// //  // ------------------------------------------------------------
// //  const double max_log = *std::max_element(loglik.begin(), loglik.end());
// //
// //  double sumw = 0.0;
// //  for (int k = 0; k < nmodels; ++k) {
// //   w[k] = std::isfinite(loglik[k]) ? std::exp(loglik[k] - max_log) : 0.0;
// //   sumw += w[k];
// //  }
// //
// //  if (!std::isfinite(sumw) || sumw <= 0.0) {
// //   throw std::runtime_error("sampleBetaCPG_Mt_latent_arma_ld: invalid model weights.");
// //  }
// //
// //  const double u = runif(gen) * sumw;
// //
// //  int mselect = nmodels - 1;
// //  double cum = 0.0;
// //
// //  for (int k = 0; k < nmodels; ++k) {
// //   cum += w[k];
// //   if (u <= cum) {
// //    mselect = k;
// //    break;
// //   }
// //  }
// //
// //  cmodel_local[mselect] += 1.0;
// //
// //  const auto& msel = models[mselect];
// //
// //  for (int t = 0; t < nt; ++t) {
// //   d(t, i) = msel[t];
// //  }
// //
// //  // ------------------------------------------------------------
// //  // 3. Sample selected model in active subspace
// //  // ------------------------------------------------------------
// //  std::vector<int> active;
// //  active.reserve(nt);
// //
// //  for (int t = 0; t < nt; ++t) {
// //   if (msel[t]) active.push_back(t);
// //  }
// //
// //  const int q = static_cast<int>(active.size());
// //
// //  beta_new.zeros();
// //  b_new.zeros();
// //
// //  if (q > 0) {
// //   arma::uvec idx(q);
// //   for (int a = 0; a < q; ++a) {
// //    idx(a) = static_cast<arma::uword>(active[a]);
// //   }
// //
// //   arma::vec rhs = rhs_base.elem(idx);
// //   arma::mat C = Bi.submat(idx, idx);
// //
// //   for (int a = 0; a < q; ++a) {
// //    const int t = active[a];
// //    C(a, a) += ww(t, i) * Ei(t, t);
// //   }
// //
// //   arma::mat L;
// //   bool ok = arma::chol(L, C, "lower");
// //
// //   if (!ok) {
// //    double jitter = 1e-10;
// //    for (int tries = 0; tries < 8 && !ok; ++tries) {
// //     C.diag() += jitter;
// //     ok = arma::chol(L, C, "lower");
// //     jitter *= 10.0;
// //    }
// //   }
// //
// //   if (!ok) {
// //    throw std::runtime_error("sampleBetaCPG_Mt_latent_arma_ld: selected model chol failed.");
// //   }
// //
// //   arma::vec y = arma::solve(arma::trimatl(L), rhs);
// //   arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
// //
// //   arma::vec z(q);
// //   for (int a = 0; a < q; ++a) {
// //    z(a) = norm(gen);
// //   }
// //
// //   arma::vec draw = mean + arma::solve(arma::trimatu(L.t()), z);
// //
// //   for (int a = 0; a < q; ++a) {
// //    const int t = active[a];
// //    beta_new(t) = draw(a);
// //    b_new(t) = draw(a);
// //   }
// //  }
// //
// //  // ------------------------------------------------------------
// //  // 4. Residual update using EFFECTIVE effects
// //  // ------------------------------------------------------------
// //  for (int t = 0; t < nt; ++t) {
// //   const double diff = b_new(t) - b(t, i);
// //
// //   if (diff != 0.0) {
// //    // Diagonal contribution
// //    r(t, i) -= ww(t, i) * diff;
// //
// //    // Off-diagonal contribution
// //    for (const LDNeighbor& nb : ld_neighbors[i]) {
// //     r(t, nb.j) -= nb.xij * diff;
// //    }
// //   }
// //  }
// //
// //  beta.col(i) = beta_new;
// //  b.col(i) = b_new;
// // }
//
// inline void sampleBetaCPG_Mt_arma_ld(
//   int i,
//   int nt,
//   int nmodels,
//   const std::vector<std::vector<int>>& models,
//   std::vector<double>& cmodel_local,
//   const std::vector<double>& pi,
//   const arma::mat& Ei,
//   const arma::mat& Bi,
//   const arma::mat& ww,
//   arma::mat& r,
//   arma::mat& beta,
//   arma::mat& b,
//   arma::Mat<int>& d,
//   const LDNeighbors& ld_neighbors,
//   std::mt19937& gen
// ) {
//  arma::vec score(nt, arma::fill::zeros);
//  arma::vec z_full(nt, arma::fill::zeros);
//  arma::vec beta_new(nt, arma::fill::zeros);
//  arma::vec b_new(nt, arma::fill::zeros);
//
//  std::vector<double> loglik(nmodels, -INFINITY);
//  std::vector<double> w(nmodels, 0.0);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm(0.0, 1.0);
//
//  if ((int)Ei.n_rows != nt || (int)Ei.n_cols != nt) {
//   throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: Ei must be nt x nt.");
//  }
//
//  if ((int)Bi.n_rows != nt || (int)Bi.n_cols != nt) {
//   throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: Bi must be nt x nt.");
//  }
//
//  // --------------------------------------------------------------------------
//  // score(t) is the marker-specific partial score:
//  //
//  //   score_t = x_i' residual_without_i,t
//  //           = r(t,i) + x_i'x_i b(t,i)
//  //
//  // z_full is the approximate univariate partial regression estimate:
//  //
//  //   z_t = score_t / ww_t
//  //
//  // For correlated residuals, the likelihood contribution is approximated as:
//  //
//  //   z_i | beta_i, E ~ N(beta_i, D^{-1} E D^{-1})
//  //
//  // where D = diag(sqrt(ww_i)).
//  //
//  // Therefore:
//  //
//  //   precision = Bi + D Ei D
//  //   rhs       = D Ei D z
//  //
//  // If all ww are equal, this reduces to:
//  //
//  //   precision = Bi + ww Ei
//  //   rhs       = Ei score
//  // --------------------------------------------------------------------------
//
//  for (int t = 0; t < nt; ++t) {
//   if (!std::isfinite(ww(t, i)) || ww(t, i) <= 0.0) {
//    throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: ww contains invalid value.");
//   }
//
//   score(t) = r(t, i) + ww(t, i) * b(t, i);
//   z_full(t) = score(t) / ww(t, i);
//  }
//
//  // --------------------------------------------------------------------------
//  // 1. Compute model log posterior weights in active subspace
//  // --------------------------------------------------------------------------
//
//  for (int k = 0; k < nmodels; ++k) {
//   if (pi[k] <= 0.0 || !std::isfinite(pi[k])) {
//    loglik[k] = -INFINITY;
//    continue;
//   }
//
//   const auto& mk = models[k];
//
//   std::vector<int> active;
//   active.reserve(nt);
//
//   for (int t = 0; t < nt; ++t) {
//    if (mk[t]) active.push_back(t);
//   }
//
//   const int q = static_cast<int>(active.size());
//
//   if (q == 0) {
//    loglik[k] = std::log(pi[k]);
//    continue;
//   }
//
//   arma::uvec idx(q);
//   arma::vec sqrtw(q, arma::fill::zeros);
//
//   for (int a = 0; a < q; ++a) {
//    const int t = active[a];
//    idx(a) = static_cast<arma::uword>(t);
//    sqrtw(a) = std::sqrt(ww(t, i));
//   }
//
//   arma::mat Ei_sub = Ei.submat(idx, idx);
//   arma::mat C = Bi.submat(idx, idx);
//   arma::vec z = z_full.elem(idx);
//
//   arma::mat D = arma::diagmat(sqrtw);
//   arma::mat P = D * Ei_sub * D;
//
//   C += P;
//   arma::vec rhs = P * z;
//
//   C = 0.5 * (C + C.t());
//
//   arma::mat L;
//   bool ok = arma::chol(L, C, "lower");
//
//   if (!ok) {
//    double jitter = 1e-10;
//    for (int tries = 0; tries < 8 && !ok; ++tries) {
//     C.diag() += jitter;
//     C = 0.5 * (C + C.t());
//     ok = arma::chol(L, C, "lower");
//     jitter *= 10.0;
//    }
//   }
//
//   if (!ok) {
//    loglik[k] = -INFINITY;
//    continue;
//   }
//
//   const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
//
//   arma::vec y = arma::solve(arma::trimatl(L), rhs);
//   arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
//
//   const double quad = arma::dot(rhs, mean);
//
//   loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
//  }
//
//  // --------------------------------------------------------------------------
//  // 2. Sample model indicator
//  // --------------------------------------------------------------------------
//
//  const double max_log = *std::max_element(loglik.begin(), loglik.end());
//
//  double sumw = 0.0;
//  for (int k = 0; k < nmodels; ++k) {
//   w[k] = std::isfinite(loglik[k]) ? std::exp(loglik[k] - max_log) : 0.0;
//   sumw += w[k];
//  }
//
//  if (!std::isfinite(sumw) || sumw <= 0.0) {
//   throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: invalid model weights.");
//  }
//
//  const double u = runif(gen) * sumw;
//
//  int mselect = nmodels - 1;
//  double cum = 0.0;
//
//  for (int k = 0; k < nmodels; ++k) {
//   cum += w[k];
//   if (u <= cum) {
//    mselect = k;
//    break;
//   }
//  }
//
//  cmodel_local[mselect] += 1.0;
//
//  const auto& msel = models[mselect];
//
//  for (int t = 0; t < nt; ++t) {
//   d(t, i) = msel[t];
//  }
//
//  // --------------------------------------------------------------------------
//  // 3. Sample selected model in active subspace
//  // --------------------------------------------------------------------------
//
//  std::vector<int> active;
//  active.reserve(nt);
//
//  for (int t = 0; t < nt; ++t) {
//   if (msel[t]) active.push_back(t);
//  }
//
//  const int q = static_cast<int>(active.size());
//
//  beta_new.zeros();
//  b_new.zeros();
//
//  if (q > 0) {
//   arma::uvec idx(q);
//   arma::vec sqrtw(q, arma::fill::zeros);
//
//   for (int a = 0; a < q; ++a) {
//    const int t = active[a];
//    idx(a) = static_cast<arma::uword>(t);
//    sqrtw(a) = std::sqrt(ww(t, i));
//   }
//
//   arma::mat Ei_sub = Ei.submat(idx, idx);
//   arma::mat C = Bi.submat(idx, idx);
//   arma::vec z = z_full.elem(idx);
//
//   arma::mat D = arma::diagmat(sqrtw);
//   arma::mat P = D * Ei_sub * D;
//
//   C += P;
//   arma::vec rhs = P * z;
//
//   C = 0.5 * (C + C.t());
//
//   arma::mat L;
//   bool ok = arma::chol(L, C, "lower");
//
//   if (!ok) {
//    double jitter = 1e-10;
//    for (int tries = 0; tries < 8 && !ok; ++tries) {
//     C.diag() += jitter;
//     C = 0.5 * (C + C.t());
//     ok = arma::chol(L, C, "lower");
//     jitter *= 10.0;
//    }
//   }
//
//   if (!ok) {
//    throw std::runtime_error("sampleBetaCPG_Mt_arma_ld: selected model chol failed.");
//   }
//
//   arma::vec y = arma::solve(arma::trimatl(L), rhs);
//   arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
//
//   arma::vec eps(q);
//   for (int a = 0; a < q; ++a) {
//    eps(a) = norm(gen);
//   }
//
//   arma::vec draw = mean + arma::solve(arma::trimatu(L.t()), eps);
//
//   for (int a = 0; a < q; ++a) {
//    const int t = active[a];
//    beta_new(t) = draw(a);
//    b_new(t) = draw(a);
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // 4. Residual update using effective effects
//  // --------------------------------------------------------------------------
//
//  for (int t = 0; t < nt; ++t) {
//   const double diff = b_new(t) - b(t, i);
//
//   if (diff != 0.0) {
//    r(t, i) -= ww(t, i) * diff;
//
//    for (const LDNeighbor& nb : ld_neighbors[i]) {
//     r(t, nb.j) -= nb.xij * diff;
//    }
//   }
//  }
//
//  beta.col(i) = beta_new;
//  b.col(i) = b_new;
// }
//
// // -----------------------------------------------------------------------------
// // Main exported function
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> mtblr_cpg_omp_csr(
//   std::vector<std::vector<double>> wy,
//   std::vector<std::vector<double>> ww,
//   std::vector<std::vector<double>> yy,
//   std::vector<std::vector<double>> b_init,
//   std::vector<double> ld_row_ptr,
//   std::vector<int> ld_col_idx,
//   std::vector<double> ld_values,
//   bool ld_col_idx_one_based,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<std::vector<int>> models,
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
//   int seed,
//   int method
// ) {
//  const int nt = wy.size();
//  const int m  = wy[0].size();
//  const int nmodels = models.size();
//
//  if (nthin <= 0) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: nthin must be > 0.");
//  }
//
//  if (nt <= 0 || m <= 0 || nmodels <= 0) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: invalid dimensions.");
//  }
//
//  double nsamples = 0.0;
//
//  arma::mat wy_mat(nt, m, arma::fill::zeros);
//  arma::mat ww_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::mat beta_mat(nt, m, arma::fill::zeros);
//  arma::mat r_mat(nt, m, arma::fill::zeros);
//
//  //arma::vec yy_vec(nt, arma::fill::zeros);
//  arma::mat YY_mat(nt, nt, arma::fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  if ((int)yy.size() != nt) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: yy must have nt rows.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)yy[t].size() != nt) {
//    throw std::runtime_error("mtblr_cpg_omp_csr: yy must be nt x nt.");
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   //yy_vec(t) = yy[t];
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    YY_mat(t, t2) = yy[t][t2];
//   }
//
//   for (int i = 0; i < m; ++i) {
//    wy_mat(t, i) = wy[t][i];
//    ww_mat(t, i) = ww[t][i];
//    b_mat(t, i)  = b_init[t][i];
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//  }
//
//  if ((int)n.size() != nt) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: n must have length nt.");
//  }
//
//  if (n[0] <= 1) {
//   throw std::runtime_error("mtblr_cpg_omp_csr: n[0] must be > 1.");
//  }
//
//  for (int t = 1; t < nt; ++t) {
//   if (n[t] != n[0]) {
//    throw std::runtime_error(
//      "mtblr_cpg_omp_csr: current shared-LD scaling assumes equal n across traits."
//    );
//   }
//  }
//
//  std::vector<double> x2(m, 0.0);
//  std::vector<int> order(m);
//
//  for (int i = 0; i < m; ++i) {
//   double best = 0.0;
//
//   for (int t = 0; t < nt; ++t) {
//    if (ww_mat(t, i) > 0.0) {
//     const double bhat = wy_mat(t, i) / ww_mat(t, i);
//     best = std::max(best, bhat * bhat);
//    }
//   }
//
//   x2[i] = best;
//   order[i] = i;
//  }
//
//  std::sort(order.begin(), order.end(),
//            [&](int a, int b) { return x2[a] > x2[b]; });
//
//  std::vector<double> xx(m, 0.0);
//
//  for (int i = 0; i < m; ++i) {
//   xx[i] = ww_mat(0, i);
//  }
//
//  LDNeighbors ld_neighbors = build_symmetric_ld_neighbors(
//   m,
//   ld_row_ptr,
//   ld_col_idx,
//   ld_values,
//   ld_col_idx_one_based,
//   xx
//  );
//
//  // Initial residual: r = wy - LD_offdiag * b
//  rebuild_residual_from_b_ld(
//   nt,
//   m,
//   wy_mat,
//   ww_mat,
//   b_mat,
//   r_mat,
//   ld_neighbors
//  );
//
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
//  arma::mat Bi, Ei, G(nt, nt, arma::fill::zeros);
//
//  if (!arma::inv_sympd(Bi, B)) {
//   throw std::runtime_error("Initial Bi inversion failed.");
//  }
//
//  if (!arma::inv_sympd(Ei, E)) {
//   throw std::runtime_error("Initial Ei inversion failed.");
//  }
//
//  std::vector<double> cmodel(nmodels);
//  std::vector<double> pis(nmodels, 0.0);
//
//  std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
//  std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
//
//  std::vector<std::vector<double>> ves(nt, std::vector<double>(nit + nburn, 0.0));
//  std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit + nburn, 0.0));
//  std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit + nburn, 0.0));
//
//  std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
//  std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
//  std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));
//
//  std::mt19937 gen(seed);
//
//  for (int it = 0; it < nit + nburn; ++it) {
//
//   std::fill(cmodel.begin(), cmodel.end(), 1.0);
//
//   arma::mat E_marker = E;
//
//   for (int t = 0; t < nt; ++t) {
//    E_marker(t, t) += adjE * std::max(G(t, t), 0.0);
//   }
//
//   arma::mat Ei_marker;
//   if (!arma::inv_sympd(Ei_marker, E_marker)) {
//    throw std::runtime_error("Adjusted marker E inversion failed.");
//   }
//
//   // ----------------------------
//   // Marker updates
//   // ----------------------------
//   // The new LD object is symmetric and shared across traits.
//   // For correctness, marker updates are currently serial.
//   if (method == 4) {
//    for (int isort = 0; isort < m; ++isort) {
//     int i = order[isort];
//     sampleBetaCPG_Mt_arma_ld(
//      i,
//      nt,
//      nmodels,
//      models,
//      cmodel,
//      pi,
//      Ei_marker,
//      //Ei,
//      Bi,
//      ww_mat,
//      r_mat,
//      beta_mat,
//      b_mat,
//      d_mat,
//      ld_neighbors,
//      gen
//     );
//    }
//   }
//
//   // ----------------------------
//   // Sample B / E
//   // ----------------------------
//   if (updateB && method == 4) {
//    sampleB_cpg_arma_omp_csr(
//     nt,
//     nub,
//     B,
//     beta_mat,
//     d_mat,
//     ssb_prior_mat,
//     gen
//    );
//
//    if (!arma::inv_sympd(Bi, B)) {
//     throw std::runtime_error("Bi inversion failed.");
//    }
//   }
//   if (updateE) {
//    // sampleE_cpg_arma_omp_csr(
//    //  nt,
//    //  nue,
//    //  E,
//    //  b_mat,
//    //  wy_mat,
//    //  r_mat,
//    //  sse_prior_mat,
//    //  yy_vec,
//    //  n,
//    //  gen
//    // );
//    sampleE_cpg_arma_omp_csr(
//     nt,
//     nue,
//     E,
//     b_mat,
//     wy_mat,
//     r_mat,
//     sse_prior_mat,
//     YY_mat,
//     n,
//     gen,
//     true
//    );
//    if (!arma::inv_sympd(Ei, E)) {
//     throw std::runtime_error("Ei inversion failed.");
//    }
//   }
//
//   // ----------------------------
//   // Sample pi
//   // ----------------------------
//   if (updatePi && method == 4) {
//    samplePi_cpg(cmodel, pi, gen);
//
//    if (it >= nburn) {
//     for (int k = 0; k < nmodels; ++k) {
//      pis[k] += pi[k];
//     }
//    }
//   }
//
//   // ----------------------------
//   // Update G
//   // ----------------------------
//   computeG_cpg_arma_omp_csr(nt, b_mat, wy_mat, r_mat, n, G);
//
//   for (int t = 0; t < nt; ++t) {
//    vbs[t][it] = B(t, t);
//    ves[t][it] = E(t, t);
//    vgs[t][it] = G(t, t);
//   }
//
//   if (it >= nburn) {
//    for (int t1 = 0; t1 < nt; ++t1) {
//     for (int t2 = 0; t2 < nt; ++t2) {
//      cvbm[t1][t2] += B(t1, t2);
//      cvem[t1][t2] += E(t1, t2);
//      cvgm[t1][t2] += G(t1, t2);
//     }
//    }
//   }
//
//   if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//    nsamples += 1.0;
//
//    for (int t = 0; t < nt; ++t) {
//     for (int i = 0; i < m; ++i) {
//      if (d_mat(t, i) > 0) {
//       dm[t][i] += 1.0;
//       bm[t][i] += b_mat(t, i);
//      }
//     }
//    }
//   }
//  }
//
//  // ----------------------------
//  // Build result
//  // ----------------------------
//  std::vector<std::vector<std::vector<double>>> result(20);
//
//  for (int k = 0; k < 20; ++k) {
//   result[k].resize(nt);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][t].resize(m);
//   result[1][t].resize(m);
//   result[2][t].resize(m);
//   result[3][t].resize(m);
//   result[4][t].resize(m);
//   result[5][t].resize(m);
//   result[6][t].resize(m);
//
//   result[7][t].resize(nit + nburn);
//   result[8][t].resize(nit + nburn);
//   result[9][t].resize(nit + nburn);
//
//   result[10][t].resize(nt);
//   result[11][t].resize(nt);
//   result[12][t].resize(nt);
//   result[13][t].resize(nt);
//   result[14][t].resize(nt);
//   result[15][t].resize(nt);
//
//   result[16][t].resize(nmodels);
//   result[17][t].resize(nmodels);
//
//   result[18][t].resize(4);
//   result[19][t].resize(2);
//  }
//
//  const double denom = std::max(nsamples, 1.0);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    result[0][t][i] = bm[t][i] / denom;
//    result[1][t][i] = dm[t][i] / denom;
//    result[2][t][i] = wy_mat(t, i);
//    result[3][t][i] = r_mat(t, i);
//    result[4][t][i] = b_mat(t, i);
//    result[5][t][i] = d_mat(t, i);
//    result[6][t][i] = i;
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < nit + nburn; ++i) {
//    result[7][t][i] = vbs[t][i];
//    result[8][t][i] = vgs[t][i];
//    result[9][t][i] = ves[t][i];
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    result[10][t1][t2] = cvbm[t1][t2] / denom;
//    result[11][t1][t2] = cvgm[t1][t2] / denom;
//    result[12][t1][t2] = cvem[t1][t2] / denom;
//    result[13][t1][t2] = B(t1, t2);
//    result[14][t1][t2] = G(t1, t2);
//    result[15][t1][t2] = E(t1, t2);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < nmodels; ++i) {
//    result[16][t][i] = pi[i];
//    result[17][t][i] = pis[i] / denom;
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < 4; ++i) {
//    result[18][t][i] = 0.0;
//   }
//   for (int i = 0; i < 2; ++i) {
//    result[19][t][i] = 0.0;
//   }
//  }
//
//  return result;
// }
