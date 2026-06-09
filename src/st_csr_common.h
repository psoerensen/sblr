#pragma once

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

struct STLDCSR {
 std::vector<uint64_t> ptr;  // length m + 1
 std::vector<int> idx;       // neighbor marker index
 std::vector<float> xij;     // pre-scaled X_i'X_j

 // Diagnostics stored as data rather than printed.
 // This keeps the core C++ wrapper-neutral for R, Python, or command-line use.
 uint64_t input_nnz = 0;
 uint64_t symmetric_nnz = 0;
 double max_abs_r = 0.0;
 double max_abs_xij = 0.0;
};

inline void read_exact_file(
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

inline uint64_t parse_uint64_from_meta(
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

inline STLDCSR read_and_build_st_ld_csr(
  const std::string& prefix,
  int m,
  const std::vector<double>& xx
) {
 const std::string row_file  = prefix + ".row_ptr.u64.bin";
 const std::string col_file  = prefix + ".col_idx.u32.0based.bin";
 const std::string val_file  = prefix + ".values.f32.bin";
 const std::string meta_file = prefix + ".meta.txt";

 if (m <= 0) {
  throw std::runtime_error("read_and_build_st_ld_csr: m must be positive.");
 }

 if (static_cast<int>(xx.size()) != m) {
  throw std::runtime_error("read_and_build_st_ld_csr: xx must have length m.");
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
   nnz_u64 = parse_uint64_from_meta(line.substr(key_nnz.size()), "nnz");
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

 read_exact_file(row_file, row_ptr.data(), row_ptr.size() * sizeof(uint64_t));
 read_exact_file(col_file, col_idx_u32.data(), col_idx_u32.size() * sizeof(uint32_t));
 read_exact_file(val_file, values_r.data(), values_r.size() * sizeof(float));

 if (row_ptr[0] != 0 || row_ptr[static_cast<std::size_t>(m)] != nnz_u64) {
  throw std::runtime_error("Invalid LD row_ptr: expected 0-based row_ptr ending at nnz.");
 }

 for (int i = 0; i < m; ++i) {
  if (row_ptr[static_cast<std::size_t>(i + 1)] < row_ptr[static_cast<std::size_t>(i)]) {
   throw std::runtime_error("Invalid LD row_ptr: row pointers are not nondecreasing.");
  }

  if (!std::isfinite(xx[static_cast<std::size_t>(i)]) ||
      xx[static_cast<std::size_t>(i)] <= 0.0) {
   throw std::runtime_error(
     "read_and_build_st_ld_csr: xx contains invalid value at marker " +
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

 STLDCSR ld;
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

 double max_abs_r = 0.0;
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

   max_abs_r = std::max(max_abs_r, std::abs(rij));

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

 // Validate fill counts.
 for (int i = 0; i < m; ++i) {
  if (offset[static_cast<std::size_t>(i)] != ld.ptr[static_cast<std::size_t>(i + 1)]) {
   throw std::runtime_error("Internal LD CSR fill-count mismatch.");
  }
 }

 ld.input_nnz = nnz_u64;
 ld.symmetric_nnz = nnz_sym;
 ld.max_abs_r = max_abs_r;
 ld.max_abs_xij = max_abs_xij;

 return ld;
}

inline void rebuild_residual_st_csr(
  int m,
  const arma::rowvec& wy,
  const arma::rowvec& ww,
  const arma::rowvec& b,
  arma::rowvec& r,
  const STLDCSR& ld
) {
 r = wy;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double bi = b(iu);

  if (bi == 0.0) continue;

  r(iu) -= ww(iu) * bi;

  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   r(static_cast<arma::uword>(j)) -=
    static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * bi;
  }
 }
}

inline double clamp_prob(double x) {
 if (!std::isfinite(x)) {
  throw std::runtime_error("clamp_prob: probability is NaN/Inf.");
 }
 return std::min(std::max(x, 1e-300), 1.0 - 1e-12);
}

inline double inv_logit_stable(double x) {
 if (x >= 0.0) {
  const double z = std::exp(-x);
  return 1.0 / (1.0 + z);
 }
 const double z = std::exp(x);
 return z / (1.0 + z);
}

inline double log1pexp_stable(double x) {
 if (x > 40.0) return x;
 if (x < -40.0) return std::exp(x);
 return std::log1p(std::exp(x));
}

inline double logit_prob(double p) {
 p = clamp_prob(p);
 return std::log(p) - std::log1p(-p);
}

inline double clamp_double(double x, double lo, double hi) {
 return std::min(std::max(x, lo), hi);
}

inline void sampleE_ST_csr(
  int m,
  double nue,
  double& ve,
  const arma::rowvec& b,
  const arma::rowvec& wy,
  const arma::rowvec& r,
  double sse_prior,
  double yy,
  int n,
  std::mt19937& gen
) {
 double b_dot_r_plus_wy = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  b_dot_r_plus_wy += b(iu) * (r(iu) + wy(iu));
 }

 const double sse = yy - b_dot_r_plus_wy;
 const double scale = sse + nue * sse_prior;

 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleE_ST_csr: invalid residual scale.");
 }

 std::chi_squared_distribution<double> rchisq(n + nue);
 const double chi2 = std::max(rchisq(gen), 1e-300);

 const double ve_new = scale / chi2;

 if (!std::isfinite(ve_new) || ve_new <= 0.0) {
  throw std::runtime_error("sampleE_ST_csr: sampled ve is invalid.");
 }

 ve = std::max(ve_new, 1e-12);
}

inline double computeG_ST_csr(
  const arma::rowvec& b,
  const arma::rowvec& wy,
  const arma::rowvec& r,
  int n
) {
 double ssg = 0.0;
 const arma::uword m = b.n_elem;

 for (arma::uword i = 0; i < m; ++i) {
  ssg += b(i) * (wy(i) - r(i));
 }

 return ssg / static_cast<double>(n);
}
