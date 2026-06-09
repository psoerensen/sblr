// =============================================================================
// mt_sufficient_statistics.cpp
// =============================================================================

#include "mt_sufficient_statistics_core.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// =============================================================================
// BED helpers
// =============================================================================
//
// BED coding assumed:
//   00 -> dosage 2
//   01 -> missing
//   10 -> dosage 1
//   11 -> dosage 0
//
// =============================================================================

inline double decode_bed_dosage_from_marker(
  const unsigned char* buf,
  int bed_row
) {
 const unsigned char x = buf[bed_row >> 2];
 const unsigned int code = (x >> (2 * (bed_row & 3))) & 3u;

 if (code == 0u) return 2.0;
 if (code == 2u) return 1.0;
 if (code == 3u) return 0.0;
 return std::numeric_limits<double>::quiet_NaN(); // missing
}

inline double compute_af_from_packed_marker(
  const unsigned char* buf,
  const int* rows0,
  int n_used
) {
 double dosage_sum = 0.0;
 int n_obs = 0;

 if (!rows0) {
  const int nbytes = (n_used + 3) / 4;

  for (int kb = 0; kb < nbytes; ++kb) {
   unsigned char x = buf[kb];
   const int jbase = kb << 2;

   for (int pos = 0; pos < 4; ++pos) {
    const int j = jbase + pos;

    if (j >= n_used) {
     break;
    }

    const unsigned int code = x & 3u;
    x >>= 2;

    if (code == 1u) {
     continue; // missing
    }

    if (code == 0u) {
     dosage_sum += 2.0;
    } else if (code == 2u) {
     dosage_sum += 1.0;
    }

    ++n_obs;
   }
  }
 } else {
  for (int k = 0; k < n_used; ++k) {
   const int bed_row = rows0[k];
   const double dosage = decode_bed_dosage_from_marker(buf, bed_row);

   if (!std::isfinite(dosage)) {
    continue;
   }

   dosage_sum += dosage;
   ++n_obs;
  }
 }

 if (n_obs <= 0) {
  return 0.0;
 }

 return dosage_sum / (2.0 * static_cast<double>(n_obs));
}

// =============================================================================
// BED sufficient-statistics core
// =============================================================================
//
// If scale = true:
//   x = (dosage - 2p) / sqrt(2p(1-p))
//   missing -> 0 after standardization, i.e. mean-imputed standardized value.
//
// If scale = false:
//   x = raw dosage, missing -> 2p.
//   This is uncentered raw dosage coding.
//
// y_flat is n_used x nt, row-major: y_flat[k * nt + t], where k indexes
// selected individuals, not necessarily BED row number.
//
// cls is 1-based marker index into the BED file.
// rows0 is 0-based individual index into the BED file. If rows0 is null, all
// rows 0:(n_bed - 1) are used in order.
//
// If use_af = false, allele frequencies are computed over the selected rows
// and stored in af_used.
//
// =============================================================================

void bed_xtx_diag_xty_core(
  const char* file,
  int n_bed,
  const int* cls,
  const double* af,
  bool use_af,
  double* af_used,
  const int* rows0,
  int n_used,
  int m,
  int nt,
  bool scale,
  const double* y_flat,
  double* xty_flat,
  double* xx,
  int nthreads,
  int MG,
  int JB,
  int TB
) {
 if (!file) {
  throw std::runtime_error("file path pointer must not be null.");
 }

 FILE* file_stream = std::fopen(file, "rb");
 if (!file_stream) {
  throw std::runtime_error("Could not open BED file.");
 }

 if (n_bed <= 0) {
  std::fclose(file_stream);
  throw std::runtime_error("n_bed must be positive.");
 }

 if (n_used <= 0) {
  std::fclose(file_stream);
  throw std::runtime_error("n_used must be positive.");
 }

 if (m <= 0) {
  std::fclose(file_stream);
  throw std::runtime_error("m must be positive.");
 }

 if (nt <= 0) {
  std::fclose(file_stream);
  throw std::runtime_error("nt must be positive.");
 }

 if (!af_used) {
  std::fclose(file_stream);
  throw std::runtime_error("af_used must not be null.");
 }

 if (use_af && !af) {
  std::fclose(file_stream);
  throw std::runtime_error("af must not be null when use_af is true.");
 }

 if (MG <= 0) MG = 64;
 if (JB <= 0) JB = 1024;
 if (TB <= 0) TB = 32;

 // Validate PLINK BED header: magic bytes + SNP-major mode.
 unsigned char magic[3];
 if (std::fread(magic, sizeof(unsigned char), 3, file_stream) != 3) {
  std::fclose(file_stream);
  throw std::runtime_error("Could not read BED header.");
 }

 if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
  std::fclose(file_stream);
  throw std::runtime_error(
    "BED file must be PLINK SNP-major format: expected magic bytes 0x6c 0x1b 0x01."
  );
 }

 const std::size_t nbytes = static_cast<std::size_t>((n_bed + 3) / 4);

 std::memset(xty_flat, 0, static_cast<std::size_t>(m) * nt * sizeof(double));
 std::memset(xx, 0, static_cast<std::size_t>(m) * sizeof(double));
 std::memset(af_used, 0, static_cast<std::size_t>(m) * sizeof(double));

 unsigned char* block_buffer =
  static_cast<unsigned char*>(std::malloc(static_cast<std::size_t>(MG) * nbytes));

 if (!block_buffer) {
  std::fclose(file_stream);
  throw std::runtime_error("Failed to allocate BED block buffer.");
 }

 std::vector<double> map0(static_cast<std::size_t>(MG));
 std::vector<double> map1(static_cast<std::size_t>(MG));
 std::vector<double> map2(static_cast<std::size_t>(MG));
 std::vector<double> map3(static_cast<std::size_t>(MG));

#ifdef _OPENMP
 if (nthreads > 0) {
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
 }
#endif

 try {
  for (int i0 = 0; i0 < m; i0 += MG) {
   const int imax = std::min(i0 + MG, m);
   const int mlen = imax - i0;

   // -------------------------------------------------------------------------
   // Read marker block, compute/use allele frequencies, and build lookup tables
   // -------------------------------------------------------------------------
   for (int ii = 0; ii < mlen; ++ii) {
    const int i = i0 + ii;

    const long long offset =
     3LL + static_cast<long long>(cls[i] - 1) * static_cast<long long>(nbytes);

    if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
     throw std::runtime_error("fseek failed while reading BED marker " + std::to_string(cls[i]) + ".");
    }

    unsigned char* buf_i = block_buffer + static_cast<std::size_t>(ii) * nbytes;

    const std::size_t nbytes_read =
     std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);

    if (nbytes_read != nbytes) {
     throw std::runtime_error(
       "Error reading BED data: nbytes_read != nbytes for marker " +
        std::to_string(cls[i]) + "."
     );
    }

    double p = 0.0;

    if (use_af) {
     p = af[i];
    } else {
     p = compute_af_from_packed_marker(buf_i, rows0, n_used);
    }

    if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
     throw std::runtime_error(
       "Invalid allele frequency for selected marker index " + std::to_string(i) +
        ", BED marker " + std::to_string(cls[i]) + "."
     );
    }

    af_used[i] = p;

    if (scale) {
     const double denom = std::sqrt(2.0 * p * (1.0 - p));

     if (denom <= 0.0 || !std::isfinite(denom)) {
      map0[static_cast<std::size_t>(ii)] = 0.0;
      map1[static_cast<std::size_t>(ii)] = 0.0;
      map2[static_cast<std::size_t>(ii)] = 0.0;
      map3[static_cast<std::size_t>(ii)] = 0.0;
     } else {
      map0[static_cast<std::size_t>(ii)] = (2.0 - 2.0 * p) / denom;  // BED 00 -> dosage 2
      map1[static_cast<std::size_t>(ii)] = 0.0;                      // BED 01 -> missing
      map2[static_cast<std::size_t>(ii)] = (1.0 - 2.0 * p) / denom;  // BED 10 -> dosage 1
      map3[static_cast<std::size_t>(ii)] = (0.0 - 2.0 * p) / denom;  // BED 11 -> dosage 0
     }
    } else {
     // Raw dosage coding with missing mean-imputed to 2p.
     // This is intentionally uncentered.
     map0[static_cast<std::size_t>(ii)] = 2.0;
     map1[static_cast<std::size_t>(ii)] = 2.0 * p;
     map2[static_cast<std::size_t>(ii)] = 1.0;
     map3[static_cast<std::size_t>(ii)] = 0.0;
    }
   }

   // Block-level accumulators.
   std::vector<double> xty_block(static_cast<std::size_t>(mlen) * nt, 0.0);
   std::vector<double> xx_block(static_cast<std::size_t>(mlen), 0.0);

   // -------------------------------------------------------------------------
   // Parallel accumulation over selected individuals
   // -------------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 std::vector<double> xty_local(static_cast<std::size_t>(mlen) * nt, 0.0);
 std::vector<double> xx_local(static_cast<std::size_t>(mlen), 0.0);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
 for (int k0 = 0; k0 < n_used; k0 += JB) {
  const int kmax = std::min(k0 + JB, n_used);

  for (int ii = 0; ii < mlen; ++ii) {
   const unsigned char* buf_i =
    block_buffer + static_cast<std::size_t>(ii) * nbytes;

   const double m0 = map0[static_cast<std::size_t>(ii)];
   const double m1 = map1[static_cast<std::size_t>(ii)];
   const double m2 = map2[static_cast<std::size_t>(ii)];
   const double m3 = map3[static_cast<std::size_t>(ii)];

   double x2_sum = 0.0;
   double* xty_i = &xty_local[static_cast<std::size_t>(ii) * nt];

   if (!rows0) {
    // Fast path: all BED rows are used in their natural order.
    // Decode one packed byte at a time, as in the original implementation.
    const int byte0 = k0 >> 2;
    const int byte1 = (kmax + 3) >> 2;

    for (int kb = byte0; kb < byte1; ++kb) {
     unsigned char buf_k = buf_i[kb];
     const int jbase = kb << 2;

     for (int pos = 0; pos < 4; ++pos) {
      const int k = jbase + pos;

      if (k < k0 || k >= kmax || k >= n_used) {
       buf_k >>= 2;
       continue;
      }

      double xj;
      switch (buf_k & 3u) {
      case 0u: xj = m0; break;
      case 1u: xj = m1; break;
      case 2u: xj = m2; break;
      default: xj = m3; break;
      }

      buf_k >>= 2;

      x2_sum += xj * xj;

      const double* yk = &y_flat[static_cast<std::size_t>(k) * nt];

      for (int t0 = 0; t0 < nt; t0 += TB) {
       const int tmax = std::min(t0 + TB, nt);

#pragma omp simd
       for (int t = t0; t < tmax; ++t) {
        xty_i[t] += xj * yk[t];
       }
      }
     }
    }
   } else {
    // Generic path: selected BED rows can be arbitrary.
    for (int k = k0; k < kmax; ++k) {
     const int bed_row = rows0[k];
     const unsigned char byte = buf_i[bed_row >> 2];
     const unsigned int code = (byte >> (2 * (bed_row & 3))) & 3u;

     double xj;
     switch (code) {
     case 0u: xj = m0; break;
     case 1u: xj = m1; break;
     case 2u: xj = m2; break;
     default: xj = m3; break;
     }

     x2_sum += xj * xj;

     const double* yk = &y_flat[static_cast<std::size_t>(k) * nt];

     for (int t0 = 0; t0 < nt; t0 += TB) {
      const int tmax = std::min(t0 + TB, nt);

#pragma omp simd
      for (int t = t0; t < tmax; ++t) {
       xty_i[t] += xj * yk[t];
      }
     }
    }
   }

   xx_local[static_cast<std::size_t>(ii)] += x2_sum;
  }
 }

#ifdef _OPENMP
#pragma omp critical
#endif
{
 for (int ii = 0; ii < mlen; ++ii) {
  xx_block[static_cast<std::size_t>(ii)] += xx_local[static_cast<std::size_t>(ii)];

  double* dst = &xty_block[static_cast<std::size_t>(ii) * nt];
  const double* src = &xty_local[static_cast<std::size_t>(ii) * nt];

  for (int t = 0; t < nt; ++t) {
   dst[t] += src[t];
  }
 }
}
}

// -------------------------------------------------------------------------
// Copy block accumulators to output
// -------------------------------------------------------------------------
for (int ii = 0; ii < mlen; ++ii) {
 const int i = i0 + ii;

 xx[i] = xx_block[static_cast<std::size_t>(ii)];

 double* dst = &xty_flat[static_cast<std::size_t>(i) * nt];
 const double* src = &xty_block[static_cast<std::size_t>(ii) * nt];

 for (int t = 0; t < nt; ++t) {
  dst[t] = src[t];
 }
}
  }
 } catch (...) {
  std::free(block_buffer);
  std::fclose(file_stream);
  throw;
 }

 std::free(block_buffer);
 std::fclose(file_stream);
}

// =============================================================================
// Wrapper-neutral C++ convenience function
// =============================================================================

BedXtXytyResult bed_xtx_xty_cpp(
  const std::string& bed_file,
  int n_bed,
  const std::vector<int>& cls,
  const std::vector<double>& af,
  bool use_af,
  const std::vector<int>& rows0,
  bool use_rows,
  const std::vector<double>& y_flat,
  int nt,
  bool scale,
  int nthreads,
  int MG,
  int JB,
  int TB
) {
 if (bed_file.empty()) {
  throw std::runtime_error("bed_file must not be empty.");
 }

 if (n_bed <= 0) {
  throw std::runtime_error("n_bed must be positive.");
 }

 const int m = static_cast<int>(cls.size());
 const int n_used = use_rows ? static_cast<int>(rows0.size()) : n_bed;

 if (m <= 0) {
  throw std::runtime_error("cls must be non-empty.");
 }

 if (n_used <= 0) {
  throw std::runtime_error("selected rows must be non-empty.");
 }

 if (nt <= 0) {
  throw std::runtime_error("nt must be positive.");
 }

 if (use_af && static_cast<int>(af.size()) != m) {
  throw std::runtime_error("af must have same length as cls when use_af is true.");
 }

 const std::size_t expected_y_size =
  static_cast<std::size_t>(n_used) * static_cast<std::size_t>(nt);

 if (y_flat.size() != expected_y_size) {
  throw std::runtime_error("y_flat must have length n_used * nt.");
 }

 if (MG <= 0) {
  throw std::runtime_error("MG must be positive.");
 }

 if (JB <= 0) {
  throw std::runtime_error("JB must be positive.");
 }

 if (TB <= 0) {
  throw std::runtime_error("TB must be positive.");
 }

 for (int i = 0; i < m; ++i) {
  if (cls[static_cast<std::size_t>(i)] <= 0) {
   throw std::runtime_error("cls must contain positive 1-based marker indices.");
  }

  if (use_af) {
   const double p = af[static_cast<std::size_t>(i)];

   if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
    throw std::runtime_error("af must contain finite values in [0, 1].");
   }
  }
 }

 if (use_rows) {
  std::vector<unsigned char> seen(static_cast<std::size_t>(n_bed), 0);

  for (int k = 0; k < n_used; ++k) {
   const int r = rows0[static_cast<std::size_t>(k)];

   if (r < 0 || r >= n_bed) {
    throw std::runtime_error("rows0 contains an out-of-range BED row index.");
   }

   if (seen[static_cast<std::size_t>(r)]) {
    throw std::runtime_error("rows0 contains duplicate row indices.");
   }

   seen[static_cast<std::size_t>(r)] = 1;
  }
 }

 for (std::size_t k = 0; k < y_flat.size(); ++k) {
  if (!std::isfinite(y_flat[k])) {
   throw std::runtime_error("y_flat contains NaN/Inf.");
  }
 }

 std::vector<double> xty_flat(static_cast<std::size_t>(m) * nt, 0.0);
 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
 std::vector<double> af_used(static_cast<std::size_t>(m), 0.0);

 bed_xtx_diag_xty_core(
  bed_file.c_str(),
  n_bed,
  cls.data(),
  use_af ? af.data() : nullptr,
  use_af,
  af_used.data(),
  use_rows ? rows0.data() : nullptr,
  n_used,
  m,
  nt,
  scale,
  y_flat.data(),
  xty_flat.data(),
  xx.data(),
  nthreads,
  MG,
  JB,
  TB
 );

 BedXtXytyResult out;
 out.n = n_used;
 out.n_bed = n_bed;
 out.m = m;
 out.nt = nt;
 out.scale = scale;
 out.rows_used = use_rows;
 out.xx = xx;
 out.af = af_used;
 out.af_computed = !use_af;
 if (use_rows) out.rows0 = rows0;

 out.wy.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
 out.ww.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
 out.yy.assign(static_cast<std::size_t>(nt), 0.0);

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   out.wy[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
    xty_flat[static_cast<std::size_t>(i) * nt + t];

   out.ww[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
    xx[static_cast<std::size_t>(i)];
  }
 }

 for (int t = 0; t < nt; ++t) {
  double s = 0.0;

  for (int k = 0; k < n_used; ++k) {
   const double v = y_flat[static_cast<std::size_t>(k) * nt + t];
   s += v * v;
  }

  out.yy[static_cast<std::size_t>(t)] = s;
 }

 out.min_xx = *std::min_element(xx.begin(), xx.end());
 out.max_xx = *std::max_element(xx.begin(), xx.end());
 out.mean_xx = std::accumulate(xx.begin(), xx.end(), 0.0) /
  static_cast<double>(xx.size());

 return out;
}

// =============================================================================
// Thin Rcpp adapter
// =============================================================================

// [[Rcpp::export]]
Rcpp::List bed_xtx_xty(
  std::string bed_file,
  int n,
  Rcpp::IntegerVector cls,
  Rcpp::Nullable<Rcpp::NumericVector> af,
  Rcpp::NumericMatrix y,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  bool scale = true,
  int nthreads = 1,
  int MG = 64,
  int JB = 1024,
  int TB = 32
) {
 const int n_bed = n;
 const int m = cls.size();
 const int nt = y.ncol();

 std::vector<int> cls_cpp(static_cast<std::size_t>(m));
 std::vector<double> af_cpp;
 std::vector<int> rows0;

 for (int i = 0; i < m; ++i) {
  if (cls[i] == NA_INTEGER) {
   Rcpp::stop("cls contains NA.");
  }

  cls_cpp[static_cast<std::size_t>(i)] = cls[i];
 }

 const bool use_af = af.isNotNull();

 if (use_af) {
  Rcpp::NumericVector af_r = Rcpp::as<Rcpp::NumericVector>(af.get());

  if (af_r.size() != m) {
   Rcpp::stop("af must have same length as cls when supplied.");
  }

  af_cpp.resize(static_cast<std::size_t>(m));

  for (int i = 0; i < m; ++i) {
   af_cpp[static_cast<std::size_t>(i)] = af_r[i];
  }
 }

 bool use_rows = false;

 if (rows.isNotNull()) {
  Rcpp::IntegerVector rows_r = Rcpp::as<Rcpp::IntegerVector>(rows.get());

  if (rows_r.size() <= 0) {
   Rcpp::stop("rows must be non-empty when supplied.");
  }

  rows0.resize(static_cast<std::size_t>(rows_r.size()));

  bool identity_rows = (rows_r.size() == n_bed);

  for (int k = 0; k < rows_r.size(); ++k) {
   if (rows_r[k] == NA_INTEGER) {
    Rcpp::stop("rows contains NA.");
   }

   rows0[static_cast<std::size_t>(k)] = rows_r[k] - 1; // R 1-based -> C++ 0-based

   if (rows_r[k] != k + 1) {
    identity_rows = false;
   }
  }

  if (identity_rows) {
   rows0.clear();
   use_rows = false;
  } else {
   use_rows = true;
  }
 }

 const int n_used = use_rows ? static_cast<int>(rows0.size()) : n_bed;

 if (y.nrow() != n_used) {
  Rcpp::stop("y must have nrow equal to length(rows) when rows is supplied, otherwise n.");
 }

 std::vector<double> y_flat(static_cast<std::size_t>(n_used) * nt);

 for (int k = 0; k < n_used; ++k) {
  for (int t = 0; t < nt; ++t) {
   y_flat[static_cast<std::size_t>(k) * nt + t] = y(k, t);
  }
 }

 BedXtXytyResult res = bed_xtx_xty_cpp(
  bed_file,
  n_bed,
  cls_cpp,
  af_cpp,
  use_af,
  rows0,
  use_rows,
  y_flat,
  nt,
  scale,
  nthreads,
  MG,
  JB,
  TB
 );

 Rcpp::IntegerVector rows_out;
 if (res.rows_used) {
  rows_out = Rcpp::IntegerVector(res.rows0.size());
  for (std::size_t k = 0; k < res.rows0.size(); ++k) {
   rows_out[static_cast<int>(k)] = res.rows0[k] + 1;
  }
 }

 return Rcpp::List::create(
  Rcpp::Named("wy") = res.wy,
  Rcpp::Named("ww") = res.ww,
  Rcpp::Named("yy") = res.yy,
  Rcpp::Named("xx") = res.xx,
  Rcpp::Named("af") = res.af,
  Rcpp::Named("af_computed") = res.af_computed,
  Rcpp::Named("rows") = rows_out,
  Rcpp::Named("rows_used") = res.rows_used,
  Rcpp::Named("n") = res.n,
  Rcpp::Named("n_bed") = res.n_bed,
  Rcpp::Named("m") = res.m,
  Rcpp::Named("nt") = res.nt,
  Rcpp::Named("scale") = res.scale,
  Rcpp::Named("mean_xx") = res.mean_xx,
  Rcpp::Named("min_xx") = res.min_xx,
  Rcpp::Named("max_xx") = res.max_xx
 );
}


// // =============================================================================
// // mt_sufficient_statistics.cpp
// // =============================================================================
//
// #include "mt_sufficient_statistics_core.h"
//
// #include <Rcpp.h>
//
// #include <algorithm>
// #include <cmath>
// #include <cstdio>
// #include <cstdlib>
// #include <cstring>
// #include <limits>
// #include <numeric>
// #include <stdexcept>
// #include <string>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// // =============================================================================
// // BED helpers
// // =============================================================================
// //
// // BED coding assumed:
// //   00 -> dosage 2
// //   01 -> missing
// //   10 -> dosage 1
// //   11 -> dosage 0
// //
// // =============================================================================
//
// inline double decode_bed_dosage_from_marker(
//   const unsigned char* buf,
//   int bed_row
// ) {
//  const unsigned char x = buf[bed_row >> 2];
//  const unsigned int code = (x >> (2 * (bed_row & 3))) & 3u;
//
//  if (code == 0u) return 2.0;
//  if (code == 2u) return 1.0;
//  if (code == 3u) return 0.0;
//  return std::numeric_limits<double>::quiet_NaN(); // missing
// }
//
// inline double compute_af_from_packed_marker(
//   const unsigned char* buf,
//   const int* rows0,
//   int n_used
// ) {
//  double dosage_sum = 0.0;
//  int n_obs = 0;
//
//  if (!rows0) {
//   const int nbytes = (n_used + 3) / 4;
//
//   for (int kb = 0; kb < nbytes; ++kb) {
//    unsigned char x = buf[kb];
//    const int jbase = kb << 2;
//
//    for (int pos = 0; pos < 4; ++pos) {
//     const int j = jbase + pos;
//
//     if (j >= n_used) {
//      break;
//     }
//
//     const unsigned int code = x & 3u;
//     x >>= 2;
//
//     if (code == 1u) {
//      continue; // missing
//     }
//
//     if (code == 0u) {
//      dosage_sum += 2.0;
//     } else if (code == 2u) {
//      dosage_sum += 1.0;
//     }
//
//     ++n_obs;
//    }
//   }
//  } else {
//   for (int k = 0; k < n_used; ++k) {
//    const int bed_row = rows0[k];
//    const double dosage = decode_bed_dosage_from_marker(buf, bed_row);
//
//    if (!std::isfinite(dosage)) {
//     continue;
//    }
//
//    dosage_sum += dosage;
//    ++n_obs;
//   }
//  }
//
//  if (n_obs <= 0) {
//   return 0.0;
//  }
//
//  return dosage_sum / (2.0 * static_cast<double>(n_obs));
// }
//
// // =============================================================================
// // BED sufficient-statistics core
// // =============================================================================
// //
// // If scale = true:
// //   x = (dosage - 2p) / sqrt(2p(1-p))
// //   missing -> 0 after standardization, i.e. mean-imputed standardized value.
// //
// // If scale = false:
// //   x = raw dosage, missing -> 2p.
// //   This is uncentered raw dosage coding.
// //
// // y_flat is n_used x nt, row-major: y_flat[k * nt + t], where k indexes
// // selected individuals, not necessarily BED row number.
// //
// // cls is 1-based marker index into the BED file.
// // rows0 is 0-based individual index into the BED file. If rows0 is null, all
// // rows 0:(n_bed - 1) are used in order.
// //
// // If use_af = false, allele frequencies are computed over the selected rows
// // and stored in af_used.
// //
// // =============================================================================
//
// void bed_xtx_diag_xty_core(
//   const char* file,
//   int n_bed,
//   const int* cls,
//   const double* af,
//   bool use_af,
//   double* af_used,
//   const int* rows0,
//   int n_used,
//   int m,
//   int nt,
//   bool scale,
//   const double* y_flat,
//   double* xty_flat,
//   double* xx,
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  if (!file) {
//   throw std::runtime_error("file path pointer must not be null.");
//  }
//
//  FILE* file_stream = std::fopen(file, "rb");
//  if (!file_stream) {
//   throw std::runtime_error("Could not open BED file.");
//  }
//
//  if (n_bed <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("n_bed must be positive.");
//  }
//
//  if (n_used <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("n_used must be positive.");
//  }
//
//  if (m <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("m must be positive.");
//  }
//
//  if (nt <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (!af_used) {
//   std::fclose(file_stream);
//   throw std::runtime_error("af_used must not be null.");
//  }
//
//  if (use_af && !af) {
//   std::fclose(file_stream);
//   throw std::runtime_error("af must not be null when use_af is true.");
//  }
//
//  if (MG <= 0) MG = 64;
//  if (JB <= 0) JB = 1024;
//  if (TB <= 0) TB = 32;
//
//  // Validate PLINK BED header: magic bytes + SNP-major mode.
//  unsigned char magic[3];
//  if (std::fread(magic, sizeof(unsigned char), 3, file_stream) != 3) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Could not read BED header.");
//  }
//
//  if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
//   std::fclose(file_stream);
//   throw std::runtime_error(
//     "BED file must be PLINK SNP-major format: expected magic bytes 0x6c 0x1b 0x01."
//   );
//  }
//
//  const std::size_t nbytes = static_cast<std::size_t>((n_bed + 3) / 4);
//
//  std::memset(xty_flat, 0, static_cast<std::size_t>(m) * nt * sizeof(double));
//  std::memset(xx, 0, static_cast<std::size_t>(m) * sizeof(double));
//  std::memset(af_used, 0, static_cast<std::size_t>(m) * sizeof(double));
//
//  unsigned char* block_buffer =
//   static_cast<unsigned char*>(std::malloc(static_cast<std::size_t>(MG) * nbytes));
//
//  if (!block_buffer) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Failed to allocate BED block buffer.");
//  }
//
//  std::vector<double> map0(static_cast<std::size_t>(MG));
//  std::vector<double> map1(static_cast<std::size_t>(MG));
//  std::vector<double> map2(static_cast<std::size_t>(MG));
//  std::vector<double> map3(static_cast<std::size_t>(MG));
//
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  try {
//   for (int i0 = 0; i0 < m; i0 += MG) {
//    const int imax = std::min(i0 + MG, m);
//    const int mlen = imax - i0;
//
//    // -------------------------------------------------------------------------
//    // Read marker block, compute/use allele frequencies, and build lookup tables
//    // -------------------------------------------------------------------------
//    for (int ii = 0; ii < mlen; ++ii) {
//     const int i = i0 + ii;
//
//     const long long offset =
//      3LL + static_cast<long long>(cls[i] - 1) * static_cast<long long>(nbytes);
//
//     if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
//      throw std::runtime_error("fseek failed while reading BED marker " + std::to_string(cls[i]) + ".");
//     }
//
//     unsigned char* buf_i = block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//     const std::size_t nbytes_read =
//      std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);
//
//     if (nbytes_read != nbytes) {
//      throw std::runtime_error(
//        "Error reading BED data: nbytes_read != nbytes for marker " +
//         std::to_string(cls[i]) + "."
//      );
//     }
//
//     double p = 0.0;
//
//     if (use_af) {
//      p = af[i];
//     } else {
//      p = compute_af_from_packed_marker(buf_i, rows0, n_used);
//     }
//
//     if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//      throw std::runtime_error(
//        "Invalid allele frequency for selected marker index " + std::to_string(i) +
//         ", BED marker " + std::to_string(cls[i]) + "."
//      );
//     }
//
//     af_used[i] = p;
//
//     if (scale) {
//      const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//      if (denom <= 0.0 || !std::isfinite(denom)) {
//       map0[static_cast<std::size_t>(ii)] = 0.0;
//       map1[static_cast<std::size_t>(ii)] = 0.0;
//       map2[static_cast<std::size_t>(ii)] = 0.0;
//       map3[static_cast<std::size_t>(ii)] = 0.0;
//      } else {
//       map0[static_cast<std::size_t>(ii)] = (2.0 - 2.0 * p) / denom;  // BED 00 -> dosage 2
//       map1[static_cast<std::size_t>(ii)] = 0.0;                      // BED 01 -> missing
//       map2[static_cast<std::size_t>(ii)] = (1.0 - 2.0 * p) / denom;  // BED 10 -> dosage 1
//       map3[static_cast<std::size_t>(ii)] = (0.0 - 2.0 * p) / denom;  // BED 11 -> dosage 0
//      }
//     } else {
//      // Raw dosage coding with missing mean-imputed to 2p.
//      // This is intentionally uncentered.
//      map0[static_cast<std::size_t>(ii)] = 2.0;
//      map1[static_cast<std::size_t>(ii)] = 2.0 * p;
//      map2[static_cast<std::size_t>(ii)] = 1.0;
//      map3[static_cast<std::size_t>(ii)] = 0.0;
//     }
//    }
//
//    // Block-level accumulators.
//    std::vector<double> xty_block(static_cast<std::size_t>(mlen) * nt, 0.0);
//    std::vector<double> xx_block(static_cast<std::size_t>(mlen), 0.0);
//
//    // -------------------------------------------------------------------------
//    // Parallel accumulation over selected individuals
//    // -------------------------------------------------------------------------
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// {
//  std::vector<double> xty_local(static_cast<std::size_t>(mlen) * nt, 0.0);
//  std::vector<double> xx_local(static_cast<std::size_t>(mlen), 0.0);
//
// #ifdef _OPENMP
// #pragma omp for schedule(static)
// #endif
//  for (int k0 = 0; k0 < n_used; k0 += JB) {
//   const int kmax = std::min(k0 + JB, n_used);
//
//   for (int ii = 0; ii < mlen; ++ii) {
//    const unsigned char* buf_i =
//     block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//    const double m0 = map0[static_cast<std::size_t>(ii)];
//    const double m1 = map1[static_cast<std::size_t>(ii)];
//    const double m2 = map2[static_cast<std::size_t>(ii)];
//    const double m3 = map3[static_cast<std::size_t>(ii)];
//
//    double x2_sum = 0.0;
//    double* xty_i = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//    if (!rows0) {
//     // Fast path: all BED rows are used in their natural order.
//     // Decode one packed byte at a time, as in the original implementation.
//     const int byte0 = k0 >> 2;
//     const int byte1 = (kmax + 3) >> 2;
//
//     for (int kb = byte0; kb < byte1; ++kb) {
//      unsigned char buf_k = buf_i[kb];
//      const int jbase = kb << 2;
//
//      for (int pos = 0; pos < 4; ++pos) {
//       const int k = jbase + pos;
//
//       if (k < k0 || k >= kmax || k >= n_used) {
//        buf_k >>= 2;
//        continue;
//       }
//
//       double xj;
//       switch (buf_k & 3u) {
//       case 0u: xj = m0; break;
//       case 1u: xj = m1; break;
//       case 2u: xj = m2; break;
//       default: xj = m3; break;
//       }
//
//       buf_k >>= 2;
//
//       x2_sum += xj * xj;
//
//       const double* yk = &y_flat[static_cast<std::size_t>(k) * nt];
//
//       for (int t0 = 0; t0 < nt; t0 += TB) {
//        const int tmax = std::min(t0 + TB, nt);
//
// #pragma omp simd
//        for (int t = t0; t < tmax; ++t) {
//         xty_i[t] += xj * yk[t];
//        }
//       }
//      }
//     }
//    } else {
//     // Generic path: selected BED rows can be arbitrary.
//     for (int k = k0; k < kmax; ++k) {
//      const int bed_row = rows0[k];
//      const unsigned char byte = buf_i[bed_row >> 2];
//      const unsigned int code = (byte >> (2 * (bed_row & 3))) & 3u;
//
//      double xj;
//      switch (code) {
//      case 0u: xj = m0; break;
//      case 1u: xj = m1; break;
//      case 2u: xj = m2; break;
//      default: xj = m3; break;
//      }
//
//      x2_sum += xj * xj;
//
//      const double* yk = &y_flat[static_cast<std::size_t>(k) * nt];
//
//      for (int t0 = 0; t0 < nt; t0 += TB) {
//       const int tmax = std::min(t0 + TB, nt);
//
// #pragma omp simd
//       for (int t = t0; t < tmax; ++t) {
//        xty_i[t] += xj * yk[t];
//       }
//      }
//     }
//    }
//
//    xx_local[static_cast<std::size_t>(ii)] += x2_sum;
//   }
//  }
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//  for (int ii = 0; ii < mlen; ++ii) {
//   xx_block[static_cast<std::size_t>(ii)] += xx_local[static_cast<std::size_t>(ii)];
//
//   double* dst = &xty_block[static_cast<std::size_t>(ii) * nt];
//   const double* src = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//   for (int t = 0; t < nt; ++t) {
//    dst[t] += src[t];
//   }
//  }
// }
// }
//
// // -------------------------------------------------------------------------
// // Copy block accumulators to output
// // -------------------------------------------------------------------------
// for (int ii = 0; ii < mlen; ++ii) {
//  const int i = i0 + ii;
//
//  xx[i] = xx_block[static_cast<std::size_t>(ii)];
//
//  double* dst = &xty_flat[static_cast<std::size_t>(i) * nt];
//  const double* src = &xty_block[static_cast<std::size_t>(ii) * nt];
//
//  for (int t = 0; t < nt; ++t) {
//   dst[t] = src[t];
//  }
// }
//   }
//  } catch (...) {
//   std::free(block_buffer);
//   std::fclose(file_stream);
//   throw;
//  }
//
//  std::free(block_buffer);
//  std::fclose(file_stream);
// }
//
// // =============================================================================
// // Wrapper-neutral C++ convenience function
// // =============================================================================
//
// BedXtXytyResult bed_xtx_xty_cpp(
//   const std::string& bed_file,
//   int n_bed,
//   const std::vector<int>& cls,
//   const std::vector<double>& af,
//   bool use_af,
//   const std::vector<int>& rows0,
//   bool use_rows,
//   const std::vector<double>& y_flat,
//   int nt,
//   bool scale,
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  if (bed_file.empty()) {
//   throw std::runtime_error("bed_file must not be empty.");
//  }
//
//  if (n_bed <= 0) {
//   throw std::runtime_error("n_bed must be positive.");
//  }
//
//  const int m = static_cast<int>(cls.size());
//  const int n_used = use_rows ? static_cast<int>(rows0.size()) : n_bed;
//
//  if (m <= 0) {
//   throw std::runtime_error("cls must be non-empty.");
//  }
//
//  if (n_used <= 0) {
//   throw std::runtime_error("selected rows must be non-empty.");
//  }
//
//  if (nt <= 0) {
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (use_af && static_cast<int>(af.size()) != m) {
//   throw std::runtime_error("af must have same length as cls when use_af is true.");
//  }
//
//  const std::size_t expected_y_size =
//   static_cast<std::size_t>(n_used) * static_cast<std::size_t>(nt);
//
//  if (y_flat.size() != expected_y_size) {
//   throw std::runtime_error("y_flat must have length n_used * nt.");
//  }
//
//  if (MG <= 0) {
//   throw std::runtime_error("MG must be positive.");
//  }
//
//  if (JB <= 0) {
//   throw std::runtime_error("JB must be positive.");
//  }
//
//  if (TB <= 0) {
//   throw std::runtime_error("TB must be positive.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   if (cls[static_cast<std::size_t>(i)] <= 0) {
//    throw std::runtime_error("cls must contain positive 1-based marker indices.");
//   }
//
//   if (use_af) {
//    const double p = af[static_cast<std::size_t>(i)];
//
//    if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//     throw std::runtime_error("af must contain finite values in [0, 1].");
//    }
//   }
//  }
//
//  if (use_rows) {
//   std::vector<unsigned char> seen(static_cast<std::size_t>(n_bed), 0);
//
//   for (int k = 0; k < n_used; ++k) {
//    const int r = rows0[static_cast<std::size_t>(k)];
//
//    if (r < 0 || r >= n_bed) {
//     throw std::runtime_error("rows0 contains an out-of-range BED row index.");
//    }
//
//    if (seen[static_cast<std::size_t>(r)]) {
//     throw std::runtime_error("rows0 contains duplicate row indices.");
//    }
//
//    seen[static_cast<std::size_t>(r)] = 1;
//   }
//  }
//
//  for (std::size_t k = 0; k < y_flat.size(); ++k) {
//   if (!std::isfinite(y_flat[k])) {
//    throw std::runtime_error("y_flat contains NaN/Inf.");
//   }
//  }
//
//  std::vector<double> xty_flat(static_cast<std::size_t>(m) * nt, 0.0);
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//  std::vector<double> af_used(static_cast<std::size_t>(m), 0.0);
//
//  bed_xtx_diag_xty_core(
//   bed_file.c_str(),
//   n_bed,
//   cls.data(),
//   use_af ? af.data() : nullptr,
//   use_af,
//   af_used.data(),
//   use_rows ? rows0.data() : nullptr,
//   n_used,
//   m,
//   nt,
//   scale,
//   y_flat.data(),
//   xty_flat.data(),
//   xx.data(),
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  BedXtXytyResult out;
//  out.n = n_used;
//  out.n_bed = n_bed;
//  out.m = m;
//  out.nt = nt;
//  out.scale = scale;
//  out.rows_used = use_rows;
//  out.xx = xx;
//  out.af = af_used;
//  out.af_computed = !use_af;
//  if (use_rows) out.rows0 = rows0;
//
//  out.wy.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
//  out.ww.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
//  out.yy.assign(static_cast<std::size_t>(nt), 0.0);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    out.wy[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
//     xty_flat[static_cast<std::size_t>(i) * nt + t];
//
//    out.ww[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
//     xx[static_cast<std::size_t>(i)];
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   double s = 0.0;
//
//   for (int k = 0; k < n_used; ++k) {
//    const double v = y_flat[static_cast<std::size_t>(k) * nt + t];
//    s += v * v;
//   }
//
//   out.yy[static_cast<std::size_t>(t)] = s;
//  }
//
//  out.min_xx = *std::min_element(xx.begin(), xx.end());
//  out.max_xx = *std::max_element(xx.begin(), xx.end());
//  out.mean_xx = std::accumulate(xx.begin(), xx.end(), 0.0) /
//   static_cast<double>(xx.size());
//
//  return out;
// }
//
// // =============================================================================
// // Thin Rcpp adapter
// // =============================================================================
//
// // [[Rcpp::export]]
// Rcpp::List bed_xtx_xty(
//   std::string bed_file,
//   int n,
//   Rcpp::IntegerVector cls,
//   Rcpp::Nullable<Rcpp::NumericVector> af,
//   Rcpp::NumericMatrix y,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   bool scale = true,
//   int nthreads = 1,
//   int MG = 64,
//   int JB = 1024,
//   int TB = 32
// ) {
//  const int n_bed = n;
//  const int m = cls.size();
//  const int nt = y.ncol();
//
//  std::vector<int> cls_cpp(static_cast<std::size_t>(m));
//  std::vector<double> af_cpp;
//  std::vector<int> rows0;
//
//  for (int i = 0; i < m; ++i) {
//   if (cls[i] == NA_INTEGER) {
//    Rcpp::stop("cls contains NA.");
//   }
//
//   cls_cpp[static_cast<std::size_t>(i)] = cls[i];
//  }
//
//  const bool use_af = af.isNotNull();
//
//  if (use_af) {
//   Rcpp::NumericVector af_r = Rcpp::as<Rcpp::NumericVector>(af.get());
//
//   if (af_r.size() != m) {
//    Rcpp::stop("af must have same length as cls when supplied.");
//   }
//
//   af_cpp.resize(static_cast<std::size_t>(m));
//
//   for (int i = 0; i < m; ++i) {
//    af_cpp[static_cast<std::size_t>(i)] = af_r[i];
//   }
//  }
//
//  bool use_rows = false;
//
//  if (rows.isNotNull()) {
//   Rcpp::IntegerVector rows_r = Rcpp::as<Rcpp::IntegerVector>(rows.get());
//
//   if (rows_r.size() <= 0) {
//    Rcpp::stop("rows must be non-empty when supplied.");
//   }
//
//   rows0.resize(static_cast<std::size_t>(rows_r.size()));
//
//   bool identity_rows = (rows_r.size() == n_bed);
//
//   for (int k = 0; k < rows_r.size(); ++k) {
//    if (rows_r[k] == NA_INTEGER) {
//     Rcpp::stop("rows contains NA.");
//    }
//
//    rows0[static_cast<std::size_t>(k)] = rows_r[k] - 1; // R 1-based -> C++ 0-based
//
//    if (rows_r[k] != k + 1) {
//     identity_rows = false;
//    }
//   }
//
//   if (identity_rows) {
//    rows0.clear();
//    use_rows = false;
//   } else {
//    use_rows = true;
//   }
//  }
//
//  const int n_used = use_rows ? static_cast<int>(rows0.size()) : n_bed;
//
//  if (y.nrow() != n_used) {
//   Rcpp::stop("y must have nrow equal to length(rows) when rows is supplied, otherwise n.");
//  }
//
//  std::vector<double> y_flat(static_cast<std::size_t>(n_used) * nt);
//
//  for (int k = 0; k < n_used; ++k) {
//   for (int t = 0; t < nt; ++t) {
//    y_flat[static_cast<std::size_t>(k) * nt + t] = y(k, t);
//   }
//  }
//
//  BedXtXytyResult res = bed_xtx_xty_cpp(
//   bed_file,
//   n_bed,
//   cls_cpp,
//   af_cpp,
//   use_af,
//   rows0,
//   use_rows,
//   y_flat,
//   nt,
//   scale,
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  Rcpp::IntegerVector rows_out;
//  if (res.rows_used) {
//   rows_out = Rcpp::IntegerVector(res.rows0.size());
//   for (std::size_t k = 0; k < res.rows0.size(); ++k) {
//    rows_out[static_cast<int>(k)] = res.rows0[k] + 1;
//   }
//  }
//
//  return Rcpp::List::create(
//   Rcpp::Named("wy") = res.wy,
//   Rcpp::Named("ww") = res.ww,
//   Rcpp::Named("yy") = res.yy,
//   Rcpp::Named("xx") = res.xx,
//   Rcpp::Named("af") = res.af,
//   Rcpp::Named("af_computed") = res.af_computed,
//   Rcpp::Named("rows") = rows_out,
//   Rcpp::Named("rows_used") = res.rows_used,
//   Rcpp::Named("n") = res.n,
//   Rcpp::Named("n_bed") = res.n_bed,
//   Rcpp::Named("m") = res.m,
//   Rcpp::Named("nt") = res.nt,
//   Rcpp::Named("scale") = res.scale,
//   Rcpp::Named("mean_xx") = res.mean_xx,
//   Rcpp::Named("min_xx") = res.min_xx,
//   Rcpp::Named("max_xx") = res.max_xx
//  );
// }



// // =============================================================================
// // mt_sufficient_statistics.cpp
// // =============================================================================
//
// #include "mt_sufficient_statistics_core.h"
//
// #include <Rcpp.h>
//
// #include <algorithm>
// #include <cmath>
// #include <cstdio>
// #include <cstdlib>
// #include <cstring>
// #include <numeric>
// #include <stdexcept>
// #include <string>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// // =============================================================================
// // BED helpers
// // =============================================================================
// //
// // BED coding assumed:
// //   00 -> dosage 2
// //   01 -> missing
// //   10 -> dosage 1
// //   11 -> dosage 0
// //
// // =============================================================================
//
// // The returned frequency is mean(dosage) / 2 over non-missing genotypes.
// inline double compute_af_from_packed_marker(
//   const unsigned char* buf,
//   int n
// ) {
//  double dosage_sum = 0.0;
//  int n_obs = 0;
//  const int nbytes = (n + 3) / 4;
//
//  for (int kb = 0; kb < nbytes; ++kb) {
//   unsigned char x = buf[kb];
//   const int jbase = kb << 2;
//
//   for (int pos = 0; pos < 4; ++pos) {
//    const int j = jbase + pos;
//
//    if (j >= n) {
//     break;
//    }
//
//    const unsigned int code = x & 3u;
//    x >>= 2;
//
//    if (code == 1u) {
//     continue; // missing
//    }
//
//    double dosage = 0.0;
//
//    if (code == 0u) {
//     dosage = 2.0;
//    } else if (code == 2u) {
//     dosage = 1.0;
//    } else {
//     dosage = 0.0;
//    }
//
//    dosage_sum += dosage;
//    ++n_obs;
//   }
//  }
//
//  if (n_obs <= 0) {
//   return 0.0;
//  }
//
//  return dosage_sum / (2.0 * static_cast<double>(n_obs));
// }
//
// // =============================================================================
// // BED sufficient-statistics core
// // =============================================================================
// //
// // If scale = true:
// //   x = (dosage - 2p) / sqrt(2p(1-p))
// //   missing -> 0 after standardization, i.e. mean-imputed standardized value.
// //
// // If scale = false:
// //   x = raw dosage, missing -> 2p.
// //   This is uncentered raw dosage coding.
// //
// // y_flat is n x nt, row-major: y_flat[j * nt + t].
// // xty_flat is m x nt, row-major: xty_flat[i * nt + t].
// // cls is 1-based marker index into the BED file.
// //
// // If use_af = false, allele frequencies are computed from the selected marker
// // bytes on the fly and stored in af_used.
// //
// // =============================================================================
//
// void bed_xtx_diag_xty_core(
//   const char* file,
//   int n,
//   const int* cls,
//   const double* af,
//   bool use_af,
//   double* af_used,
//   int m,
//   int nt,
//   bool scale,
//   const double* y_flat,
//   double* xty_flat,
//   double* xx,
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  if (!file) {
//   throw std::runtime_error("file path pointer must not be null.");
//  }
//
//  FILE* file_stream = std::fopen(file, "rb");
//  if (!file_stream) {
//   throw std::runtime_error("Could not open BED file.");
//  }
//
//  if (n <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("n must be positive.");
//  }
//
//  if (m <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("m must be positive.");
//  }
//
//  if (nt <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (!af_used) {
//   std::fclose(file_stream);
//   throw std::runtime_error("af_used must not be null.");
//  }
//
//  if (use_af && !af) {
//   std::fclose(file_stream);
//   throw std::runtime_error("af must not be null when use_af is true.");
//  }
//
//  if (MG <= 0) MG = 64;
//  if (JB <= 0) JB = 1024;
//  if (TB <= 0) TB = 32;
//
//  // Validate PLINK BED header: magic bytes + SNP-major mode.
//  unsigned char magic[3];
//  if (std::fread(magic, sizeof(unsigned char), 3, file_stream) != 3) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Could not read BED header.");
//  }
//
//  if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
//   std::fclose(file_stream);
//   throw std::runtime_error(
//     "BED file must be PLINK SNP-major format: expected magic bytes 0x6c 0x1b 0x01."
//   );
//  }
//
//  const std::size_t nbytes = static_cast<std::size_t>((n + 3) / 4);
//
//  std::memset(xty_flat, 0, static_cast<std::size_t>(m) * nt * sizeof(double));
//  std::memset(xx, 0, static_cast<std::size_t>(m) * sizeof(double));
//  std::memset(af_used, 0, static_cast<std::size_t>(m) * sizeof(double));
//
//  unsigned char* block_buffer =
//   static_cast<unsigned char*>(std::malloc(static_cast<std::size_t>(MG) * nbytes));
//
//  if (!block_buffer) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Failed to allocate BED block buffer.");
//  }
//
//  std::vector<double> map0(static_cast<std::size_t>(MG));
//  std::vector<double> map1(static_cast<std::size_t>(MG));
//  std::vector<double> map2(static_cast<std::size_t>(MG));
//  std::vector<double> map3(static_cast<std::size_t>(MG));
//
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  try {
//   for (int i0 = 0; i0 < m; i0 += MG) {
//    const int imax = std::min(i0 + MG, m);
//    const int mlen = imax - i0;
//
//    // -------------------------------------------------------------------------
//    // Read marker block, compute/use allele frequencies, and build lookup tables
//    // -------------------------------------------------------------------------
//    for (int ii = 0; ii < mlen; ++ii) {
//     const int i = i0 + ii;
//
//     const long long offset =
//      3LL + static_cast<long long>(cls[i] - 1) * static_cast<long long>(nbytes);
//
//     if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
//      throw std::runtime_error("fseek failed while reading BED marker " + std::to_string(cls[i]) + ".");
//     }
//
//     unsigned char* buf_i = block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//     const std::size_t nbytes_read =
//      std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);
//
//     if (nbytes_read != nbytes) {
//      throw std::runtime_error(
//        "Error reading BED data: nbytes_read != nbytes for marker " +
//         std::to_string(cls[i]) + "."
//      );
//     }
//
//     double p = 0.0;
//
//     if (use_af) {
//      p = af[i];
//     } else {
//      p = compute_af_from_packed_marker(buf_i, n);
//     }
//
//     if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//      throw std::runtime_error(
//        "Invalid allele frequency for selected marker index " + std::to_string(i) +
//         ", BED marker " + std::to_string(cls[i]) + "."
//      );
//     }
//
//     af_used[i] = p;
//
//     if (scale) {
//      const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//      if (denom <= 0.0 || !std::isfinite(denom)) {
//       map0[static_cast<std::size_t>(ii)] = 0.0;
//       map1[static_cast<std::size_t>(ii)] = 0.0;
//       map2[static_cast<std::size_t>(ii)] = 0.0;
//       map3[static_cast<std::size_t>(ii)] = 0.0;
//      } else {
//       map0[static_cast<std::size_t>(ii)] = (2.0 - 2.0 * p) / denom;  // BED 00 -> dosage 2
//       map1[static_cast<std::size_t>(ii)] = 0.0;                      // BED 01 -> missing
//       map2[static_cast<std::size_t>(ii)] = (1.0 - 2.0 * p) / denom;  // BED 10 -> dosage 1
//       map3[static_cast<std::size_t>(ii)] = (0.0 - 2.0 * p) / denom;  // BED 11 -> dosage 0
//      }
//     } else {
//      // Raw dosage coding with missing mean-imputed to 2p.
//      // This is intentionally uncentered.
//      map0[static_cast<std::size_t>(ii)] = 2.0;
//      map1[static_cast<std::size_t>(ii)] = 2.0 * p;
//      map2[static_cast<std::size_t>(ii)] = 1.0;
//      map3[static_cast<std::size_t>(ii)] = 0.0;
//     }
//    }
//
//    // Block-level accumulators.
//    std::vector<double> xty_block(static_cast<std::size_t>(mlen) * nt, 0.0);
//    std::vector<double> xx_block(static_cast<std::size_t>(mlen), 0.0);
//
//    // -------------------------------------------------------------------------
//    // Parallel accumulation over individual blocks
//    // -------------------------------------------------------------------------
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// {
//  std::vector<double> xty_local(static_cast<std::size_t>(mlen) * nt, 0.0);
//  std::vector<double> xx_local(static_cast<std::size_t>(mlen), 0.0);
//
// #ifdef _OPENMP
// #pragma omp for schedule(static)
// #endif
//  for (int j0 = 0; j0 < n; j0 += JB) {
//   const int jmax = std::min(j0 + JB, n);
//
//   const int byte0 = j0 >> 2;
//   const int byte1 = (jmax + 3) >> 2;
//
//   for (int ii = 0; ii < mlen; ++ii) {
//    const unsigned char* buf_i =
//     block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//    const double m0 = map0[static_cast<std::size_t>(ii)];
//    const double m1 = map1[static_cast<std::size_t>(ii)];
//    const double m2 = map2[static_cast<std::size_t>(ii)];
//    const double m3 = map3[static_cast<std::size_t>(ii)];
//
//    double x2_sum = 0.0;
//    double* xty_i = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//    for (int kb = byte0; kb < byte1; ++kb) {
//     unsigned char buf_k = buf_i[kb];
//     const int jbase = kb << 2;
//
//     for (int pos = 0; pos < 4; ++pos) {
//      const int j = jbase + pos;
//
//      if (j < j0 || j >= jmax || j >= n) {
//       buf_k >>= 2;
//       continue;
//      }
//
//      double xj;
//      switch (buf_k & 3u) {
//      case 0u: xj = m0; break;
//      case 1u: xj = m1; break;
//      case 2u: xj = m2; break;
//      default: xj = m3; break;
//      }
//
//      buf_k >>= 2;
//
//      x2_sum += xj * xj;
//
//      const double* yj = &y_flat[static_cast<std::size_t>(j) * nt];
//
//      for (int t0 = 0; t0 < nt; t0 += TB) {
//       const int tmax = std::min(t0 + TB, nt);
//
// #pragma omp simd
//       for (int t = t0; t < tmax; ++t) {
//        xty_i[t] += xj * yj[t];
//       }
//      }
//     }
//    }
//
//    xx_local[static_cast<std::size_t>(ii)] += x2_sum;
//   }
//  }
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//  for (int ii = 0; ii < mlen; ++ii) {
//   xx_block[static_cast<std::size_t>(ii)] += xx_local[static_cast<std::size_t>(ii)];
//
//   double* dst = &xty_block[static_cast<std::size_t>(ii) * nt];
//   const double* src = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//   for (int t = 0; t < nt; ++t) {
//    dst[t] += src[t];
//   }
//  }
// }
// }
//
// // -------------------------------------------------------------------------
// // Copy block accumulators to output
// // -------------------------------------------------------------------------
// for (int ii = 0; ii < mlen; ++ii) {
//  const int i = i0 + ii;
//
//  xx[i] = xx_block[static_cast<std::size_t>(ii)];
//
//  double* dst = &xty_flat[static_cast<std::size_t>(i) * nt];
//  const double* src = &xty_block[static_cast<std::size_t>(ii) * nt];
//
//  for (int t = 0; t < nt; ++t) {
//   dst[t] = src[t];
//  }
// }
//   }
//  } catch (...) {
//   std::free(block_buffer);
//   std::fclose(file_stream);
//   throw;
//  }
//
//  std::free(block_buffer);
//  std::fclose(file_stream);
// }
//
// // =============================================================================
// // Wrapper-neutral C++ convenience function
// // =============================================================================
//
// BedXtXytyResult bed_xtx_xty_cpp(
//   const std::string& bed_file,
//   int n,
//   const std::vector<int>& cls,
//   const std::vector<double>& af,
//   bool use_af,
//   const std::vector<double>& y_flat,
//   int nt,
//   bool scale,
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  if (bed_file.empty()) {
//   throw std::runtime_error("bed_file must not be empty.");
//  }
//
//  if (n <= 0) {
//   throw std::runtime_error("n must be positive.");
//  }
//
//  const int m = static_cast<int>(cls.size());
//
//  if (m <= 0) {
//   throw std::runtime_error("cls must be non-empty.");
//  }
//
//  if (nt <= 0) {
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (use_af && static_cast<int>(af.size()) != m) {
//   throw std::runtime_error("af must have same length as cls when use_af is true.");
//  }
//
//  const std::size_t expected_y_size =
//   static_cast<std::size_t>(n) * static_cast<std::size_t>(nt);
//
//  if (y_flat.size() != expected_y_size) {
//   throw std::runtime_error("y_flat must have length n * nt.");
//  }
//
//  if (MG <= 0) {
//   throw std::runtime_error("MG must be positive.");
//  }
//
//  if (JB <= 0) {
//   throw std::runtime_error("JB must be positive.");
//  }
//
//  if (TB <= 0) {
//   throw std::runtime_error("TB must be positive.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   if (cls[static_cast<std::size_t>(i)] <= 0) {
//    throw std::runtime_error("cls must contain positive 1-based marker indices.");
//   }
//
//   if (use_af) {
//    const double p = af[static_cast<std::size_t>(i)];
//
//    if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//     throw std::runtime_error("af must contain finite values in [0, 1].");
//    }
//   }
//  }
//
//  for (std::size_t k = 0; k < y_flat.size(); ++k) {
//   if (!std::isfinite(y_flat[k])) {
//    throw std::runtime_error("y_flat contains NaN/Inf.");
//   }
//  }
//
//  std::vector<double> xty_flat(static_cast<std::size_t>(m) * nt, 0.0);
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//  std::vector<double> af_used(static_cast<std::size_t>(m), 0.0);
//
//  bed_xtx_diag_xty_core(
//   bed_file.c_str(),
//   n,
//   cls.data(),
//   use_af ? af.data() : nullptr,
//   use_af,
//   af_used.data(),
//   m,
//   nt,
//   scale,
//   y_flat.data(),
//   xty_flat.data(),
//   xx.data(),
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  BedXtXytyResult out;
//  out.n = n;
//  out.m = m;
//  out.nt = nt;
//  out.scale = scale;
//  out.xx = xx;
//  out.af = af_used;
//  out.af_computed = !use_af;
//
//  out.wy.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
//  out.ww.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
//  out.yy.assign(static_cast<std::size_t>(nt), 0.0);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    out.wy[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
//     xty_flat[static_cast<std::size_t>(i) * nt + t];
//
//    out.ww[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
//     xx[static_cast<std::size_t>(i)];
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   double s = 0.0;
//
//   for (int j = 0; j < n; ++j) {
//    const double v = y_flat[static_cast<std::size_t>(j) * nt + t];
//    s += v * v;
//   }
//
//   out.yy[static_cast<std::size_t>(t)] = s;
//  }
//
//  out.min_xx = *std::min_element(xx.begin(), xx.end());
//  out.max_xx = *std::max_element(xx.begin(), xx.end());
//  out.mean_xx = std::accumulate(xx.begin(), xx.end(), 0.0) /
//   static_cast<double>(xx.size());
//
//  return out;
// }
//
// // =============================================================================
// // Thin Rcpp adapter
// // =============================================================================
//
// // [[Rcpp::export]]
// Rcpp::List bed_xtx_xty(
//   std::string bed_file,
//   int n,
//   Rcpp::IntegerVector cls,
//   Rcpp::Nullable<Rcpp::NumericVector> af,
//   Rcpp::NumericMatrix y,
//   bool scale = true,
//   int nthreads = 1,
//   int MG = 64,
//   int JB = 1024,
//   int TB = 32
// ) {
//  const int m = cls.size();
//  const int nt = y.ncol();
//
//  if (y.nrow() != n) {
//   Rcpp::stop("y must have n rows.");
//  }
//
//  std::vector<int> cls_cpp(static_cast<std::size_t>(m));
//  std::vector<double> af_cpp;
//  std::vector<double> y_flat(static_cast<std::size_t>(n) * nt);
//
//  for (int i = 0; i < m; ++i) {
//   if (cls[i] == NA_INTEGER) {
//    Rcpp::stop("cls contains NA.");
//   }
//
//   cls_cpp[static_cast<std::size_t>(i)] = cls[i];
//  }
//
//  const bool use_af = af.isNotNull();
//
//  if (use_af) {
//   Rcpp::NumericVector af_r = Rcpp::as<Rcpp::NumericVector>(af.get());
//
//   if (af_r.size() != m) {
//    Rcpp::stop("af must have same length as cls when supplied.");
//   }
//
//   af_cpp.resize(static_cast<std::size_t>(m));
//
//   for (int i = 0; i < m; ++i) {
//    af_cpp[static_cast<std::size_t>(i)] = af_r[i];
//   }
//  }
//
//  for (int j = 0; j < n; ++j) {
//   for (int t = 0; t < nt; ++t) {
//    y_flat[static_cast<std::size_t>(j) * nt + t] = y(j, t);
//   }
//  }
//
//  BedXtXytyResult res = bed_xtx_xty_cpp(
//   bed_file,
//   n,
//   cls_cpp,
//   af_cpp,
//   use_af,
//   y_flat,
//   nt,
//   scale,
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("wy") = res.wy,
//   Rcpp::Named("ww") = res.ww,
//   Rcpp::Named("yy") = res.yy,
//   Rcpp::Named("xx") = res.xx,
//   Rcpp::Named("af") = res.af,
//   Rcpp::Named("af_computed") = res.af_computed,
//   Rcpp::Named("n") = res.n,
//   Rcpp::Named("m") = res.m,
//   Rcpp::Named("nt") = res.nt,
//   Rcpp::Named("scale") = res.scale,
//   Rcpp::Named("mean_xx") = res.mean_xx,
//   Rcpp::Named("min_xx") = res.min_xx,
//   Rcpp::Named("max_xx") = res.max_xx
//  );
// }
//

// // =============================================================================
// // mt_sufficient_statistics.cpp
// // =============================================================================
//
// #include "mt_sufficient_statistics_core.h"
//
// #include <Rcpp.h>
//
// #include <algorithm>
// #include <cmath>
// #include <cstdio>
// #include <cstdlib>
// #include <cstring>
// #include <numeric>
// #include <stdexcept>
// #include <string>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// // =============================================================================
// // BED sufficient-statistics core
// // =============================================================================
// //
// // BED coding assumed:
// //   00 -> dosage 2
// //   01 -> missing
// //   10 -> dosage 1
// //   11 -> dosage 0
// //
// // If scale = true:
// //   x = (dosage - 2p) / sqrt(2p(1-p))
// //   missing -> 0 after standardization, i.e. mean-imputed standardized value.
// //
// // If scale = false:
// //   x = raw dosage, missing -> 2p.
// //   This is uncentered raw dosage coding.
// //
// // y_flat is n x nt, row-major: y_flat[j * nt + t].
// // xty_flat is m x nt, row-major: xty_flat[i * nt + t].
// // cls is 1-based marker index into the BED file.
// //
// // =============================================================================
//
//
// inline double compute_af_from_packed_marker(
//   const unsigned char* buf,
//   int n
// ) {
//  double dosage_sum = 0.0;
//  int n_obs = 0;
//  const int nbytes = (n + 3) / 4;
//
//  for (int kb = 0; kb < nbytes; ++kb) {
//   unsigned char x = buf[kb];
//   const int jbase = kb << 2;
//
//   for (int pos = 0; pos < 4; ++pos) {
//    const int j = jbase + pos;
//
//    if (j >= n) {
//     break;
//    }
//
//    const unsigned int code = x & 3u;
//    x >>= 2;
//
//    if (code == 1u) {
//     continue; // missing
//    }
//
//    double dosage = 0.0;
//
//    if (code == 0u) {
//     dosage = 2.0;
//    } else if (code == 2u) {
//     dosage = 1.0;
//    } else {
//     dosage = 0.0;
//    }
//
//    dosage_sum += dosage;
//    ++n_obs;
//   }
//  }
//
//  if (n_obs <= 0) {
//   return 0.0;
//  }
//
//  return dosage_sum / (2.0 * static_cast<double>(n_obs));
// }
//
// // =============================================================================
// // BED sufficient-statistics core
// // =============================================================================
// //
// // If scale = true:
// //   x = (dosage - 2p) / sqrt(2p(1-p))
// //   missing -> 0 after standardization, i.e. mean-imputed standardized value.
// //
// // If scale = false:
// //   x = raw dosage, missing -> 2p.
// //   This is uncentered raw dosage coding.
// //
// // y_flat is n x nt, row-major: y_flat[j * nt + t].
// // xty_flat is m x nt, row-major: xty_flat[i * nt + t].
// // cls is 1-based marker index into the BED file.
// //
// // If use_af = false, allele frequencies are computed from the selected marker
// // bytes on the fly and stored in af_used.
// //
// // =============================================================================
//
// void bed_xtx_diag_xty_core(
//   const char* file,
//   int n,
//   const int* cls,
//   const double* af,
//   bool use_af,
//   double* af_used,
//   int m,
//   int nt,
//   bool scale,
//   const double* y_flat,
//   double* xty_flat,
//   double* xx,
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  FILE* file_stream = std::fopen(file, "rb");
//  if (!file_stream) {
//   throw std::runtime_error("Could not open BED file.");
//  }
//
//  if (n <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("n must be positive.");
//  }
//
//  if (m <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("m must be positive.");
//  }
//
//  if (nt <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (!af_used) {
//   std::fclose(file_stream);
//   throw std::runtime_error("af_used must not be null.");
//  }
//
//  if (use_af && !af) {
//   std::fclose(file_stream);
//   throw std::runtime_error("af must not be null when use_af is true.");
//  }
//
//  if (MG <= 0) MG = 64;
//  if (JB <= 0) JB = 1024;
//  if (TB <= 0) TB = 32;
//
//  // Validate PLINK BED header: magic bytes + SNP-major mode.
//  unsigned char magic[3];
//  if (std::fread(magic, sizeof(unsigned char), 3, file_stream) != 3) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Could not read BED header.");
//  }
//
//  if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
//   std::fclose(file_stream);
//   throw std::runtime_error(
//     "BED file must be PLINK SNP-major format: expected magic bytes 0x6c 0x1b 0x01."
//   );
//  }
//
//  const std::size_t nbytes = static_cast<std::size_t>((n + 3) / 4);
//
//  std::memset(xty_flat, 0, static_cast<std::size_t>(m) * nt * sizeof(double));
//  std::memset(xx, 0, static_cast<std::size_t>(m) * sizeof(double));
//  std::memset(af_used, 0, static_cast<std::size_t>(m) * sizeof(double));
//
//  unsigned char* block_buffer =
//   static_cast<unsigned char*>(std::malloc(static_cast<std::size_t>(MG) * nbytes));
//
//  if (!block_buffer) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Failed to allocate BED block buffer.");
//  }
//
//  std::vector<double> map0(static_cast<std::size_t>(MG));
//  std::vector<double> map1(static_cast<std::size_t>(MG));
//  std::vector<double> map2(static_cast<std::size_t>(MG));
//  std::vector<double> map3(static_cast<std::size_t>(MG));
//
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  try {
//   for (int i0 = 0; i0 < m; i0 += MG) {
//    const int imax = std::min(i0 + MG, m);
//    const int mlen = imax - i0;
//
//    // -------------------------------------------------------------------------
//    // Read marker block, compute/use allele frequencies, and build lookup tables
//    // -------------------------------------------------------------------------
//    for (int ii = 0; ii < mlen; ++ii) {
//     const int i = i0 + ii;
//
//     const long long offset =
//      3LL + static_cast<long long>(cls[i] - 1) * static_cast<long long>(nbytes);
//
//     if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
//      throw std::runtime_error("fseek failed while reading BED marker " + std::to_string(cls[i]) + ".");
//     }
//
//     unsigned char* buf_i = block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//     const std::size_t nbytes_read =
//      std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);
//
//     if (nbytes_read != nbytes) {
//      throw std::runtime_error(
//        "Error reading BED data: nbytes_read != nbytes for marker " +
//         std::to_string(cls[i]) + "."
//      );
//     }
//
//     double p = 0.0;
//
//     if (use_af) {
//      p = af[i];
//     } else {
//      p = compute_af_from_packed_marker(buf_i, n);
//     }
//
//     if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//      throw std::runtime_error(
//        "Invalid allele frequency for selected marker index " + std::to_string(i) +
//         ", BED marker " + std::to_string(cls[i]) + "."
//      );
//     }
//
//     af_used[i] = p;
//
//     if (scale) {
//      const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//      if (denom <= 0.0 || !std::isfinite(denom)) {
//       map0[static_cast<std::size_t>(ii)] = 0.0;
//       map1[static_cast<std::size_t>(ii)] = 0.0;
//       map2[static_cast<std::size_t>(ii)] = 0.0;
//       map3[static_cast<std::size_t>(ii)] = 0.0;
//      } else {
//       map0[static_cast<std::size_t>(ii)] = (2.0 - 2.0 * p) / denom;  // BED 00 -> dosage 2
//       map1[static_cast<std::size_t>(ii)] = 0.0;                      // BED 01 -> missing
//       map2[static_cast<std::size_t>(ii)] = (1.0 - 2.0 * p) / denom;  // BED 10 -> dosage 1
//       map3[static_cast<std::size_t>(ii)] = (0.0 - 2.0 * p) / denom;  // BED 11 -> dosage 0
//      }
//     } else {
//      // Raw dosage coding with missing mean-imputed to 2p.
//      // This is intentionally uncentered.
//      map0[static_cast<std::size_t>(ii)] = 2.0;
//      map1[static_cast<std::size_t>(ii)] = 2.0 * p;
//      map2[static_cast<std::size_t>(ii)] = 1.0;
//      map3[static_cast<std::size_t>(ii)] = 0.0;
//     }
//    }
//
//    // Block-level accumulators.
//    std::vector<double> xty_block(static_cast<std::size_t>(mlen) * nt, 0.0);
//    std::vector<double> xx_block(static_cast<std::size_t>(mlen), 0.0);
//
//    // -------------------------------------------------------------------------
//    // Parallel accumulation over individual blocks
//    // -------------------------------------------------------------------------
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// {
//  std::vector<double> xty_local(static_cast<std::size_t>(mlen) * nt, 0.0);
//  std::vector<double> xx_local(static_cast<std::size_t>(mlen), 0.0);
//
// #ifdef _OPENMP
// #pragma omp for schedule(static)
// #endif
//  for (int j0 = 0; j0 < n; j0 += JB) {
//   const int jmax = std::min(j0 + JB, n);
//
//   const int byte0 = j0 >> 2;
//   const int byte1 = (jmax + 3) >> 2;
//
//   for (int ii = 0; ii < mlen; ++ii) {
//    const unsigned char* buf_i =
//     block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//    const double m0 = map0[static_cast<std::size_t>(ii)];
//    const double m1 = map1[static_cast<std::size_t>(ii)];
//    const double m2 = map2[static_cast<std::size_t>(ii)];
//    const double m3 = map3[static_cast<std::size_t>(ii)];
//
//    double x2_sum = 0.0;
//    double* xty_i = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//    for (int kb = byte0; kb < byte1; ++kb) {
//     unsigned char buf_k = buf_i[kb];
//     const int jbase = kb << 2;
//
//     for (int pos = 0; pos < 4; ++pos) {
//      const int j = jbase + pos;
//
//      if (j < j0 || j >= jmax || j >= n) {
//       buf_k >>= 2;
//       continue;
//      }
//
//      double xj;
//      switch (buf_k & 3u) {
//      case 0u: xj = m0; break;
//      case 1u: xj = m1; break;
//      case 2u: xj = m2; break;
//      default: xj = m3; break;
//      }
//
//      buf_k >>= 2;
//
//      x2_sum += xj * xj;
//
//      const double* yj = &y_flat[static_cast<std::size_t>(j) * nt];
//
//      for (int t0 = 0; t0 < nt; t0 += TB) {
//       const int tmax = std::min(t0 + TB, nt);
//
// #pragma omp simd
//       for (int t = t0; t < tmax; ++t) {
//        xty_i[t] += xj * yj[t];
//       }
//      }
//     }
//    }
//
//    xx_local[static_cast<std::size_t>(ii)] += x2_sum;
//   }
//  }
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//  for (int ii = 0; ii < mlen; ++ii) {
//   xx_block[static_cast<std::size_t>(ii)] += xx_local[static_cast<std::size_t>(ii)];
//
//   double* dst = &xty_block[static_cast<std::size_t>(ii) * nt];
//   const double* src = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//   for (int t = 0; t < nt; ++t) {
//    dst[t] += src[t];
//   }
//  }
// }
// }
//
// // -------------------------------------------------------------------------
// // Copy block accumulators to output
// // -------------------------------------------------------------------------
// for (int ii = 0; ii < mlen; ++ii) {
//  const int i = i0 + ii;
//
//  xx[i] = xx_block[static_cast<std::size_t>(ii)];
//
//  double* dst = &xty_flat[static_cast<std::size_t>(i) * nt];
//  const double* src = &xty_block[static_cast<std::size_t>(ii) * nt];
//
//  for (int t = 0; t < nt; ++t) {
//   dst[t] = src[t];
//  }
// }
//   }
//  } catch (...) {
//   std::free(block_buffer);
//   std::fclose(file_stream);
//   throw;
//  }
//
//  std::free(block_buffer);
//  std::fclose(file_stream);
// }
//
// // =============================================================================
// // Wrapper-neutral C++ convenience function
// // =============================================================================
//
// BedXtXytyResult bed_xtx_xty_cpp(
//   const std::string& bed_file,
//   int n,
//   const std::vector<int>& cls,
//   const std::vector<double>& af,
//   bool use_af,
//   const std::vector<double>& y_flat,
//   int nt,
//   bool scale,
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  if (bed_file.empty()) {
//   throw std::runtime_error("bed_file must not be empty.");
//  }
//
//  if (n <= 0) {
//   throw std::runtime_error("n must be positive.");
//  }
//
//  const int m = static_cast<int>(cls.size());
//
//  if (m <= 0) {
//   throw std::runtime_error("cls must be non-empty.");
//  }
//
//  if (nt <= 0) {
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (use_af && static_cast<int>(af.size()) != m) {
//   throw std::runtime_error("af must have same length as cls when use_af is true.");
//  }
//
//  if (static_cast<int>(y_flat.size()) != n * nt) {
//   throw std::runtime_error("y_flat must have length n * nt.");
//  }
//
//  if (MG <= 0) {
//   throw std::runtime_error("MG must be positive.");
//  }
//
//  if (JB <= 0) {
//   throw std::runtime_error("JB must be positive.");
//  }
//
//  if (TB <= 0) {
//   throw std::runtime_error("TB must be positive.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   if (cls[static_cast<std::size_t>(i)] <= 0) {
//    throw std::runtime_error("cls must contain positive 1-based marker indices.");
//   }
//
//   if (use_af) {
//    const double p = af[static_cast<std::size_t>(i)];
//
//    if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//     throw std::runtime_error("af must contain finite values in [0, 1].");
//    }
//   }
//  }
//
//  for (std::size_t k = 0; k < y_flat.size(); ++k) {
//   if (!std::isfinite(y_flat[k])) {
//    throw std::runtime_error("y_flat contains NaN/Inf.");
//   }
//  }
//
//  std::vector<double> xty_flat(static_cast<std::size_t>(m) * nt, 0.0);
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//  std::vector<double> af_used(static_cast<std::size_t>(m), 0.0);
//
//  bed_xtx_diag_xty_core(
//   bed_file.c_str(),
//   n,
//   cls.data(),
//   use_af ? af.data() : nullptr,
//   use_af,
//   af_used.data(),
//   m,
//   nt,
//   scale,
//   y_flat.data(),
//   xty_flat.data(),
//   xx.data(),
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  BedXtXytyResult out;
//  out.n = n;
//  out.m = m;
//  out.nt = nt;
//  out.scale = scale;
//  out.xx = xx;
//  out.af = af_used;
//  out.af_computed = !use_af;
//
//  out.wy.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
//  out.ww.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
//  out.yy.assign(static_cast<std::size_t>(nt), 0.0);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    out.wy[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
//     xty_flat[static_cast<std::size_t>(i) * nt + t];
//
//    out.ww[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
//     xx[static_cast<std::size_t>(i)];
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   double s = 0.0;
//
//   for (int j = 0; j < n; ++j) {
//    const double v = y_flat[static_cast<std::size_t>(j) * nt + t];
//    s += v * v;
//   }
//
//   out.yy[static_cast<std::size_t>(t)] = s;
//  }
//
//  out.min_xx = *std::min_element(xx.begin(), xx.end());
//  out.max_xx = *std::max_element(xx.begin(), xx.end());
//  out.mean_xx = std::accumulate(xx.begin(), xx.end(), 0.0) /
//   static_cast<double>(xx.size());
//
//  return out;
// }
//
// // =============================================================================
// // Thin Rcpp adapter
// // =============================================================================
//
// // [[Rcpp::export]]
// Rcpp::List bed_xtx_xty(
//   std::string bed_file,
//   int n,
//   Rcpp::IntegerVector cls,
//   Rcpp::Nullable<Rcpp::NumericVector> af,
//   Rcpp::NumericMatrix y,
//   bool scale = true,
//   int nthreads = 1,
//   int MG = 64,
//   int JB = 1024,
//   int TB = 32
// ) {
//  const int m = cls.size();
//  const int nt = y.ncol();
//
//  if (y.nrow() != n) {
//   Rcpp::stop("y must have n rows.");
//  }
//
//  std::vector<int> cls_cpp(static_cast<std::size_t>(m));
//  std::vector<double> af_cpp;
//  std::vector<double> y_flat(static_cast<std::size_t>(n) * nt);
//
//  for (int i = 0; i < m; ++i) {
//   if (cls[i] == NA_INTEGER) {
//    Rcpp::stop("cls contains NA.");
//   }
//
//   cls_cpp[static_cast<std::size_t>(i)] = cls[i];
//  }
//
//  const bool use_af = af.isNotNull();
//
//  if (use_af) {
//   Rcpp::NumericVector af_r = Rcpp::as<Rcpp::NumericVector>(af.get());
//
//   if (af_r.size() != m) {
//    Rcpp::stop("af must have same length as cls when supplied.");
//   }
//
//   af_cpp.resize(static_cast<std::size_t>(m));
//
//   for (int i = 0; i < m; ++i) {
//    af_cpp[static_cast<std::size_t>(i)] = af_r[i];
//   }
//  }
//
//  for (int j = 0; j < n; ++j) {
//   for (int t = 0; t < nt; ++t) {
//    y_flat[static_cast<std::size_t>(j) * nt + t] = y(j, t);
//   }
//  }
//
//  BedXtXytyResult res = bed_xtx_xty_cpp(
//   bed_file,
//   n,
//   cls_cpp,
//   af_cpp,
//   use_af,
//   y_flat,
//   nt,
//   scale,
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("wy") = res.wy,
//   Rcpp::Named("ww") = res.ww,
//   Rcpp::Named("yy") = res.yy,
//   Rcpp::Named("xx") = res.xx,
//   Rcpp::Named("af") = res.af,
//   Rcpp::Named("af_computed") = res.af_computed,
//   Rcpp::Named("n") = res.n,
//   Rcpp::Named("m") = res.m,
//   Rcpp::Named("nt") = res.nt,
//   Rcpp::Named("scale") = res.scale,
//   Rcpp::Named("mean_xx") = res.mean_xx,
//   Rcpp::Named("min_xx") = res.min_xx,
//   Rcpp::Named("max_xx") = res.max_xx
//  );
// }




// void bed_xtx_diag_xty_core(
//   const char* file,
//   int n,
//   const int* cls,
//   const double* af,
//   int m,
//   int nt,
//   bool scale,
//   const double* y_flat,
//   double* xty_flat,
//   double* xx,
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  FILE* file_stream = std::fopen(file, "rb");
//  if (!file_stream) {
//   throw std::runtime_error("Could not open BED file.");
//  }
//
//  if (n <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("n must be positive.");
//  }
//
//  if (m <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("m must be positive.");
//  }
//
//  if (nt <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (MG <= 0) MG = 64;
//  if (JB <= 0) JB = 1024;
//  if (TB <= 0) TB = 32;
//
//  // Validate PLINK BED header: magic bytes + SNP-major mode.
//  unsigned char magic[3];
//  if (std::fread(magic, sizeof(unsigned char), 3, file_stream) != 3) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Could not read BED header.");
//  }
//
//  if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
//   std::fclose(file_stream);
//   throw std::runtime_error(
//     "BED file must be PLINK SNP-major format: expected magic bytes 0x6c 0x1b 0x01."
//   );
//  }
//
//  const std::size_t nbytes = static_cast<std::size_t>((n + 3) / 4);
//
//  std::memset(xty_flat, 0, static_cast<std::size_t>(m) * nt * sizeof(double));
//  std::memset(xx, 0, static_cast<std::size_t>(m) * sizeof(double));
//
//  unsigned char* block_buffer =
//   static_cast<unsigned char*>(std::malloc(static_cast<std::size_t>(MG) * nbytes));
//
//  if (!block_buffer) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Failed to allocate BED block buffer.");
//  }
//
//  std::vector<double> map0(static_cast<std::size_t>(MG));
//  std::vector<double> map1(static_cast<std::size_t>(MG));
//  std::vector<double> map2(static_cast<std::size_t>(MG));
//  std::vector<double> map3(static_cast<std::size_t>(MG));
//
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  try {
//   for (int i0 = 0; i0 < m; i0 += MG) {
//    const int imax = std::min(i0 + MG, m);
//    const int mlen = imax - i0;
//
//    // -------------------------------------------------------------------------
//    // Read marker block and build genotype lookup tables
//    // -------------------------------------------------------------------------
//    for (int ii = 0; ii < mlen; ++ii) {
//     const int i = i0 + ii;
//
//     const long long offset =
//      3LL + static_cast<long long>(cls[i] - 1) * static_cast<long long>(nbytes);
//
//     if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
//      throw std::runtime_error("fseek failed while reading BED marker " + std::to_string(cls[i]) + ".");
//     }
//
//     unsigned char* buf_i = block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//     const std::size_t nbytes_read =
//      std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);
//
//     if (nbytes_read != nbytes) {
//      throw std::runtime_error(
//        "Error reading BED data: nbytes_read != nbytes for marker " +
//         std::to_string(cls[i]) + "."
//      );
//     }
//
//     const double p = af[i];
//
//     if (scale) {
//      const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//      if (denom <= 0.0 || !std::isfinite(denom)) {
//       map0[static_cast<std::size_t>(ii)] = 0.0;
//       map1[static_cast<std::size_t>(ii)] = 0.0;
//       map2[static_cast<std::size_t>(ii)] = 0.0;
//       map3[static_cast<std::size_t>(ii)] = 0.0;
//      } else {
//       map0[static_cast<std::size_t>(ii)] = (2.0 - 2.0 * p) / denom;  // BED 00 -> dosage 2
//       map1[static_cast<std::size_t>(ii)] = 0.0;                      // BED 01 -> missing
//       map2[static_cast<std::size_t>(ii)] = (1.0 - 2.0 * p) / denom;  // BED 10 -> dosage 1
//       map3[static_cast<std::size_t>(ii)] = (0.0 - 2.0 * p) / denom;  // BED 11 -> dosage 0
//      }
//     } else {
//      // Raw dosage coding with missing mean-imputed to 2p.
//      // This is intentionally uncentered.
//      map0[static_cast<std::size_t>(ii)] = 2.0;
//      map1[static_cast<std::size_t>(ii)] = 2.0 * p;
//      map2[static_cast<std::size_t>(ii)] = 1.0;
//      map3[static_cast<std::size_t>(ii)] = 0.0;
//     }
//    }
//
//    // Block-level accumulators.
//    std::vector<double> xty_block(static_cast<std::size_t>(mlen) * nt, 0.0);
//    std::vector<double> xx_block(static_cast<std::size_t>(mlen), 0.0);
//
//    // -------------------------------------------------------------------------
//    // Parallel accumulation over individual blocks
//    // -------------------------------------------------------------------------
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// {
//  std::vector<double> xty_local(static_cast<std::size_t>(mlen) * nt, 0.0);
//  std::vector<double> xx_local(static_cast<std::size_t>(mlen), 0.0);
//
// #ifdef _OPENMP
// #pragma omp for schedule(static)
// #endif
//  for (int j0 = 0; j0 < n; j0 += JB) {
//   const int jmax = std::min(j0 + JB, n);
//
//   const int byte0 = j0 >> 2;
//   const int byte1 = (jmax + 3) >> 2;
//
//   for (int ii = 0; ii < mlen; ++ii) {
//    const unsigned char* buf_i =
//     block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//    const double m0 = map0[static_cast<std::size_t>(ii)];
//    const double m1 = map1[static_cast<std::size_t>(ii)];
//    const double m2 = map2[static_cast<std::size_t>(ii)];
//    const double m3 = map3[static_cast<std::size_t>(ii)];
//
//    double x2_sum = 0.0;
//    double* xty_i = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//    for (int kb = byte0; kb < byte1; ++kb) {
//     unsigned char buf_k = buf_i[kb];
//     const int jbase = kb << 2;
//
//     for (int pos = 0; pos < 4; ++pos) {
//      const int j = jbase + pos;
//
//      if (j < j0 || j >= jmax || j >= n) {
//       buf_k >>= 2;
//       continue;
//      }
//
//      double xj;
//      switch (buf_k & 3u) {
//      case 0u: xj = m0; break;
//      case 1u: xj = m1; break;
//      case 2u: xj = m2; break;
//      default: xj = m3; break;
//      }
//
//      buf_k >>= 2;
//
//      x2_sum += xj * xj;
//
//      const double* yj = &y_flat[static_cast<std::size_t>(j) * nt];
//
//      for (int t0 = 0; t0 < nt; t0 += TB) {
//       const int tmax = std::min(t0 + TB, nt);
//
// #pragma omp simd
//       for (int t = t0; t < tmax; ++t) {
//        xty_i[t] += xj * yj[t];
//       }
//      }
//     }
//    }
//
//    xx_local[static_cast<std::size_t>(ii)] += x2_sum;
//   }
//  }
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//  for (int ii = 0; ii < mlen; ++ii) {
//   xx_block[static_cast<std::size_t>(ii)] += xx_local[static_cast<std::size_t>(ii)];
//
//   double* dst = &xty_block[static_cast<std::size_t>(ii) * nt];
//   const double* src = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//   for (int t = 0; t < nt; ++t) {
//    dst[t] += src[t];
//   }
//  }
// }
// }
//
// // -------------------------------------------------------------------------
// // Copy block accumulators to output
// // -------------------------------------------------------------------------
// for (int ii = 0; ii < mlen; ++ii) {
//  const int i = i0 + ii;
//
//  xx[i] = xx_block[static_cast<std::size_t>(ii)];
//
//  double* dst = &xty_flat[static_cast<std::size_t>(i) * nt];
//  const double* src = &xty_block[static_cast<std::size_t>(ii) * nt];
//
//  for (int t = 0; t < nt; ++t) {
//   dst[t] = src[t];
//  }
// }
//   }
//  } catch (...) {
//   std::free(block_buffer);
//   std::fclose(file_stream);
//   throw;
//  }
//
//  std::free(block_buffer);
//  std::fclose(file_stream);
// }
//
// // =============================================================================
// // Wrapper-neutral C++ convenience function
// // =============================================================================
//
// BedXtXytyResult bed_xtx_xty_cpp(
//   const std::string& bed_file,
//   int n,
//   const std::vector<int>& cls,
//   const std::vector<double>& af,
//   const std::vector<double>& y_flat,
//   int nt,
//   bool scale,
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  if (bed_file.empty()) {
//   throw std::runtime_error("bed_file must not be empty.");
//  }
//
//  if (n <= 0) {
//   throw std::runtime_error("n must be positive.");
//  }
//
//  const int m = static_cast<int>(cls.size());
//
//  if (m <= 0) {
//   throw std::runtime_error("cls must be non-empty.");
//  }
//
//  if (nt <= 0) {
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (static_cast<int>(af.size()) != m) {
//   throw std::runtime_error("af must have same length as cls.");
//  }
//
//  if (static_cast<int>(y_flat.size()) != n * nt) {
//   throw std::runtime_error("y_flat must have length n * nt.");
//  }
//
//  if (MG <= 0) {
//   throw std::runtime_error("MG must be positive.");
//  }
//
//  if (JB <= 0) {
//   throw std::runtime_error("JB must be positive.");
//  }
//
//  if (TB <= 0) {
//   throw std::runtime_error("TB must be positive.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   if (cls[static_cast<std::size_t>(i)] <= 0) {
//    throw std::runtime_error("cls must contain positive 1-based marker indices.");
//   }
//
//   const double p = af[static_cast<std::size_t>(i)];
//
//   if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//    throw std::runtime_error("af must contain finite values in [0, 1].");
//   }
//  }
//
//  for (std::size_t k = 0; k < y_flat.size(); ++k) {
//   if (!std::isfinite(y_flat[k])) {
//    throw std::runtime_error("y_flat contains NaN/Inf.");
//   }
//  }
//
//  std::vector<double> xty_flat(static_cast<std::size_t>(m) * nt, 0.0);
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//
//  bed_xtx_diag_xty_core(
//   bed_file.c_str(),
//   n,
//   cls.data(),
//   af.data(),
//   m,
//   nt,
//   scale,
//   y_flat.data(),
//   xty_flat.data(),
//   xx.data(),
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  BedXtXytyResult out;
//  out.n = n;
//  out.m = m;
//  out.nt = nt;
//  out.scale = scale;
//  out.xx = xx;
//
//  out.wy.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
//  out.ww.assign(static_cast<std::size_t>(nt), std::vector<double>(static_cast<std::size_t>(m), 0.0));
//  out.yy.assign(static_cast<std::size_t>(nt), 0.0);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    out.wy[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
//     xty_flat[static_cast<std::size_t>(i) * nt + t];
//
//    out.ww[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)] =
//     xx[static_cast<std::size_t>(i)];
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   double s = 0.0;
//
//   for (int j = 0; j < n; ++j) {
//    const double v = y_flat[static_cast<std::size_t>(j) * nt + t];
//    s += v * v;
//   }
//
//   out.yy[static_cast<std::size_t>(t)] = s;
//  }
//
//  out.min_xx = *std::min_element(xx.begin(), xx.end());
//  out.max_xx = *std::max_element(xx.begin(), xx.end());
//  out.mean_xx = std::accumulate(xx.begin(), xx.end(), 0.0) /
//   static_cast<double>(xx.size());
//
//  return out;
// }
//
// // =============================================================================
// // Thin Rcpp adapter
// // =============================================================================
//
// // [[Rcpp::export]]
// Rcpp::List bed_xtx_xty(
//   std::string bed_file,
//   int n,
//   Rcpp::IntegerVector cls,
//   Rcpp::NumericVector af,
//   Rcpp::NumericMatrix y,
//   bool scale = true,
//   int nthreads = 1,
//   int MG = 64,
//   int JB = 1024,
//   int TB = 32
// ) {
//  const int m = cls.size();
//  const int nt = y.ncol();
//
//  if (y.nrow() != n) {
//   Rcpp::stop("y must have n rows.");
//  }
//
//  if (af.size() != m) {
//   Rcpp::stop("af must have same length as cls.");
//  }
//
//  std::vector<int> cls_cpp(static_cast<std::size_t>(m));
//  std::vector<double> af_cpp(static_cast<std::size_t>(m));
//  std::vector<double> y_flat(static_cast<std::size_t>(n) * nt);
//
//  for (int i = 0; i < m; ++i) {
//   if (cls[i] == NA_INTEGER) {
//    Rcpp::stop("cls contains NA.");
//   }
//
//   cls_cpp[static_cast<std::size_t>(i)] = cls[i];
//   af_cpp[static_cast<std::size_t>(i)] = af[i];
//  }
//
//  for (int j = 0; j < n; ++j) {
//   for (int t = 0; t < nt; ++t) {
//    y_flat[static_cast<std::size_t>(j) * nt + t] = y(j, t);
//   }
//  }
//
//  BedXtXytyResult res = bed_xtx_xty_cpp(
//   bed_file,
//   n,
//   cls_cpp,
//   af_cpp,
//   y_flat,
//   nt,
//   scale,
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("wy") = res.wy,
//   Rcpp::Named("ww") = res.ww,
//   Rcpp::Named("yy") = res.yy,
//   Rcpp::Named("xx") = res.xx,
//   Rcpp::Named("n") = res.n,
//   Rcpp::Named("m") = res.m,
//   Rcpp::Named("nt") = res.nt,
//   Rcpp::Named("scale") = res.scale,
//   Rcpp::Named("mean_xx") = res.mean_xx,
//   Rcpp::Named("min_xx") = res.min_xx,
//   Rcpp::Named("max_xx") = res.max_xx
//  );
// }
//

// #include "mt_sufficient_statistics_core.h"
//
// #include <Rcpp.h>
// #include <vector>
// #include <string>
//
// #include <cmath>
// #include <cstdio>
// #include <vector>
// #include <stdexcept>
// #include <algorithm>
// #include <cstdlib>
// #include <string>
// #include <cstring>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// void bed_xtx_diag_xty_core(
//   const char* file,
//   int n,
//   const int* cls,
//   const double* af,
//   int m,
//   int nt,
//   bool scale,
//   const double* y_flat,     // n x nt, row-major: y_flat[j * nt + t]
//   double* xty_flat,         // m x nt, row-major: xty_flat[i * nt + t]
//   double* xx,               // length m
//   int nthreads,
//   int MG,
//   int JB,
//   int TB
// ) {
//  FILE* file_stream = std::fopen(file, "rb");
//  if (!file_stream) {
//   throw std::runtime_error("Could not open BED file.");
//  }
//
//  if (n <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("n must be positive.");
//  }
//
//  if (m <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("m must be positive.");
//  }
//
//  if (nt <= 0) {
//   std::fclose(file_stream);
//   throw std::runtime_error("nt must be positive.");
//  }
//
//  if (MG <= 0) MG = 64;
//  if (JB <= 0) JB = 1024;
//  if (TB <= 0) TB = 32;
//
//  const std::size_t nbytes = static_cast<std::size_t>((n + 3) / 4);
//
//  std::memset(xty_flat, 0, static_cast<std::size_t>(m) * nt * sizeof(double));
//  std::memset(xx, 0, static_cast<std::size_t>(m) * sizeof(double));
//
//  unsigned char* block_buffer =
//   static_cast<unsigned char*>(std::malloc(static_cast<std::size_t>(MG) * nbytes));
//
//  if (!block_buffer) {
//   std::fclose(file_stream);
//   throw std::runtime_error("Failed to allocate BED block buffer.");
//  }
//
//  std::vector<double> map0(MG), map1(MG), map2(MG), map3(MG);
//
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  try {
//   for (int i0 = 0; i0 < m; i0 += MG) {
//    const int imax = std::min(i0 + MG, m);
//    const int mlen = imax - i0;
//
//    // ------------------------------------------------------------
//    // Read marker block and build genotype lookup tables
//    // ------------------------------------------------------------
//    for (int ii = 0; ii < mlen; ++ii) {
//     const int i = i0 + ii;
//
//     const long long offset =
//      3LL + static_cast<long long>(cls[i] - 1) * static_cast<long long>(nbytes);
//
//     if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
//      throw std::runtime_error("fseek failed.");
//     }
//
//     unsigned char* buf_i = block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//     const std::size_t nbytes_read =
//      std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);
//
//     if (nbytes_read != nbytes) {
//      throw std::runtime_error("Error reading BED data: nbytes_read != nbytes.");
//     }
//
//     const double p = af[i];
//
//     if (scale) {
//      const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//      if (denom <= 0.0 || !std::isfinite(denom)) {
//       map0[ii] = 0.0;
//       map1[ii] = 0.0;
//       map2[ii] = 0.0;
//       map3[ii] = 0.0;
//      } else {
//       map0[ii] = (2.0 - 2.0 * p) / denom;  // BED 00 -> dosage 2
//       map1[ii] = 0.0;                      // BED 01 -> missing
//       map2[ii] = (1.0 - 2.0 * p) / denom;  // BED 10 -> dosage 1
//       map3[ii] = (0.0 - 2.0 * p) / denom;  // BED 11 -> dosage 0
//      }
//     } else {
//      // Raw dosage coding with missing mean-imputed to 2p.
//      // If you instead want centered unscaled coding, use:
//      //   map0 = 2 - 2p, map1 = 0, map2 = 1 - 2p, map3 = -2p
//      map0[ii] = 2.0;
//      map1[ii] = 2.0 * p;   // missing -> mean dosage
//      map2[ii] = 1.0;
//      map3[ii] = 0.0;
//     }
//    }
//
//    // Block-level accumulators
//    std::vector<double> xty_block(static_cast<std::size_t>(mlen) * nt, 0.0);
//    std::vector<double> xx_block(static_cast<std::size_t>(mlen), 0.0);
//
//    // ------------------------------------------------------------
//    // Parallel accumulation over individual blocks
//    // ------------------------------------------------------------
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// {
//  std::vector<double> xty_local(static_cast<std::size_t>(mlen) * nt, 0.0);
//  std::vector<double> xx_local(static_cast<std::size_t>(mlen), 0.0);
//
// #ifdef _OPENMP
// #pragma omp for schedule(static)
// #endif
//  for (int j0 = 0; j0 < n; j0 += JB) {
//   const int jmax = std::min(j0 + JB, n);
//
//   const int byte0 = j0 >> 2;
//   const int byte1 = (jmax + 3) >> 2;
//
//   for (int ii = 0; ii < mlen; ++ii) {
//    const unsigned char* buf_i =
//     block_buffer + static_cast<std::size_t>(ii) * nbytes;
//
//    const double m0 = map0[ii];
//    const double m1 = map1[ii];
//    const double m2 = map2[ii];
//    const double m3 = map3[ii];
//
//    double x2_sum = 0.0;
//    double* xty_i = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//    for (int kb = byte0; kb < byte1; ++kb) {
//     unsigned char buf_k = buf_i[kb];
//     const int jbase = kb << 2;
//
//     for (int pos = 0; pos < 4; ++pos) {
//      const int j = jbase + pos;
//
//      if (j < j0 || j >= jmax || j >= n) {
//       buf_k >>= 2;
//       continue;
//      }
//
//      double xj;
//      switch (buf_k & 3u) {
//      case 0u: xj = m0; break;
//      case 1u: xj = m1; break;
//      case 2u: xj = m2; break;
//      default: xj = m3; break;
//      }
//
//      buf_k >>= 2;
//
//      x2_sum += xj * xj;
//
//      const double* yj = &y_flat[static_cast<std::size_t>(j) * nt];
//
//      for (int t0 = 0; t0 < nt; t0 += TB) {
//       const int tmax = std::min(t0 + TB, nt);
//
// #pragma omp simd
//       for (int t = t0; t < tmax; ++t) {
//        xty_i[t] += xj * yj[t];
//       }
//      }
//     }
//    }
//
//    xx_local[ii] += x2_sum;
//   }
//  }
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//  for (int ii = 0; ii < mlen; ++ii) {
//   xx_block[ii] += xx_local[ii];
//
//   double* dst = &xty_block[static_cast<std::size_t>(ii) * nt];
//   const double* src = &xty_local[static_cast<std::size_t>(ii) * nt];
//
//   for (int t = 0; t < nt; ++t) {
//    dst[t] += src[t];
//   }
//  }
// }
// }
//
// // ------------------------------------------------------------
// // Copy block accumulators to output
// // ------------------------------------------------------------
// for (int ii = 0; ii < mlen; ++ii) {
//  const int i = i0 + ii;
//
//  xx[i] = xx_block[ii];
//
//  double* dst = &xty_flat[static_cast<std::size_t>(i) * nt];
//  const double* src = &xty_block[static_cast<std::size_t>(ii) * nt];
//
//  for (int t = 0; t < nt; ++t) {
//   dst[t] = src[t];
//  }
// }
//   }
//  } catch (...) {
//   std::free(block_buffer);
//   std::fclose(file_stream);
//   throw;
//  }
//
//  std::free(block_buffer);
//  std::fclose(file_stream);
// }
//
//
// // [[Rcpp::export]]
// Rcpp::List bed_xtx_xty(
//   std::string bed_file,
//   int n,
//   Rcpp::IntegerVector cls,
//   Rcpp::NumericVector af,
//   Rcpp::NumericMatrix y,
//   bool scale = true,
//   int nthreads = 1,
//   int MG = 64,
//   int JB = 1024,
//   int TB = 32
// ) {
//  const int m = cls.size();
//  const int nt = y.ncol();
//
//  if (y.nrow() != n) {
//   Rcpp::stop("y must have n rows.");
//  }
//  if (af.size() != m) {
//   Rcpp::stop("af must have same length as cls.");
//  }
//
//  std::vector<int> cls_cpp(m);
//  std::vector<double> af_cpp(m);
//
//  for (int i = 0; i < m; ++i) {
//   cls_cpp[i] = cls[i];
//   af_cpp[i] = af[i];
//  }
//
//  // y_flat is n x nt, row-major: y_flat[j * nt + t]
//  std::vector<double> y_flat(static_cast<std::size_t>(n) * nt);
//
//  for (int j = 0; j < n; ++j) {
//   for (int t = 0; t < nt; ++t) {
//    y_flat[static_cast<std::size_t>(j) * nt + t] = y(j, t);
//   }
//  }
//
//  // Core output:
//  // xty_flat is m x nt, row-major: xty_flat[i * nt + t]
//  std::vector<double> xty_flat(static_cast<std::size_t>(m) * nt, 0.0);
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//
//  bed_xtx_diag_xty_core(
//   bed_file.c_str(),
//   n,
//   cls_cpp.data(),
//   af_cpp.data(),
//   m,
//   nt,
//   scale,
//   y_flat.data(),
//   xty_flat.data(),
//   xx.data(),
//   nthreads,
//   MG,
//   JB,
//   TB
//  );
//
//  // Return as std::vector<std::vector<double>> with trait-major layout:
//  // wy[t][i], ww[t][i]
//  std::vector<std::vector<double>> wy(
//    nt,
//    std::vector<double>(m, 0.0)
//  );
//
//  std::vector<std::vector<double>> ww(
//    nt,
//    std::vector<double>(m, 0.0)
//  );
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    wy[t][i] = xty_flat[static_cast<std::size_t>(i) * nt + t];
//    ww[t][i] = xx[i];
//   }
//  }
//
//  // yy[t] = y_t' y_t
//  std::vector<double> yy(nt, 0.0);
//
//  for (int t = 0; t < nt; ++t) {
//   double s = 0.0;
//   for (int j = 0; j < n; ++j) {
//    const double v = y(j, t);
//    s += v * v;
//   }
//   yy[t] = s;
//  }
//
//  return Rcpp::List::create(
//   Rcpp::Named("wy") = wy,
//   Rcpp::Named("ww") = ww,
//   Rcpp::Named("yy") = yy,
//   Rcpp::Named("n") = n,
//   Rcpp::Named("m") = m,
//   Rcpp::Named("nt") = nt,
//   Rcpp::Named("scale") = scale
//  );
// }
