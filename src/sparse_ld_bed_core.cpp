// sparse_ld_bed_core.cpp
//
// Compute sparse LD from one or more PLINK BED-style 2-bit genotype files
// and write CSR to disk.
//
// BED coding assumed:
//   00 -> dosage 2
//   01 -> missing
//   10 -> dosage 1
//   11 -> dosage 0
//
// Missing genotypes are mean-imputed after standardization, i.e. z = 0.
// LD values stored are true correlations r = x_i'x_j / sqrt(x_i'x_i x_j'x_j), not r^2.
// Streaming CSR output stores upper triangle only, with implicit diagonal 1.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Rcpp.h>
#include <R_ext/BLAS.h>

#include "packed_bed.h"
#include "st_bed_decode.h"

// extern "C" {
// #include <cblas.h>
// }

// [[Rcpp::export]]
Rcpp::List sparseLD_thread_info(int nthreads = 0) {
#ifdef _OPENMP
 int actual_threads_default = 1;
 int actual_threads_requested = 1;

#pragma omp parallel
{
#pragma omp single
{
 actual_threads_default = omp_get_num_threads();
}
}

if (nthreads > 0) {
#pragma omp parallel num_threads(nthreads)
{
#pragma omp single
{
 actual_threads_requested = omp_get_num_threads();
}
}
} else {
 actual_threads_requested = actual_threads_default;
}

return Rcpp::List::create(
 Rcpp::Named("openmp") = true,
 Rcpp::Named("requested_nthreads") = nthreads,
 Rcpp::Named("actual_threads_default_region") = actual_threads_default,
 Rcpp::Named("actual_threads_requested_region") = actual_threads_requested,
 Rcpp::Named("omp_get_max_threads") = omp_get_max_threads(),
 Rcpp::Named("omp_get_num_procs") = omp_get_num_procs()
);
#else
return Rcpp::List::create(
 Rcpp::Named("openmp") = false,
 Rcpp::Named("requested_nthreads") = nthreads,
 Rcpp::Named("actual_threads_default_region") = 1,
 Rcpp::Named("actual_threads_requested_region") = 1,
 Rcpp::Named("omp_get_max_threads") = 1,
 Rcpp::Named("omp_get_num_procs") = 1
);
#endif
}

// -----------------------------------------------------------------------------
// Small helpers
// -----------------------------------------------------------------------------

static inline void checked_fwrite(
  const void* data,
  std::size_t elem_size,
  std::size_t n_elem,
  FILE* fs,
  const std::string& what
) {
 if (n_elem == 0) return;
 const std::size_t got = std::fwrite(data, elem_size, n_elem, fs);
 if (got != n_elem) throw std::runtime_error("Short write while writing " + what);
}

static inline void validate_sparse_ld_args(
  int n,
  int m,
  int block_size,
  float r2_threshold
) {
 if (n <= 1) {
  throw std::runtime_error("Need at least two samples to compute LD.");
 }

 if (m <= 0) {
  throw std::runtime_error("Need at least one marker to compute LD.");
 }

 if (block_size <= 0) {
  throw std::runtime_error("block_size must be positive.");
 }

 if (!std::isfinite(r2_threshold) || r2_threshold < 0.0f) {
  throw std::runtime_error("r2_threshold must be finite and non-negative.");
 }
}

static inline void validate_af_and_pos(
  const std::vector<double>& af_cpp,
  const std::vector<int>& pos_cpp,
  int m
) {
 if (static_cast<int>(af_cpp.size()) != m) {
  throw std::runtime_error("af must have one value per selected marker.");
 }

 for (int i = 0; i < m; ++i) {
  const double p = af_cpp[static_cast<std::size_t>(i)];

  if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
   throw std::runtime_error("af contains invalid values; expected finite values in [0, 1].");
  }
 }

 if (!pos_cpp.empty() && static_cast<int>(pos_cpp.size()) != m) {
  throw std::runtime_error("pos_bp must have one value per selected marker.");
 }
}

static inline double xtx_to_corr(
  float xij,
  const std::vector<double>& xx,
  int i,
  int j
) {
 const double denom = std::sqrt(
  xx[static_cast<std::size_t>(i)] *
   xx[static_cast<std::size_t>(j)]
 );

 if (denom <= 0.0 || !std::isfinite(denom)) {
  return 0.0;
 }

 const double r = static_cast<double>(xij) / denom;

 if (!std::isfinite(r)) {
  return 0.0;
 }

 return r;
}

// -----------------------------------------------------------------------------
// CSR objects and LD helpers
// -----------------------------------------------------------------------------

struct SparseLDCSR {
 int m = 0;
 std::vector<uint64_t> row_ptr;
 std::vector<uint32_t> col_idx;
 std::vector<float> values;
};

struct StreamingCSRResult {
 int m = 0;
 int n_used = 0;
 int n_bed = 0;
 uint64_t nnz = 0;
 std::string row_file;
 std::string col_file;
 std::string val_file;
 std::string meta_file;
};

// static void compute_ld_block_sgemm(
//   const float* ZA,
//   const float* ZB,
//   int n,
//   int ma,
//   int mb,
//   float* R
// ) {
//  // Raw crossproduct: R_ij = x_i' x_j.
//  // Correlations are computed explicitly as R_ij / sqrt(xx_i * xx_j)
//  // before thresholding and storage.
//  const float alpha = 1.0f;
//  const float beta = 0.0f;
// 
//  cblas_sgemm(
//   CblasRowMajor,
//   CblasNoTrans,
//   CblasTrans,
//   ma,
//   mb,
//   n,
//   alpha,
//   ZA,
//   n,
//   ZB,
//   n,
//   beta,
//   R,
//   mb
//  );
// }
static void compute_ld_block_sgemm(
    const float* ZA,
    const float* ZB,
    int n,
    int ma,
    int mb,
    float* R
) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int ia = 0; ia < ma; ++ia) {
    const float* za = ZA + static_cast<std::size_t>(ia) * n;
    
    for (int ib = 0; ib < mb; ++ib) {
      const float* zb = ZB + static_cast<std::size_t>(ib) * n;
      
      double acc = 0.0;
      for (int k = 0; k < n; ++k) {
        acc += static_cast<double>(za[k]) * static_cast<double>(zb[k]);
      }
      
      R[static_cast<std::size_t>(ia) * mb + ib] = static_cast<float>(acc);
    }
  }
}

static inline bool within_window(
  int i,
  int j,
  const int* pos_bp,
  int max_distance_bp,
  int max_distance_variants
) {
 if (j <= i) return false;
 if (max_distance_variants > 0 && j - i > max_distance_variants) return false;
 if (pos_bp && max_distance_bp > 0 && pos_bp[j] - pos_bp[i] > max_distance_bp) return false;
 return true;
}

static inline bool block_may_overlap_window(
  int a0,
  int ma,
  int b0,
  const int* pos_bp,
  int max_distance_bp,
  int max_distance_variants
) {
 const int last_a = a0 + ma - 1;
 if (max_distance_variants > 0 && b0 - last_a > max_distance_variants) return false;
 if (pos_bp && max_distance_bp > 0 && pos_bp[b0] - pos_bp[last_a] > max_distance_bp) return false;
 return true;
}

// -----------------------------------------------------------------------------
// In-memory CSR core
// -----------------------------------------------------------------------------

static SparseLDCSR sparse_ld_packed_core(
  const PackedBedMatrix& G,
  const double* af,
  const int* pos_bp,
  int max_distance_bp,
  int max_distance_variants,
  float r2_threshold,
  int block_size,
  int nthreads
) {
#ifdef _OPENMP
 if (nthreads > 0) {
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
 }
#endif

 const int n = G.n;
 const int m = G.m;
 validate_sparse_ld_args(n, m, block_size, r2_threshold);

 const int nth = effective_nthreads(nthreads);
 const std::vector<double> xx = compute_xx_from_packed_standardized(G, af, nth);
 const int nb = (m + block_size - 1) / block_size;
 std::vector<uint64_t> row_counts(static_cast<std::size_t>(m), 0);

 auto count_or_fill = [&](SparseLDCSR* csr_ptr) {
  std::vector<uint64_t> write_pos;
  if (csr_ptr) write_pos = csr_ptr->row_ptr;

  std::vector<float> ZA(static_cast<std::size_t>(block_size) * n);
  std::vector<float> ZB(static_cast<std::size_t>(block_size) * n);
  std::vector<float> R(static_cast<std::size_t>(block_size) * block_size);

  for (int a_block = 0; a_block < nb; ++a_block) {
   const int a0 = a_block * block_size;
   const int ma = std::min(block_size, m - a0);
   decode_packed_block_float(G, a0, ma, af, ZA.data(), nth);

   for (int b_block = a_block; b_block < nb; ++b_block) {
    const int b0 = b_block * block_size;
    const int mb = std::min(block_size, m - b0);

    if (!block_may_overlap_window(a0, ma, b0, pos_bp, max_distance_bp, max_distance_variants)) break;

    if (b_block == a_block) {
     compute_ld_block_sgemm(ZA.data(), ZA.data(), n, ma, ma, R.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
     for (int ia = 0; ia < ma; ++ia) {
      const int i = a0 + ia;
      uint64_t cnt = 0;
      uint64_t pos = csr_ptr ? write_pos[static_cast<std::size_t>(i)] : 0;

      for (int ib = ia + 1; ib < ma; ++ib) {
       const int j = b0 + ib;
       if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;

       const float xij = R[static_cast<std::size_t>(ia) * ma + ib];
       const double r = xtx_to_corr(xij, xx, i, j);
       const double r2 = r * r;

       if (r2 >= static_cast<double>(r2_threshold)) {
        if (csr_ptr) {
         csr_ptr->col_idx[static_cast<std::size_t>(pos)] = static_cast<uint32_t>(j);
         csr_ptr->values[static_cast<std::size_t>(pos)] = static_cast<float>(r);
         ++pos;
        } else {
         ++cnt;
        }
       }
      }

      if (csr_ptr) write_pos[static_cast<std::size_t>(i)] = pos;
      else row_counts[static_cast<std::size_t>(i)] += cnt;
     }
    } else {
     decode_packed_block_float(G, b0, mb, af, ZB.data(), nth);
     compute_ld_block_sgemm(ZA.data(), ZB.data(), n, ma, mb, R.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
     for (int ia = 0; ia < ma; ++ia) {
      const int i = a0 + ia;
      uint64_t cnt = 0;
      uint64_t pos = csr_ptr ? write_pos[static_cast<std::size_t>(i)] : 0;

      for (int ib = 0; ib < mb; ++ib) {
       const int j = b0 + ib;
       if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;

       const float xij = R[static_cast<std::size_t>(ia) * mb + ib];
       const double r = xtx_to_corr(xij, xx, i, j);
       const double r2 = r * r;

       if (r2 >= static_cast<double>(r2_threshold)) {
        if (csr_ptr) {
         csr_ptr->col_idx[static_cast<std::size_t>(pos)] = static_cast<uint32_t>(j);
         csr_ptr->values[static_cast<std::size_t>(pos)] = static_cast<float>(r);
         ++pos;
        } else {
         ++cnt;
        }
       }
      }

      if (csr_ptr) write_pos[static_cast<std::size_t>(i)] = pos;
      else row_counts[static_cast<std::size_t>(i)] += cnt;
     }
    }
   }
  }

#ifndef NDEBUG
  if (csr_ptr) {
   for (int i = 0; i < m; ++i) {
    assert(write_pos[static_cast<std::size_t>(i)] == csr_ptr->row_ptr[static_cast<std::size_t>(i) + 1]);
   }
  }
#endif
 };

 count_or_fill(nullptr);

 SparseLDCSR csr;
 csr.m = m;
 csr.row_ptr.resize(static_cast<std::size_t>(m) + 1);
 csr.row_ptr[0] = 0;
 for (int i = 0; i < m; ++i) {
  csr.row_ptr[static_cast<std::size_t>(i) + 1] =
   csr.row_ptr[static_cast<std::size_t>(i)] + row_counts[static_cast<std::size_t>(i)];
 }

 const uint64_t nnz = csr.row_ptr[static_cast<std::size_t>(m)];
 csr.col_idx.resize(static_cast<std::size_t>(nnz));
 csr.values.resize(static_cast<std::size_t>(nnz));

 count_or_fill(&csr);
 return csr;
}

// -----------------------------------------------------------------------------
// Streaming CSR core
// -----------------------------------------------------------------------------

static StreamingCSRResult sparse_ld_packed_stream_to_disk_core(
  const PackedBedMatrix& G,
  int n_bed,
  const double* af,
  const int* pos_bp,
  int max_distance_bp,
  int max_distance_variants,
  float r2_threshold,
  int block_size,
  int nthreads,
  const std::string& out_prefix
) {
#ifdef _OPENMP
 if (nthreads > 0) {
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
 }
#endif

 const int n = G.n;
 const int m = G.m;
 validate_sparse_ld_args(n, m, block_size, r2_threshold);

 const int nth = effective_nthreads(nthreads);
 const std::vector<double> xx = compute_xx_from_packed_standardized(G, af, nth);
 const int nb = (m + block_size - 1) / block_size;

 StreamingCSRResult res;
 res.m = m;
 res.n_used = n;
 res.n_bed = n_bed;
 res.nnz = 0;
 res.row_file = out_prefix + ".row_ptr.u64.bin";
 res.col_file = out_prefix + ".col_idx.u32.0based.bin";
 res.val_file = out_prefix + ".values.f32.bin";
 res.meta_file = out_prefix + ".meta.txt";

 FILE* row_fs = std::fopen(res.row_file.c_str(), "wb");
 FILE* col_fs = std::fopen(res.col_file.c_str(), "wb");
 FILE* val_fs = std::fopen(res.val_file.c_str(), "wb");

 if (!row_fs || !col_fs || !val_fs) {
  if (row_fs) std::fclose(row_fs);
  if (col_fs) std::fclose(col_fs);
  if (val_fs) std::fclose(val_fs);
  throw std::runtime_error("Could not open one or more CSR output files.");
 }

 try {
  std::vector<float> ZA(static_cast<std::size_t>(block_size) * n);
  std::vector<float> ZB(static_cast<std::size_t>(block_size) * n);
  std::vector<float> R(static_cast<std::size_t>(block_size) * block_size);

  uint64_t cumulative_nnz = 0;
  checked_fwrite(&cumulative_nnz, sizeof(uint64_t), 1, row_fs, "row_ptr[0]");

  for (int a_block = 0; a_block < nb; ++a_block) {
   const int a0 = a_block * block_size;
   const int ma = std::min(block_size, m - a0);

   decode_packed_block_float(G, a0, ma, af, ZA.data(), nth);

   std::vector<std::vector<uint32_t>> block_cols(static_cast<std::size_t>(ma));
   std::vector<std::vector<float>> block_vals(static_cast<std::size_t>(ma));

   for (int ia = 0; ia < ma; ++ia) {
    block_cols[static_cast<std::size_t>(ia)].reserve(128);
    block_vals[static_cast<std::size_t>(ia)].reserve(128);
   }

   for (int b_block = a_block; b_block < nb; ++b_block) {
    const int b0 = b_block * block_size;
    const int mb = std::min(block_size, m - b0);

    if (!block_may_overlap_window(a0, ma, b0, pos_bp, max_distance_bp, max_distance_variants)) break;

    if (b_block == a_block) {
     compute_ld_block_sgemm(ZA.data(), ZA.data(), n, ma, ma, R.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
     for (int ia = 0; ia < ma; ++ia) {
      const int i = a0 + ia;
      std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
      std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];

      for (int ib = ia + 1; ib < ma; ++ib) {
       const int j = b0 + ib;
       if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;

       const float xij = R[static_cast<std::size_t>(ia) * ma + ib];
       const double r = xtx_to_corr(xij, xx, i, j);
       const double r2 = r * r;

       if (r2 >= static_cast<double>(r2_threshold)) {
        cols.push_back(static_cast<uint32_t>(j));
        vals.push_back(static_cast<float>(r));
       }
      }
     }
    } else {
     decode_packed_block_float(G, b0, mb, af, ZB.data(), nth);
     compute_ld_block_sgemm(ZA.data(), ZB.data(), n, ma, mb, R.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
     for (int ia = 0; ia < ma; ++ia) {
      const int i = a0 + ia;
      std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
      std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];

      for (int ib = 0; ib < mb; ++ib) {
       const int j = b0 + ib;
       if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;

       const float xij = R[static_cast<std::size_t>(ia) * mb + ib];
       const double r = xtx_to_corr(xij, xx, i, j);
       const double r2 = r * r;

       if (r2 >= static_cast<double>(r2_threshold)) {
        cols.push_back(static_cast<uint32_t>(j));
        vals.push_back(static_cast<float>(r));
       }
      }
     }
    }
   }

   for (int ia = 0; ia < ma; ++ia) {
    const std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
    const std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];

    checked_fwrite(cols.data(), sizeof(uint32_t), cols.size(), col_fs, "col_idx");
    checked_fwrite(vals.data(), sizeof(float), vals.size(), val_fs, "values");

    cumulative_nnz += static_cast<uint64_t>(cols.size());
    checked_fwrite(&cumulative_nnz, sizeof(uint64_t), 1, row_fs, "row_ptr");
   }
  }

  res.nnz = cumulative_nnz;

  std::fclose(row_fs);
  std::fclose(col_fs);
  std::fclose(val_fs);
  row_fs = nullptr;
  col_fs = nullptr;
  val_fs = nullptr;

  std::ofstream meta(res.meta_file.c_str());
  if (!meta.is_open()) throw std::runtime_error("Could not open metadata output file.");

  meta << "format=sparse_ld_csr" << std::endl;
  meta << "storage=streamed_upper_triangle" << std::endl;
  meta << "n_bed=" << n_bed << std::endl;
  meta << "n_used=" << n << std::endl;
  meta << "n_samples_used=" << n << std::endl;
  meta << "n_variants=" << m << std::endl;
  meta << "nnz=" << static_cast<unsigned long long>(res.nnz) << std::endl;
  meta << "triangle=upper" << std::endl;
  meta << "diagonal=implicit_1" << std::endl;
  meta << "row_ptr_file=" << res.row_file << std::endl;
  meta << "col_idx_file=" << res.col_file << std::endl;
  meta << "values_file=" << res.val_file << std::endl;
  meta << "row_ptr_type=uint64" << std::endl;
  meta << "col_idx_type=uint32" << std::endl;
  meta << "values_type=float32" << std::endl;
  meta << "index_base=0" << std::endl;
  meta << "value=r" << std::endl;
  meta << "ld_normalization=sqrt_xx" << std::endl;
  meta << "sgemm_value_before_normalization=xij" << std::endl;
  meta << "r2_threshold=" << static_cast<double>(r2_threshold) << std::endl;
  meta << "max_distance_bp=" << max_distance_bp << std::endl;
  meta << "max_distance_variants=" << max_distance_variants << std::endl;
  meta << "block_size=" << block_size << std::endl;
  meta << "nthreads=" << nthreads << std::endl;
  meta.close();

 } catch (...) {
  if (row_fs) std::fclose(row_fs);
  if (col_fs) std::fclose(col_fs);
  if (val_fs) std::fclose(val_fs);
  throw;
 }

 return res;
}

// -----------------------------------------------------------------------------
// Rcpp helpers
// -----------------------------------------------------------------------------

static Rcpp::NumericVector u64_to_numeric(const std::vector<uint64_t>& x) {
 Rcpp::NumericVector out(x.size());
 for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<double>(x[i]);
 return out;
}

static Rcpp::IntegerVector u32_to_integer_1based(const std::vector<uint32_t>& x) {
 Rcpp::IntegerVector out(x.size());
 for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<int>(x[i]) + 1;
 return out;
}

static Rcpp::NumericVector f32_to_numeric(const std::vector<float>& x) {
 Rcpp::NumericVector out(x.size());
 for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<double>(x[i]);
 return out;
}

static std::vector<std::string> copy_bed_files(Rcpp::CharacterVector bed_files) {
 std::vector<std::string> out(static_cast<std::size_t>(bed_files.size()));
 for (int i = 0; i < bed_files.size(); ++i) out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(bed_files[i]);
 return out;
}

static std::vector<std::vector<int>> copy_int_list(Rcpp::List xlist) {
 std::vector<std::vector<int>> out(static_cast<std::size_t>(xlist.size()));
 for (int f = 0; f < xlist.size(); ++f) {
  Rcpp::IntegerVector x = xlist[f];
  out[static_cast<std::size_t>(f)].resize(static_cast<std::size_t>(x.size()));
  for (int i = 0; i < x.size(); ++i) out[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)] = x[i];
 }
 return out;
}

static std::vector<double> flatten_af_list(Rcpp::Nullable<Rcpp::List> af) {
 std::vector<double> out;
 if (af.isNotNull()) {
  Rcpp::List af_list = Rcpp::as<Rcpp::List>(af.get());
  for (int f = 0; f < af_list.size(); ++f) {
   Rcpp::NumericVector x = af_list[f];
   for (int i = 0; i < x.size(); ++i) out.push_back(x[i]);
  }
 }
 return out;
}

static std::vector<int> flatten_pos_list_or_empty(Rcpp::Nullable<Rcpp::List> pos_bp) {
 std::vector<int> out;
 if (pos_bp.isNotNull()) {
  Rcpp::List pos_list = Rcpp::as<Rcpp::List>(pos_bp.get());
  for (int f = 0; f < pos_list.size(); ++f) {
   Rcpp::IntegerVector x = pos_list[f];
   for (int i = 0; i < x.size(); ++i) out.push_back(x[i]);
  }
 }
 return out;
}

static std::vector<int> copy_rows0_or_empty(Rcpp::Nullable<Rcpp::IntegerVector> rows, int n) {
 std::vector<int> out;
 if (rows.isNotNull()) {
  Rcpp::IntegerVector r(rows);
  out.resize(static_cast<std::size_t>(r.size()));

  for (int i = 0; i < r.size(); ++i) {
   if (r[i] == NA_INTEGER) throw std::runtime_error("rows contains NA.");
   if (r[i] <= 0) throw std::runtime_error("rows must contain positive 1-based indices.");
   out[static_cast<std::size_t>(i)] = r[i] - 1;
  }

  if (static_cast<int>(out.size()) == n) {
   bool identity = true;
   for (int i = 0; i < n; ++i) {
    if (out[static_cast<std::size_t>(i)] != i) {
     identity = false;
     break;
    }
   }
   if (identity) out.clear();
  }
 }
 return out;
}

// -----------------------------------------------------------------------------
// Rcpp exports
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List sparseLD_stream_CSR(
  Rcpp::CharacterVector bed_files,
  int n,
  Rcpp::List cls,
  std::string out_prefix,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::Nullable<Rcpp::List> af = R_NilValue,
  Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
  int max_distance_bp = 1000000,
  int max_distance_variants = 1000,
  double r2_threshold = 0.01,
  int block_size = 1024,
  int nthreads = 1
) {
 if (out_prefix.empty()) Rcpp::stop("out_prefix must not be empty.");

 std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
 std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
 std::vector<int> rows0 = copy_rows0_or_empty(rows, n);

 const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
 const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());

 PackedBedMatrix G = read_bedfiles_to_packed_matrix(
  bed_files_cpp,
  n,
  rows0_ptr,
  n_rows,
  cls_by_file
 );

 std::vector<double> af_cpp = flatten_af_list(af);
 const bool af_computed = af_cpp.empty();
 if (af_computed) af_cpp = compute_af_from_packed(G);

 std::vector<int> pos_cpp = flatten_pos_list_or_empty(pos_bp);
 validate_af_and_pos(af_cpp, pos_cpp, G.m);
 const int* pos_ptr = pos_cpp.empty() ? nullptr : pos_cpp.data();

 StreamingCSRResult res = sparse_ld_packed_stream_to_disk_core(
  G,
  n,
  af_cpp.data(),
  pos_ptr,
  max_distance_bp,
  max_distance_variants,
  static_cast<float>(r2_threshold),
  block_size,
  nthreads,
  out_prefix
 );

 return Rcpp::List::create(
  Rcpp::Named("row_ptr_file") = res.row_file,
  Rcpp::Named("col_idx_file") = res.col_file,
  Rcpp::Named("values_file") = res.val_file,
  Rcpp::Named("meta_file") = res.meta_file,
  Rcpp::Named("nrow") = res.m,
  Rcpp::Named("ncol") = res.m,
  Rcpp::Named("max_distance_bp") = max_distance_bp,
  Rcpp::Named("max_distance_variants") = max_distance_variants,
  Rcpp::Named("r2_threshold") = r2_threshold,
  Rcpp::Named("block_size") = block_size,
  Rcpp::Named("nthreads") = nthreads,
  Rcpp::Named("nnz") = static_cast<double>(res.nnz),
  Rcpp::Named("n_bed") = n,
  Rcpp::Named("n_used") = G.n,
  Rcpp::Named("af_computed") = af_computed,
  Rcpp::Named("upper_triangle") = true,
  Rcpp::Named("diag") = "implicit_1",
  Rcpp::Named("index_base") = 0,
  Rcpp::Named("value") = "r",
  Rcpp::Named("ld_normalization") = "sqrt_xx"
 );
}

// [[Rcpp::export]]
Rcpp::List sparseLD_to_CSR(
  Rcpp::CharacterVector bed_files,
  int n,
  Rcpp::List cls,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::Nullable<Rcpp::List> af = R_NilValue,
  Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
  int max_distance_bp = 1000000,
  int max_distance_variants = 1000,
  double r2_threshold = 0.01,
  int block_size = 1024,
  int nthreads = 1
) {
 std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
 std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
 std::vector<int> rows0 = copy_rows0_or_empty(rows, n);

 const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
 const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());

 PackedBedMatrix G = read_bedfiles_to_packed_matrix(
  bed_files_cpp,
  n,
  rows0_ptr,
  n_rows,
  cls_by_file
 );

 std::vector<double> af_cpp = flatten_af_list(af);
 const bool af_computed = af_cpp.empty();
 if (af_computed) af_cpp = compute_af_from_packed(G);

 std::vector<int> pos_cpp = flatten_pos_list_or_empty(pos_bp);
 validate_af_and_pos(af_cpp, pos_cpp, G.m);
 const int* pos_ptr = pos_cpp.empty() ? nullptr : pos_cpp.data();

 SparseLDCSR csr = sparse_ld_packed_core(
  G,
  af_cpp.data(),
  pos_ptr,
  max_distance_bp,
  max_distance_variants,
  static_cast<float>(r2_threshold),
  block_size,
  nthreads
 );

 return Rcpp::List::create(
  Rcpp::Named("row_ptr") = u64_to_numeric(csr.row_ptr),
  Rcpp::Named("col_idx") = u32_to_integer_1based(csr.col_idx),
  Rcpp::Named("values") = f32_to_numeric(csr.values),
  Rcpp::Named("nrow") = G.m,
  Rcpp::Named("ncol") = G.m,
  Rcpp::Named("nnz") = static_cast<double>(csr.row_ptr.back()),
  Rcpp::Named("n_bed") = n,
  Rcpp::Named("n_used") = G.n,
  Rcpp::Named("max_distance_bp") = max_distance_bp,
  Rcpp::Named("max_distance_variants") = max_distance_variants,
  Rcpp::Named("r2_threshold") = r2_threshold,
  Rcpp::Named("block_size") = block_size,
  Rcpp::Named("nthreads") = nthreads,
  Rcpp::Named("af_computed") = af_computed,
  Rcpp::Named("upper_triangle") = true,
  Rcpp::Named("diag") = "implicit_1",
  Rcpp::Named("index_base") = 1,
  Rcpp::Named("value") = "r",
  Rcpp::Named("ld_normalization") = "sqrt_xx"
 );
}

// [[Rcpp::export]]
Rcpp::List sparseLD_write_CSR(
  Rcpp::CharacterVector bed_files,
  int n,
  Rcpp::List cls,
  std::string out_prefix,
  Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
  Rcpp::Nullable<Rcpp::List> af = R_NilValue,
  Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
  int max_distance_bp = 1000000,
  int max_distance_variants = 1000,
  double r2_threshold = 0.01,
  int block_size = 1024,
  int nthreads = 1
) {
 return sparseLD_stream_CSR(
  bed_files,
  n,
  cls,
  out_prefix,
  rows,
  af,
  pos_bp,
  max_distance_bp,
  max_distance_variants,
  r2_threshold,
  block_size,
  nthreads
 );
}

// [[Rcpp::export]]
Rcpp::List sparseLD_read_CSR(
  std::string prefix,
  bool one_based = true
) {
 const std::string meta_file = prefix + ".meta.txt";
 const std::string row_file = prefix + ".row_ptr.u64.bin";
 const std::string col_file = prefix + ".col_idx.u32.0based.bin";
 const std::string val_file = prefix + ".values.f32.bin";

 std::ifstream meta(meta_file.c_str());
 if (!meta.is_open()) Rcpp::stop("Could not open metadata file: %s", meta_file);

 int m = -1;
 double nnz_d = -1.0;

 int max_distance_bp = -1;
 int max_distance_variants = -1;
 double r2_threshold = R_NaN;
 int block_size = -1;
 int nthreads = -1;
 int n_bed = -1;
 int n_used = -1;
 std::string value_type = "r";
 std::string ld_normalization = "";

 std::string line;
 while (std::getline(meta, line)) {
  const std::string key_m = "n_variants=";
  const std::string key_nnz = "nnz=";
  const std::string key_max_bp = "max_distance_bp=";
  const std::string key_max_var = "max_distance_variants=";
  const std::string key_r2 = "r2_threshold=";
  const std::string key_block = "block_size=";
  const std::string key_threads = "nthreads=";
  const std::string key_n_bed = "n_bed=";
  const std::string key_n_used = "n_used=";
  const std::string key_n_samples_used = "n_samples_used=";
  const std::string key_value = "value=";
  const std::string key_ld_norm = "ld_normalization=";

  if (line.rfind(key_m, 0) == 0) {
   m = std::stoi(line.substr(key_m.size()));
  } else if (line.rfind(key_nnz, 0) == 0) {
   nnz_d = std::stod(line.substr(key_nnz.size()));
  } else if (line.rfind(key_max_bp, 0) == 0) {
   max_distance_bp = std::stoi(line.substr(key_max_bp.size()));
  } else if (line.rfind(key_max_var, 0) == 0) {
   max_distance_variants = std::stoi(line.substr(key_max_var.size()));
  } else if (line.rfind(key_r2, 0) == 0) {
   r2_threshold = std::stod(line.substr(key_r2.size()));
  } else if (line.rfind(key_block, 0) == 0) {
   block_size = std::stoi(line.substr(key_block.size()));
  } else if (line.rfind(key_threads, 0) == 0) {
   nthreads = std::stoi(line.substr(key_threads.size()));
  } else if (line.rfind(key_n_bed, 0) == 0) {
   n_bed = std::stoi(line.substr(key_n_bed.size()));
  } else if (line.rfind(key_n_used, 0) == 0) {
   n_used = std::stoi(line.substr(key_n_used.size()));
  } else if (line.rfind(key_n_samples_used, 0) == 0 && n_used < 0) {
   n_used = std::stoi(line.substr(key_n_samples_used.size()));
  } else if (line.rfind(key_value, 0) == 0) {
   value_type = line.substr(key_value.size());
  } else if (line.rfind(key_ld_norm, 0) == 0) {
   ld_normalization = line.substr(key_ld_norm.size());
  }
 }
 meta.close();

 if (m <= 0) Rcpp::stop("Could not read n_variants from metadata.");
 if (nnz_d < 0.0) Rcpp::stop("Could not read nnz from metadata.");

 const std::size_t nnz = static_cast<std::size_t>(nnz_d);

 std::vector<uint64_t> row_ptr_u64(static_cast<std::size_t>(m) + 1);
 std::vector<uint32_t> col_idx_u32(nnz);
 std::vector<float> values_f32(nnz);

 auto read_file = [](const std::string& path, void* data, std::size_t nbytes) {
  FILE* fs = std::fopen(path.c_str(), "rb");
  if (!fs) throw std::runtime_error("Could not open file: " + path);
  const std::size_t got = std::fread(data, 1, nbytes, fs);
  std::fclose(fs);
  if (got != nbytes) throw std::runtime_error("Short read from file: " + path);
 };

 read_file(row_file, row_ptr_u64.data(), row_ptr_u64.size() * sizeof(uint64_t));
 read_file(col_file, col_idx_u32.data(), col_idx_u32.size() * sizeof(uint32_t));
 read_file(val_file, values_f32.data(), values_f32.size() * sizeof(float));

 if (row_ptr_u64[0] != 0 || row_ptr_u64[static_cast<std::size_t>(m)] != static_cast<uint64_t>(nnz)) {
  Rcpp::stop("Invalid row_ptr: expected 0-based pointer ending at nnz.");
 }

 for (int i = 0; i < m; ++i) {
  if (row_ptr_u64[static_cast<std::size_t>(i + 1)] < row_ptr_u64[static_cast<std::size_t>(i)]) {
   Rcpp::stop("Invalid row_ptr: row pointers are not nondecreasing.");
  }
 }

 Rcpp::NumericVector row_ptr(row_ptr_u64.size());
 for (std::size_t i = 0; i < row_ptr_u64.size(); ++i) {
  row_ptr[i] = static_cast<double>(row_ptr_u64[i]);
 }

 Rcpp::IntegerVector col_idx(col_idx_u32.size());
 for (std::size_t i = 0; i < col_idx_u32.size(); ++i) {
  if (col_idx_u32[i] >= static_cast<uint32_t>(m)) {
   Rcpp::stop("LD column index out of range.");
  }
  col_idx[i] = static_cast<int>(col_idx_u32[i]) + (one_based ? 1 : 0);
 }

 Rcpp::NumericVector values(values_f32.size());
 for (std::size_t i = 0; i < values_f32.size(); ++i) {
  if (!std::isfinite(values_f32[i])) {
   Rcpp::stop("LD value contains NaN/Inf.");
  }
  values[i] = static_cast<double>(values_f32[i]);
 }

 return Rcpp::List::create(
  Rcpp::Named("row_ptr") = row_ptr,
  Rcpp::Named("col_idx") = col_idx,
  Rcpp::Named("values") = values,
  Rcpp::Named("nrow") = m,
  Rcpp::Named("ncol") = m,
  Rcpp::Named("nnz") = static_cast<double>(nnz),
  Rcpp::Named("upper_triangle") = true,
  Rcpp::Named("diag") = "implicit_1",
  Rcpp::Named("index_base") = one_based ? 1 : 0,
  Rcpp::Named("value") = value_type,
  Rcpp::Named("ld_normalization") = ld_normalization,
  Rcpp::Named("max_distance_bp") =
   max_distance_bp >= 0 ? Rcpp::wrap(max_distance_bp) : Rcpp::wrap(NA_INTEGER),
    Rcpp::Named("max_distance_variants") =
     max_distance_variants >= 0 ? Rcpp::wrap(max_distance_variants) : Rcpp::wrap(NA_INTEGER),
      Rcpp::Named("r2_threshold") =
       std::isfinite(r2_threshold) ? Rcpp::wrap(r2_threshold) : Rcpp::wrap(NA_REAL),
       Rcpp::Named("block_size") =
        block_size >= 0 ? Rcpp::wrap(block_size) : Rcpp::wrap(NA_INTEGER),
         Rcpp::Named("nthreads") =
          nthreads >= 0 ? Rcpp::wrap(nthreads) : Rcpp::wrap(NA_INTEGER),
           Rcpp::Named("n_bed") =
            n_bed >= 0 ? Rcpp::wrap(n_bed) : Rcpp::wrap(NA_INTEGER),
             Rcpp::Named("n_used") =
              n_used >= 0 ? Rcpp::wrap(n_used) : Rcpp::wrap(NA_INTEGER)
 );
}


// // sparse_ld_bed_core.cpp
// //
// // Compute sparse LD from one or more PLINK BED-style 2-bit genotype files
// // and write CSR to disk.
// //
// // BED coding assumed:
// //   00 -> dosage 2
// //   01 -> missing
// //   10 -> dosage 1
// //   11 -> dosage 0
// //
// // Missing genotypes are mean-imputed after standardization, i.e. z = 0.
// // LD values stored are correlations r, not r^2.
// // Streaming CSR output stores upper triangle only, with implicit diagonal 1.
//
// #include <algorithm>
// #include <cassert>
// #include <cmath>
// #include <cstdint>
// #include <cstdio>
// #include <fstream>
// #include <numeric>
// #include <stdexcept>
// #include <string>
// #include <unordered_map>
// #include <unordered_set>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// #include <Rcpp.h>
//
// #include "packed_bed.h"
//
// extern "C" {
// #include <cblas.h>
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_thread_info(int nthreads = 0) {
// #ifdef _OPENMP
//  int actual_threads_default = 1;
//  int actual_threads_requested = 1;
//
// #pragma omp parallel
// {
// #pragma omp single
// {
//  actual_threads_default = omp_get_num_threads();
// }
// }
//
// if (nthreads > 0) {
// #pragma omp parallel num_threads(nthreads)
// {
// #pragma omp single
// {
//  actual_threads_requested = omp_get_num_threads();
// }
// }
// } else {
//  actual_threads_requested = actual_threads_default;
// }
//
// return Rcpp::List::create(
//  Rcpp::Named("openmp") = true,
//  Rcpp::Named("requested_nthreads") = nthreads,
//  Rcpp::Named("actual_threads_default_region") = actual_threads_default,
//  Rcpp::Named("actual_threads_requested_region") = actual_threads_requested,
//  Rcpp::Named("omp_get_max_threads") = omp_get_max_threads(),
//  Rcpp::Named("omp_get_num_procs") = omp_get_num_procs()
// );
// #else
// return Rcpp::List::create(
//  Rcpp::Named("openmp") = false,
//  Rcpp::Named("requested_nthreads") = nthreads,
//  Rcpp::Named("actual_threads_default_region") = 1,
//  Rcpp::Named("actual_threads_requested_region") = 1,
//  Rcpp::Named("omp_get_max_threads") = 1,
//  Rcpp::Named("omp_get_num_procs") = 1
// );
// #endif
// }
//
// // -----------------------------------------------------------------------------
// // Small helpers
// // -----------------------------------------------------------------------------
//
// static inline void checked_fwrite(
//   const void* data,
//   std::size_t elem_size,
//   std::size_t n_elem,
//   FILE* fs,
//   const std::string& what
// ) {
//  if (n_elem == 0) return;
//  const std::size_t got = std::fwrite(data, elem_size, n_elem, fs);
//  if (got != n_elem) throw std::runtime_error("Short write while writing " + what);
// }
//
// static inline void validate_sparse_ld_args(
//   int n,
//   int m,
//   int block_size,
//   float r2_threshold
// ) {
//  if (n <= 1) {
//   throw std::runtime_error("Need at least two samples to compute LD.");
//  }
//
//  if (m <= 0) {
//   throw std::runtime_error("Need at least one marker to compute LD.");
//  }
//
//  if (block_size <= 0) {
//   throw std::runtime_error("block_size must be positive.");
//  }
//
//  if (!std::isfinite(r2_threshold) || r2_threshold < 0.0f) {
//   throw std::runtime_error("r2_threshold must be finite and non-negative.");
//  }
// }
//
// static inline void validate_af_and_pos(
//   const std::vector<double>& af_cpp,
//   const std::vector<int>& pos_cpp,
//   int m
// ) {
//  if (static_cast<int>(af_cpp.size()) != m) {
//   throw std::runtime_error("af must have one value per selected marker.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   const double p = af_cpp[static_cast<std::size_t>(i)];
//
//   if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//    throw std::runtime_error("af contains invalid values; expected finite values in [0, 1].");
//   }
//  }
//
//  if (!pos_cpp.empty() && static_cast<int>(pos_cpp.size()) != m) {
//   throw std::runtime_error("pos_bp must have one value per selected marker.");
//  }
// }
//
// // -----------------------------------------------------------------------------
// // Decode packed block to standardized float matrix
// // -----------------------------------------------------------------------------
//
// static void decode_packed_block_float(
//   const PackedBedMatrix& G,
//   int marker_start,
//   int marker_len,
//   const double* af,
//   float* Z,
//   int nthreads
// ) {
//  const int n = G.n;
//  const std::size_t nbytes = G.nbytes;
//
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static) num_threads(nthreads)
// #endif
//  for (int ii = 0; ii < marker_len; ++ii) {
//   const int global_i = marker_start + ii;
//   const uint8_t* packed = G.row(global_i);
//   float* z = Z + static_cast<std::size_t>(ii) * n;
//
//   const double p = af[global_i];
//   const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//   float lut[4];
//   if (denom <= 0.0 || !std::isfinite(denom)) {
//    lut[0] = lut[1] = lut[2] = lut[3] = 0.0f;
//   } else {
//    lut[0] = static_cast<float>((2.0 - 2.0 * p) / denom);
//    lut[1] = 0.0f;
//    lut[2] = static_cast<float>((1.0 - 2.0 * p) / denom);
//    lut[3] = static_cast<float>((0.0 - 2.0 * p) / denom);
//   }
//
//   for (std::size_t kb = 0; kb < nbytes; ++kb) {
//    const uint8_t x = packed[kb];
//    const int jbase = static_cast<int>(kb << 2);
//
//    if (jbase + 0 < n) z[jbase + 0] = lut[(x >> 0) & 3u];
//    if (jbase + 1 < n) z[jbase + 1] = lut[(x >> 2) & 3u];
//    if (jbase + 2 < n) z[jbase + 2] = lut[(x >> 4) & 3u];
//    if (jbase + 3 < n) z[jbase + 3] = lut[(x >> 6) & 3u];
//   }
//  }
// }
//
// // -----------------------------------------------------------------------------
// // CSR objects and LD helpers
// // -----------------------------------------------------------------------------
//
// struct SparseLDCSR {
//  int m = 0;
//  std::vector<uint64_t> row_ptr;
//  std::vector<uint32_t> col_idx;
//  std::vector<float> values;
// };
//
// struct StreamingCSRResult {
//  int m = 0;
//  int n_used = 0;
//  int n_bed = 0;
//  uint64_t nnz = 0;
//  std::string row_file;
//  std::string col_file;
//  std::string val_file;
//  std::string meta_file;
// };
//
// static void compute_ld_block_sgemm(
//   const float* ZA,
//   const float* ZB,
//   int n,
//   int ma,
//   int mb,
//   float* R
// ) {
//  const float alpha = 1.0f / static_cast<float>(n - 1);
//  const float beta = 0.0f;
//
//  cblas_sgemm(
//   CblasRowMajor,
//   CblasNoTrans,
//   CblasTrans,
//   ma,
//   mb,
//   n,
//   alpha,
//   ZA,
//   n,
//   ZB,
//   n,
//   beta,
//   R,
//   mb
//  );
// }
//
// static inline bool within_window(
//   int i,
//   int j,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants
// ) {
//  if (j <= i) return false;
//  if (max_distance_variants > 0 && j - i > max_distance_variants) return false;
//  if (pos_bp && max_distance_bp > 0 && pos_bp[j] - pos_bp[i] > max_distance_bp) return false;
//  return true;
// }
//
// static inline bool block_may_overlap_window(
//   int a0,
//   int ma,
//   int b0,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants
// ) {
//  const int last_a = a0 + ma - 1;
//  if (max_distance_variants > 0 && b0 - last_a > max_distance_variants) return false;
//  if (pos_bp && max_distance_bp > 0 && pos_bp[b0] - pos_bp[last_a] > max_distance_bp) return false;
//  return true;
// }
//
// // -----------------------------------------------------------------------------
// // In-memory CSR core
// // -----------------------------------------------------------------------------
//
// static SparseLDCSR sparse_ld_packed_core(
//   const PackedBedMatrix& G,
//   const double* af,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants,
//   float r2_threshold,
//   int block_size,
//   int nthreads
// ) {
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  const int n = G.n;
//  const int m = G.m;
//  validate_sparse_ld_args(n, m, block_size, r2_threshold);
//
//  const int nb = (m + block_size - 1) / block_size;
//  std::vector<uint64_t> row_counts(static_cast<std::size_t>(m), 0);
//
//  auto count_or_fill = [&](SparseLDCSR* csr_ptr) {
//   std::vector<uint64_t> write_pos;
//   if (csr_ptr) write_pos = csr_ptr->row_ptr;
//
//   std::vector<float> ZA(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> ZB(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> R(static_cast<std::size_t>(block_size) * block_size);
//
//   for (int a_block = 0; a_block < nb; ++a_block) {
//    const int a0 = a_block * block_size;
//    const int ma = std::min(block_size, m - a0);
//    decode_packed_block_float(G, a0, ma, af, ZA.data(), nthreads);
//
//    for (int b_block = a_block; b_block < nb; ++b_block) {
//     const int b0 = b_block * block_size;
//     const int mb = std::min(block_size, m - b0);
//
//     if (!block_may_overlap_window(a0, ma, b0, pos_bp, max_distance_bp, max_distance_variants)) break;
//
//     if (b_block == a_block) {
//      compute_ld_block_sgemm(ZA.data(), ZA.data(), n, ma, ma, R.data());
// #pragma omp parallel for schedule(static) num_threads(nthreads)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       uint64_t cnt = 0;
//       uint64_t pos = csr_ptr ? write_pos[static_cast<std::size_t>(i)] : 0;
//
//       for (int ib = ia + 1; ib < ma; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * ma + ib];
//        if (r * r >= r2_threshold) {
//         if (csr_ptr) {
//          csr_ptr->col_idx[static_cast<std::size_t>(pos)] = static_cast<uint32_t>(j);
//          csr_ptr->values[static_cast<std::size_t>(pos)] = r;
//          ++pos;
//         } else {
//          ++cnt;
//         }
//        }
//       }
//
//       if (csr_ptr) write_pos[static_cast<std::size_t>(i)] = pos;
//       else row_counts[static_cast<std::size_t>(i)] += cnt;
//      }
//     } else {
//      decode_packed_block_float(G, b0, mb, af, ZB.data(), nthreads);
//      compute_ld_block_sgemm(ZA.data(), ZB.data(), n, ma, mb, R.data());
// #pragma omp parallel for schedule(static) num_threads(nthreads)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       uint64_t cnt = 0;
//       uint64_t pos = csr_ptr ? write_pos[static_cast<std::size_t>(i)] : 0;
//
//       for (int ib = 0; ib < mb; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * mb + ib];
//        if (r * r >= r2_threshold) {
//         if (csr_ptr) {
//          csr_ptr->col_idx[static_cast<std::size_t>(pos)] = static_cast<uint32_t>(j);
//          csr_ptr->values[static_cast<std::size_t>(pos)] = r;
//          ++pos;
//         } else {
//          ++cnt;
//         }
//        }
//       }
//
//       if (csr_ptr) write_pos[static_cast<std::size_t>(i)] = pos;
//       else row_counts[static_cast<std::size_t>(i)] += cnt;
//      }
//     }
//    }
//   }
//
// #ifndef NDEBUG
//   if (csr_ptr) {
//    for (int i = 0; i < m; ++i) {
//     assert(write_pos[static_cast<std::size_t>(i)] == csr_ptr->row_ptr[static_cast<std::size_t>(i) + 1]);
//    }
//   }
// #endif
//  };
//
//  count_or_fill(nullptr);
//
//  SparseLDCSR csr;
//  csr.m = m;
//  csr.row_ptr.resize(static_cast<std::size_t>(m) + 1);
//  csr.row_ptr[0] = 0;
//  for (int i = 0; i < m; ++i) {
//   csr.row_ptr[static_cast<std::size_t>(i) + 1] =
//    csr.row_ptr[static_cast<std::size_t>(i)] + row_counts[static_cast<std::size_t>(i)];
//  }
//
//  const uint64_t nnz = csr.row_ptr[static_cast<std::size_t>(m)];
//  csr.col_idx.resize(static_cast<std::size_t>(nnz));
//  csr.values.resize(static_cast<std::size_t>(nnz));
//
//  count_or_fill(&csr);
//  return csr;
// }
//
// // -----------------------------------------------------------------------------
// // Streaming CSR core
// // -----------------------------------------------------------------------------
//
// static StreamingCSRResult sparse_ld_packed_stream_to_disk_core(
//   const PackedBedMatrix& G,
//   int n_bed,
//   const double* af,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants,
//   float r2_threshold,
//   int block_size,
//   int nthreads,
//   const std::string& out_prefix
// ) {
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  const int n = G.n;
//  const int m = G.m;
//  validate_sparse_ld_args(n, m, block_size, r2_threshold);
//
//  const int nb = (m + block_size - 1) / block_size;
//
//  StreamingCSRResult res;
//  res.m = m;
//  res.n_used = n;
//  res.n_bed = n_bed;
//  res.nnz = 0;
//  res.row_file = out_prefix + ".row_ptr.u64.bin";
//  res.col_file = out_prefix + ".col_idx.u32.0based.bin";
//  res.val_file = out_prefix + ".values.f32.bin";
//  res.meta_file = out_prefix + ".meta.txt";
//
//  FILE* row_fs = std::fopen(res.row_file.c_str(), "wb");
//  FILE* col_fs = std::fopen(res.col_file.c_str(), "wb");
//  FILE* val_fs = std::fopen(res.val_file.c_str(), "wb");
//
//  if (!row_fs || !col_fs || !val_fs) {
//   if (row_fs) std::fclose(row_fs);
//   if (col_fs) std::fclose(col_fs);
//   if (val_fs) std::fclose(val_fs);
//   throw std::runtime_error("Could not open one or more CSR output files.");
//  }
//
//  try {
//   std::vector<float> ZA(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> ZB(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> R(static_cast<std::size_t>(block_size) * block_size);
//
//   uint64_t cumulative_nnz = 0;
//   checked_fwrite(&cumulative_nnz, sizeof(uint64_t), 1, row_fs, "row_ptr[0]");
//
//   for (int a_block = 0; a_block < nb; ++a_block) {
//    const int a0 = a_block * block_size;
//    const int ma = std::min(block_size, m - a0);
//
//    decode_packed_block_float(G, a0, ma, af, ZA.data(), nthreads);
//
//    std::vector<std::vector<uint32_t>> block_cols(static_cast<std::size_t>(ma));
//    std::vector<std::vector<float>> block_vals(static_cast<std::size_t>(ma));
//    for (int ia = 0; ia < ma; ++ia) {
//     block_cols[static_cast<std::size_t>(ia)].reserve(128);
//     block_vals[static_cast<std::size_t>(ia)].reserve(128);
//    }
//
//    for (int b_block = a_block; b_block < nb; ++b_block) {
//     const int b0 = b_block * block_size;
//     const int mb = std::min(block_size, m - b0);
//
//     if (!block_may_overlap_window(a0, ma, b0, pos_bp, max_distance_bp, max_distance_variants)) break;
//
//     if (b_block == a_block) {
//      compute_ld_block_sgemm(ZA.data(), ZA.data(), n, ma, ma, R.data());
// #pragma omp parallel for schedule(static) num_threads(nthreads)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//       std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//       for (int ib = ia + 1; ib < ma; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * ma + ib];
//        if (r * r >= r2_threshold) {
//         cols.push_back(static_cast<uint32_t>(j));
//         vals.push_back(r);
//        }
//       }
//      }
//     } else {
//      decode_packed_block_float(G, b0, mb, af, ZB.data(), nthreads);
//      compute_ld_block_sgemm(ZA.data(), ZB.data(), n, ma, mb, R.data());
// #pragma omp parallel for schedule(static) num_threads(nthreads)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//       std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//       for (int ib = 0; ib < mb; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * mb + ib];
//        if (r * r >= r2_threshold) {
//         cols.push_back(static_cast<uint32_t>(j));
//         vals.push_back(r);
//        }
//       }
//      }
//     }
//    }
//
//    for (int ia = 0; ia < ma; ++ia) {
//     const std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//     const std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//     checked_fwrite(cols.data(), sizeof(uint32_t), cols.size(), col_fs, "col_idx");
//     checked_fwrite(vals.data(), sizeof(float), vals.size(), val_fs, "values");
//
//     cumulative_nnz += static_cast<uint64_t>(cols.size());
//     checked_fwrite(&cumulative_nnz, sizeof(uint64_t), 1, row_fs, "row_ptr");
//    }
//   }
//
//   res.nnz = cumulative_nnz;
//
//   std::fclose(row_fs);
//   std::fclose(col_fs);
//   std::fclose(val_fs);
//   row_fs = nullptr;
//   col_fs = nullptr;
//   val_fs = nullptr;
//
//   std::ofstream meta(res.meta_file.c_str());
//   if (!meta.is_open()) throw std::runtime_error("Could not open metadata output file.");
//
//   meta << "format=sparse_ld_csr" << std::endl;
//   meta << "storage=streamed_upper_triangle" << std::endl;
//   meta << "n_bed=" << n_bed << std::endl;
//   meta << "n_used=" << n << std::endl;
//   meta << "n_samples_used=" << n << std::endl;
//   meta << "n_variants=" << m << std::endl;
//   meta << "nnz=" << static_cast<unsigned long long>(res.nnz) << std::endl;
//   meta << "triangle=upper" << std::endl;
//   meta << "diagonal=implicit_1" << std::endl;
//   meta << "row_ptr_file=" << res.row_file << std::endl;
//   meta << "col_idx_file=" << res.col_file << std::endl;
//   meta << "values_file=" << res.val_file << std::endl;
//   meta << "row_ptr_type=uint64" << std::endl;
//   meta << "col_idx_type=uint32" << std::endl;
//   meta << "values_type=float32" << std::endl;
//   meta << "index_base=0" << std::endl;
//   meta << "value=r" << std::endl;
//   meta << "r2_threshold=" << static_cast<double>(r2_threshold) << std::endl;
//   meta << "max_distance_bp=" << max_distance_bp << std::endl;
//   meta << "max_distance_variants=" << max_distance_variants << std::endl;
//   meta << "block_size=" << block_size << std::endl;
//   meta << "nthreads=" << nthreads << std::endl;
//   meta.close();
//
//  } catch (...) {
//   if (row_fs) std::fclose(row_fs);
//   if (col_fs) std::fclose(col_fs);
//   if (val_fs) std::fclose(val_fs);
//   throw;
//  }
//
//  return res;
// }
//
// // -----------------------------------------------------------------------------
// // Rcpp helpers
// // -----------------------------------------------------------------------------
//
// static Rcpp::NumericVector u64_to_numeric(const std::vector<uint64_t>& x) {
//  Rcpp::NumericVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<double>(x[i]);
//  return out;
// }
//
// static Rcpp::IntegerVector u32_to_integer_1based(const std::vector<uint32_t>& x) {
//  Rcpp::IntegerVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<int>(x[i]) + 1;
//  return out;
// }
//
// static Rcpp::NumericVector f32_to_numeric(const std::vector<float>& x) {
//  Rcpp::NumericVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<double>(x[i]);
//  return out;
// }
//
// static std::vector<std::string> copy_bed_files(Rcpp::CharacterVector bed_files) {
//  std::vector<std::string> out(bed_files.size());
//  for (int i = 0; i < bed_files.size(); ++i) out[i] = Rcpp::as<std::string>(bed_files[i]);
//  return out;
// }
//
// static std::vector<std::vector<int>> copy_int_list(Rcpp::List xlist) {
//  std::vector<std::vector<int>> out(xlist.size());
//  for (int f = 0; f < xlist.size(); ++f) {
//   Rcpp::IntegerVector x = xlist[f];
//   out[static_cast<std::size_t>(f)].resize(x.size());
//   for (int i = 0; i < x.size(); ++i) out[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)] = x[i];
//  }
//  return out;
// }
//
// static std::vector<double> flatten_af_list(Rcpp::Nullable<Rcpp::List> af) {
//  std::vector<double> out;
//  if (af.isNotNull()) {
//   Rcpp::List af_list = Rcpp::as<Rcpp::List>(af.get());
//   for (int f = 0; f < af_list.size(); ++f) {
//    Rcpp::NumericVector x = af_list[f];
//    for (int i = 0; i < x.size(); ++i) out.push_back(x[i]);
//   }
//  }
//  return out;
// }
//
// static std::vector<int> flatten_pos_list_or_empty(Rcpp::Nullable<Rcpp::List> pos_bp) {
//  std::vector<int> out;
//  if (pos_bp.isNotNull()) {
//   Rcpp::List pos_list = Rcpp::as<Rcpp::List>(pos_bp.get());
//   for (int f = 0; f < pos_list.size(); ++f) {
//    Rcpp::IntegerVector x = pos_list[f];
//    for (int i = 0; i < x.size(); ++i) out.push_back(x[i]);
//   }
//  }
//  return out;
// }
//
// static std::vector<int> copy_rows0_or_empty(Rcpp::Nullable<Rcpp::IntegerVector> rows, int n) {
//  std::vector<int> out;
//  if (rows.isNotNull()) {
//   Rcpp::IntegerVector r(rows);
//   out.resize(r.size());
//   for (int i = 0; i < r.size(); ++i) {
//    if (r[i] == NA_INTEGER) throw std::runtime_error("rows contains NA.");
//    out[static_cast<std::size_t>(i)] = r[i] - 1;
//   }
//
//   if (static_cast<int>(out.size()) == n) {
//    bool identity = true;
//    for (int i = 0; i < n; ++i) {
//     if (out[static_cast<std::size_t>(i)] != i) {
//      identity = false;
//      break;
//     }
//    }
//    if (identity) out.clear();
//   }
//  }
//  return out;
// }
//
// // -----------------------------------------------------------------------------
// // Rcpp exports
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_stream_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   std::string out_prefix,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  if (out_prefix.empty()) Rcpp::stop("out_prefix must not be empty.");
//
//  std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
//  std::vector<int> rows0 = copy_rows0_or_empty(rows, n);
//
//  const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
//  const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bed_files_cpp,
//   n,
//   rows0_ptr,
//   n_rows,
//   cls_by_file
//  );
//
//  std::vector<double> af_cpp = flatten_af_list(af);
//  const bool af_computed = af_cpp.empty();
//  if (af_computed) af_cpp = compute_af_from_packed(G);
//
//  std::vector<int> pos_cpp = flatten_pos_list_or_empty(pos_bp);
//  validate_af_and_pos(af_cpp, pos_cpp, G.m);
//  const int* pos_ptr = pos_cpp.empty() ? nullptr : pos_cpp.data();
//
//  StreamingCSRResult res = sparse_ld_packed_stream_to_disk_core(
//   G,
//   n,
//   af_cpp.data(),
//   pos_ptr,
//   max_distance_bp,
//   max_distance_variants,
//   static_cast<float>(r2_threshold),
//   block_size,
//   nthreads,
//   out_prefix
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr_file") = res.row_file,
//   Rcpp::Named("col_idx_file") = res.col_file,
//   Rcpp::Named("values_file") = res.val_file,
//   Rcpp::Named("meta_file") = res.meta_file,
//   Rcpp::Named("nrow") = res.m,
//   Rcpp::Named("ncol") = res.m,
//   Rcpp::Named("max_distance_bp") = max_distance_bp,
//   Rcpp::Named("max_distance_variants") = max_distance_variants,
//   Rcpp::Named("r2_threshold") = r2_threshold,
//   Rcpp::Named("block_size") = block_size,
//   Rcpp::Named("nthreads") = nthreads,
//   Rcpp::Named("nnz") = static_cast<double>(res.nnz),
//   Rcpp::Named("n_bed") = n,
//   Rcpp::Named("n_used") = G.n,
//   Rcpp::Named("af_computed") = af_computed,
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = 0,
//   Rcpp::Named("value") = "r"
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_to_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
//  std::vector<int> rows0 = copy_rows0_or_empty(rows, n);
//
//  const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
//  const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bed_files_cpp,
//   n,
//   rows0_ptr,
//   n_rows,
//   cls_by_file
//  );
//
//  std::vector<double> af_cpp = flatten_af_list(af);
//  const bool af_computed = af_cpp.empty();
//  if (af_computed) af_cpp = compute_af_from_packed(G);
//
//  std::vector<int> pos_cpp = flatten_pos_list_or_empty(pos_bp);
//  validate_af_and_pos(af_cpp, pos_cpp, G.m);
//  const int* pos_ptr = pos_cpp.empty() ? nullptr : pos_cpp.data();
//
//  SparseLDCSR csr = sparse_ld_packed_core(
//   G,
//   af_cpp.data(),
//   pos_ptr,
//   max_distance_bp,
//   max_distance_variants,
//   static_cast<float>(r2_threshold),
//   block_size,
//   nthreads
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr") = u64_to_numeric(csr.row_ptr),
//   Rcpp::Named("col_idx") = u32_to_integer_1based(csr.col_idx),
//   Rcpp::Named("values") = f32_to_numeric(csr.values),
//   Rcpp::Named("nrow") = G.m,
//   Rcpp::Named("ncol") = G.m,
//   Rcpp::Named("nnz") = static_cast<double>(csr.row_ptr.back()),
//   Rcpp::Named("n_bed") = n,
//   Rcpp::Named("n_used") = G.n,
//   Rcpp::Named("max_distance_bp") = max_distance_bp,
//   Rcpp::Named("max_distance_variants") = max_distance_variants,
//   Rcpp::Named("r2_threshold") = r2_threshold,
//   Rcpp::Named("block_size") = block_size,
//   Rcpp::Named("nthreads") = nthreads,
//   Rcpp::Named("af_computed") = af_computed,
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = 1,
//   Rcpp::Named("value") = "r"
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_write_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   std::string out_prefix,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  return sparseLD_stream_CSR(
//   bed_files,
//   n,
//   cls,
//   out_prefix,
//   rows,
//   af,
//   pos_bp,
//   max_distance_bp,
//   max_distance_variants,
//   r2_threshold,
//   block_size,
//   nthreads
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_read_CSR(
//   std::string prefix,
//   bool one_based = true
// ) {
//  const std::string meta_file = prefix + ".meta.txt";
//  const std::string row_file = prefix + ".row_ptr.u64.bin";
//  const std::string col_file = prefix + ".col_idx.u32.0based.bin";
//  const std::string val_file = prefix + ".values.f32.bin";
//
//  std::ifstream meta(meta_file.c_str());
//  if (!meta.is_open()) Rcpp::stop("Could not open metadata file: %s", meta_file);
//
//  int m = -1;
//  double nnz_d = -1.0;
//
//  int max_distance_bp = -1;
//  int max_distance_variants = -1;
//  double r2_threshold = R_NaN;
//  int block_size = -1;
//  int nthreads = -1;
//  int n_bed = -1;
//  int n_used = -1;
//
//  std::string line;
//  while (std::getline(meta, line)) {
//   const std::string key_m = "n_variants=";
//   const std::string key_nnz = "nnz=";
//   const std::string key_max_bp = "max_distance_bp=";
//   const std::string key_max_var = "max_distance_variants=";
//   const std::string key_r2 = "r2_threshold=";
//   const std::string key_block = "block_size=";
//   const std::string key_threads = "nthreads=";
//   const std::string key_n_bed = "n_bed=";
//   const std::string key_n_used = "n_used=";
//   const std::string key_n_samples_used = "n_samples_used=";
//
//   if (line.rfind(key_m, 0) == 0) {
//    m = std::stoi(line.substr(key_m.size()));
//   } else if (line.rfind(key_nnz, 0) == 0) {
//    nnz_d = std::stod(line.substr(key_nnz.size()));
//   } else if (line.rfind(key_max_bp, 0) == 0) {
//    max_distance_bp = std::stoi(line.substr(key_max_bp.size()));
//   } else if (line.rfind(key_max_var, 0) == 0) {
//    max_distance_variants = std::stoi(line.substr(key_max_var.size()));
//   } else if (line.rfind(key_r2, 0) == 0) {
//    r2_threshold = std::stod(line.substr(key_r2.size()));
//   } else if (line.rfind(key_block, 0) == 0) {
//    block_size = std::stoi(line.substr(key_block.size()));
//   } else if (line.rfind(key_threads, 0) == 0) {
//    nthreads = std::stoi(line.substr(key_threads.size()));
//   } else if (line.rfind(key_n_bed, 0) == 0) {
//    n_bed = std::stoi(line.substr(key_n_bed.size()));
//   } else if (line.rfind(key_n_used, 0) == 0) {
//    n_used = std::stoi(line.substr(key_n_used.size()));
//   } else if (line.rfind(key_n_samples_used, 0) == 0 && n_used < 0) {
//    n_used = std::stoi(line.substr(key_n_samples_used.size()));
//   }
//  }
//  meta.close();
//
//  if (m <= 0) Rcpp::stop("Could not read n_variants from metadata.");
//  if (nnz_d < 0.0) Rcpp::stop("Could not read nnz from metadata.");
//
//  const std::size_t nnz = static_cast<std::size_t>(nnz_d);
//
//  std::vector<uint64_t> row_ptr_u64(static_cast<std::size_t>(m) + 1);
//  std::vector<uint32_t> col_idx_u32(nnz);
//  std::vector<float> values_f32(nnz);
//
//  auto read_file = [](const std::string& path, void* data, std::size_t nbytes) {
//   FILE* fs = std::fopen(path.c_str(), "rb");
//   if (!fs) throw std::runtime_error("Could not open file: " + path);
//   const std::size_t got = std::fread(data, 1, nbytes, fs);
//   std::fclose(fs);
//   if (got != nbytes) throw std::runtime_error("Short read from file: " + path);
//  };
//
//  read_file(row_file, row_ptr_u64.data(), row_ptr_u64.size() * sizeof(uint64_t));
//  read_file(col_file, col_idx_u32.data(), col_idx_u32.size() * sizeof(uint32_t));
//  read_file(val_file, values_f32.data(), values_f32.size() * sizeof(float));
//
//  if (row_ptr_u64[0] != 0 || row_ptr_u64[static_cast<std::size_t>(m)] != static_cast<uint64_t>(nnz)) {
//   Rcpp::stop("Invalid row_ptr: expected 0-based pointer ending at nnz.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   if (row_ptr_u64[static_cast<std::size_t>(i + 1)] < row_ptr_u64[static_cast<std::size_t>(i)]) {
//    Rcpp::stop("Invalid row_ptr: row pointers are not nondecreasing.");
//   }
//  }
//
//  Rcpp::NumericVector row_ptr(row_ptr_u64.size());
//  for (std::size_t i = 0; i < row_ptr_u64.size(); ++i) {
//   row_ptr[i] = static_cast<double>(row_ptr_u64[i]);
//  }
//
//  Rcpp::IntegerVector col_idx(col_idx_u32.size());
//  for (std::size_t i = 0; i < col_idx_u32.size(); ++i) {
//   if (col_idx_u32[i] >= static_cast<uint32_t>(m)) {
//    Rcpp::stop("LD column index out of range.");
//   }
//   col_idx[i] = static_cast<int>(col_idx_u32[i]) + (one_based ? 1 : 0);
//  }
//
//  Rcpp::NumericVector values(values_f32.size());
//  for (std::size_t i = 0; i < values_f32.size(); ++i) {
//   if (!std::isfinite(values_f32[i])) {
//    Rcpp::stop("LD value contains NaN/Inf.");
//   }
//   values[i] = static_cast<double>(values_f32[i]);
//  }
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr") = row_ptr,
//   Rcpp::Named("col_idx") = col_idx,
//   Rcpp::Named("values") = values,
//   Rcpp::Named("nrow") = m,
//   Rcpp::Named("ncol") = m,
//   Rcpp::Named("nnz") = static_cast<double>(nnz),
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = one_based ? 1 : 0,
//   Rcpp::Named("value") = "r",
//   Rcpp::Named("max_distance_bp") =
//    max_distance_bp >= 0 ? Rcpp::wrap(max_distance_bp) : Rcpp::wrap(NA_INTEGER),
//     Rcpp::Named("max_distance_variants") =
//      max_distance_variants >= 0 ? Rcpp::wrap(max_distance_variants) : Rcpp::wrap(NA_INTEGER),
//       Rcpp::Named("r2_threshold") =
//        std::isfinite(r2_threshold) ? Rcpp::wrap(r2_threshold) : Rcpp::wrap(NA_REAL),
//        Rcpp::Named("block_size") =
//         block_size >= 0 ? Rcpp::wrap(block_size) : Rcpp::wrap(NA_INTEGER),
//          Rcpp::Named("nthreads") =
//           nthreads >= 0 ? Rcpp::wrap(nthreads) : Rcpp::wrap(NA_INTEGER),
//            Rcpp::Named("n_bed") =
//             n_bed >= 0 ? Rcpp::wrap(n_bed) : Rcpp::wrap(NA_INTEGER),
//              Rcpp::Named("n_used") =
//               n_used >= 0 ? Rcpp::wrap(n_used) : Rcpp::wrap(NA_INTEGER)
//  );
// }
//

// // sparse_ld_bed_core.cpp
// //
// // Compute sparse LD from one or more PLINK BED-style 2-bit genotype files
// // and write CSR to disk.
// //
// // BED coding assumed:
// //   00 -> dosage 2
// //   01 -> missing
// //   10 -> dosage 1
// //   11 -> dosage 0
// //
// // Missing genotypes are mean-imputed after standardization, i.e. z = 0.
// // LD values stored are correlations r, not r^2.
// // Streaming CSR output stores upper triangle only, with implicit diagonal 1.
//
// #include <algorithm>
// #include <cassert>
// #include <cmath>
// #include <cstdint>
// #include <cstdio>
// #include <fstream>
// #include <numeric>
// #include <stdexcept>
// #include <string>
// #include <unordered_map>
// #include <unordered_set>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// #include <Rcpp.h>
//
// #include "packed_bed.h"
//
// extern "C" {
// #include <cblas.h>
// }
//
// // -----------------------------------------------------------------------------
// // Small helpers
// // -----------------------------------------------------------------------------
//
// static inline void checked_fwrite(
//   const void* data,
//   std::size_t elem_size,
//   std::size_t n_elem,
//   FILE* fs,
//   const std::string& what
// ) {
//  if (n_elem == 0) return;
//  const std::size_t got = std::fwrite(data, elem_size, n_elem, fs);
//  if (got != n_elem) throw std::runtime_error("Short write while writing " + what);
// }
//
// static inline void validate_sparse_ld_args(
//   int n,
//   int m,
//   int block_size,
//   float r2_threshold
// ) {
//  if (n <= 1) {
//   throw std::runtime_error("Need at least two samples to compute LD.");
//  }
//
//  if (m <= 0) {
//   throw std::runtime_error("Need at least one marker to compute LD.");
//  }
//
//  if (block_size <= 0) {
//   throw std::runtime_error("block_size must be positive.");
//  }
//
//  if (!std::isfinite(r2_threshold) || r2_threshold < 0.0f) {
//   throw std::runtime_error("r2_threshold must be finite and non-negative.");
//  }
// }
//
// static inline void validate_af_and_pos(
//   const std::vector<double>& af_cpp,
//   const std::vector<int>& pos_cpp,
//   int m
// ) {
//  if (static_cast<int>(af_cpp.size()) != m) {
//   throw std::runtime_error("af must have one value per selected marker.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   const double p = af_cpp[static_cast<std::size_t>(i)];
//
//   if (!std::isfinite(p) || p < 0.0 || p > 1.0) {
//    throw std::runtime_error("af contains invalid values; expected finite values in [0, 1].");
//   }
//  }
//
//  if (!pos_cpp.empty() && static_cast<int>(pos_cpp.size()) != m) {
//   throw std::runtime_error("pos_bp must have one value per selected marker.");
//  }
// }
//
// // -----------------------------------------------------------------------------
// // Decode packed block to standardized float matrix
// // -----------------------------------------------------------------------------
//
// static void decode_packed_block_float(
//   const PackedBedMatrix& G,
//   int marker_start,
//   int marker_len,
//   const double* af,
//   float* Z
// ) {
//  const int n = G.n;
//  const std::size_t nbytes = G.nbytes;
//
// #pragma omp parallel for schedule(static)
//  for (int ii = 0; ii < marker_len; ++ii) {
//   const int global_i = marker_start + ii;
//   const uint8_t* packed = G.row(global_i);
//   float* z = Z + static_cast<std::size_t>(ii) * n;
//
//   const double p = af[global_i];
//   const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//   float lut[4];
//   if (denom <= 0.0 || !std::isfinite(denom)) {
//    lut[0] = lut[1] = lut[2] = lut[3] = 0.0f;
//   } else {
//    lut[0] = static_cast<float>((2.0 - 2.0 * p) / denom);
//    lut[1] = 0.0f;
//    lut[2] = static_cast<float>((1.0 - 2.0 * p) / denom);
//    lut[3] = static_cast<float>((0.0 - 2.0 * p) / denom);
//   }
//
//   for (std::size_t kb = 0; kb < nbytes; ++kb) {
//    const uint8_t x = packed[kb];
//    const int jbase = static_cast<int>(kb << 2);
//
//    if (jbase + 0 < n) z[jbase + 0] = lut[(x >> 0) & 3u];
//    if (jbase + 1 < n) z[jbase + 1] = lut[(x >> 2) & 3u];
//    if (jbase + 2 < n) z[jbase + 2] = lut[(x >> 4) & 3u];
//    if (jbase + 3 < n) z[jbase + 3] = lut[(x >> 6) & 3u];
//   }
//  }
// }
//
// // -----------------------------------------------------------------------------
// // CSR objects and LD helpers
// // -----------------------------------------------------------------------------
//
// struct SparseLDCSR {
//  int m = 0;
//  std::vector<uint64_t> row_ptr;
//  std::vector<uint32_t> col_idx;
//  std::vector<float> values;
// };
//
// struct StreamingCSRResult {
//  int m = 0;
//  int n_used = 0;
//  int n_bed = 0;
//  uint64_t nnz = 0;
//  std::string row_file;
//  std::string col_file;
//  std::string val_file;
//  std::string meta_file;
// };
//
// static void compute_ld_block_sgemm(
//   const float* ZA,
//   const float* ZB,
//   int n,
//   int ma,
//   int mb,
//   float* R
// ) {
//  const float alpha = 1.0f / static_cast<float>(n - 1);
//  const float beta = 0.0f;
//
//  cblas_sgemm(
//   CblasRowMajor,
//   CblasNoTrans,
//   CblasTrans,
//   ma,
//   mb,
//   n,
//   alpha,
//   ZA,
//   n,
//   ZB,
//   n,
//   beta,
//   R,
//   mb
//  );
// }
//
// static inline bool within_window(
//   int i,
//   int j,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants
// ) {
//  if (j <= i) return false;
//  if (max_distance_variants > 0 && j - i > max_distance_variants) return false;
//  if (pos_bp && max_distance_bp > 0 && pos_bp[j] - pos_bp[i] > max_distance_bp) return false;
//  return true;
// }
//
// static inline bool block_may_overlap_window(
//   int a0,
//   int ma,
//   int b0,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants
// ) {
//  const int last_a = a0 + ma - 1;
//  if (max_distance_variants > 0 && b0 - last_a > max_distance_variants) return false;
//  if (pos_bp && max_distance_bp > 0 && pos_bp[b0] - pos_bp[last_a] > max_distance_bp) return false;
//  return true;
// }
//
// // -----------------------------------------------------------------------------
// // In-memory CSR core
// // -----------------------------------------------------------------------------
//
// static SparseLDCSR sparse_ld_packed_core(
//   const PackedBedMatrix& G,
//   const double* af,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants,
//   float r2_threshold,
//   int block_size,
//   int nthreads
// ) {
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  const int n = G.n;
//  const int m = G.m;
//  validate_sparse_ld_args(n, m, block_size, r2_threshold);
//
//  const int nb = (m + block_size - 1) / block_size;
//  std::vector<uint64_t> row_counts(static_cast<std::size_t>(m), 0);
//
//  auto count_or_fill = [&](SparseLDCSR* csr_ptr) {
//   std::vector<uint64_t> write_pos;
//   if (csr_ptr) write_pos = csr_ptr->row_ptr;
//
//   std::vector<float> ZA(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> ZB(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> R(static_cast<std::size_t>(block_size) * block_size);
//
//   for (int a_block = 0; a_block < nb; ++a_block) {
//    const int a0 = a_block * block_size;
//    const int ma = std::min(block_size, m - a0);
//    decode_packed_block_float(G, a0, ma, af, ZA.data());
//
//    for (int b_block = a_block; b_block < nb; ++b_block) {
//     const int b0 = b_block * block_size;
//     const int mb = std::min(block_size, m - b0);
//
//     if (!block_may_overlap_window(a0, ma, b0, pos_bp, max_distance_bp, max_distance_variants)) break;
//
//     if (b_block == a_block) {
//      compute_ld_block_sgemm(ZA.data(), ZA.data(), n, ma, ma, R.data());
// #pragma omp parallel for schedule(static)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       uint64_t cnt = 0;
//       uint64_t pos = csr_ptr ? write_pos[static_cast<std::size_t>(i)] : 0;
//
//       for (int ib = ia + 1; ib < ma; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * ma + ib];
//        if (r * r >= r2_threshold) {
//         if (csr_ptr) {
//          csr_ptr->col_idx[static_cast<std::size_t>(pos)] = static_cast<uint32_t>(j);
//          csr_ptr->values[static_cast<std::size_t>(pos)] = r;
//          ++pos;
//         } else {
//          ++cnt;
//         }
//        }
//       }
//
//       if (csr_ptr) write_pos[static_cast<std::size_t>(i)] = pos;
//       else row_counts[static_cast<std::size_t>(i)] += cnt;
//      }
//     } else {
//      decode_packed_block_float(G, b0, mb, af, ZB.data());
//      compute_ld_block_sgemm(ZA.data(), ZB.data(), n, ma, mb, R.data());
// #pragma omp parallel for schedule(static)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       uint64_t cnt = 0;
//       uint64_t pos = csr_ptr ? write_pos[static_cast<std::size_t>(i)] : 0;
//
//       for (int ib = 0; ib < mb; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * mb + ib];
//        if (r * r >= r2_threshold) {
//         if (csr_ptr) {
//          csr_ptr->col_idx[static_cast<std::size_t>(pos)] = static_cast<uint32_t>(j);
//          csr_ptr->values[static_cast<std::size_t>(pos)] = r;
//          ++pos;
//         } else {
//          ++cnt;
//         }
//        }
//       }
//
//       if (csr_ptr) write_pos[static_cast<std::size_t>(i)] = pos;
//       else row_counts[static_cast<std::size_t>(i)] += cnt;
//      }
//     }
//    }
//   }
//
// #ifndef NDEBUG
//   if (csr_ptr) {
//    for (int i = 0; i < m; ++i) {
//     assert(write_pos[static_cast<std::size_t>(i)] == csr_ptr->row_ptr[static_cast<std::size_t>(i) + 1]);
//    }
//   }
// #endif
//  };
//
//  count_or_fill(nullptr);
//
//  SparseLDCSR csr;
//  csr.m = m;
//  csr.row_ptr.resize(static_cast<std::size_t>(m) + 1);
//  csr.row_ptr[0] = 0;
//  for (int i = 0; i < m; ++i) {
//   csr.row_ptr[static_cast<std::size_t>(i) + 1] =
//    csr.row_ptr[static_cast<std::size_t>(i)] + row_counts[static_cast<std::size_t>(i)];
//  }
//
//  const uint64_t nnz = csr.row_ptr[static_cast<std::size_t>(m)];
//  csr.col_idx.resize(static_cast<std::size_t>(nnz));
//  csr.values.resize(static_cast<std::size_t>(nnz));
//
//  count_or_fill(&csr);
//  return csr;
// }
//
// // -----------------------------------------------------------------------------
// // Streaming CSR core
// // -----------------------------------------------------------------------------
//
// static StreamingCSRResult sparse_ld_packed_stream_to_disk_core(
//   const PackedBedMatrix& G,
//   int n_bed,
//   const double* af,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants,
//   float r2_threshold,
//   int block_size,
//   int nthreads,
//   const std::string& out_prefix
// ) {
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  const int n = G.n;
//  const int m = G.m;
//  validate_sparse_ld_args(n, m, block_size, r2_threshold);
//
//  const int nb = (m + block_size - 1) / block_size;
//
//  StreamingCSRResult res;
//  res.m = m;
//  res.n_used = n;
//  res.n_bed = n_bed;
//  res.nnz = 0;
//  res.row_file = out_prefix + ".row_ptr.u64.bin";
//  res.col_file = out_prefix + ".col_idx.u32.0based.bin";
//  res.val_file = out_prefix + ".values.f32.bin";
//  res.meta_file = out_prefix + ".meta.txt";
//
//  FILE* row_fs = std::fopen(res.row_file.c_str(), "wb");
//  FILE* col_fs = std::fopen(res.col_file.c_str(), "wb");
//  FILE* val_fs = std::fopen(res.val_file.c_str(), "wb");
//
//  if (!row_fs || !col_fs || !val_fs) {
//   if (row_fs) std::fclose(row_fs);
//   if (col_fs) std::fclose(col_fs);
//   if (val_fs) std::fclose(val_fs);
//   throw std::runtime_error("Could not open one or more CSR output files.");
//  }
//
//  try {
//   std::vector<float> ZA(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> ZB(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> R(static_cast<std::size_t>(block_size) * block_size);
//
//   uint64_t cumulative_nnz = 0;
//   checked_fwrite(&cumulative_nnz, sizeof(uint64_t), 1, row_fs, "row_ptr[0]");
//
//   for (int a_block = 0; a_block < nb; ++a_block) {
//    const int a0 = a_block * block_size;
//    const int ma = std::min(block_size, m - a0);
//
//    decode_packed_block_float(G, a0, ma, af, ZA.data());
//
//    std::vector<std::vector<uint32_t>> block_cols(static_cast<std::size_t>(ma));
//    std::vector<std::vector<float>> block_vals(static_cast<std::size_t>(ma));
//    for (int ia = 0; ia < ma; ++ia) {
//     block_cols[static_cast<std::size_t>(ia)].reserve(128);
//     block_vals[static_cast<std::size_t>(ia)].reserve(128);
//    }
//
//    for (int b_block = a_block; b_block < nb; ++b_block) {
//     const int b0 = b_block * block_size;
//     const int mb = std::min(block_size, m - b0);
//
//     if (!block_may_overlap_window(a0, ma, b0, pos_bp, max_distance_bp, max_distance_variants)) break;
//
//     if (b_block == a_block) {
//      compute_ld_block_sgemm(ZA.data(), ZA.data(), n, ma, ma, R.data());
// #pragma omp parallel for schedule(static)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//       std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//       for (int ib = ia + 1; ib < ma; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * ma + ib];
//        if (r * r >= r2_threshold) {
//         cols.push_back(static_cast<uint32_t>(j));
//         vals.push_back(r);
//        }
//       }
//      }
//     } else {
//      decode_packed_block_float(G, b0, mb, af, ZB.data());
//      compute_ld_block_sgemm(ZA.data(), ZB.data(), n, ma, mb, R.data());
// #pragma omp parallel for schedule(static)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//       std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//       for (int ib = 0; ib < mb; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * mb + ib];
//        if (r * r >= r2_threshold) {
//         cols.push_back(static_cast<uint32_t>(j));
//         vals.push_back(r);
//        }
//       }
//      }
//     }
//    }
//
//    for (int ia = 0; ia < ma; ++ia) {
//     const std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//     const std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//     checked_fwrite(cols.data(), sizeof(uint32_t), cols.size(), col_fs, "col_idx");
//     checked_fwrite(vals.data(), sizeof(float), vals.size(), val_fs, "values");
//
//     cumulative_nnz += static_cast<uint64_t>(cols.size());
//     checked_fwrite(&cumulative_nnz, sizeof(uint64_t), 1, row_fs, "row_ptr");
//    }
//   }
//
//   res.nnz = cumulative_nnz;
//
//   std::fclose(row_fs);
//   std::fclose(col_fs);
//   std::fclose(val_fs);
//   row_fs = nullptr;
//   col_fs = nullptr;
//   val_fs = nullptr;
//
//   std::ofstream meta(res.meta_file.c_str());
//   if (!meta.is_open()) throw std::runtime_error("Could not open metadata output file.");
//
//   meta << "format=sparse_ld_csr" << std::endl;
//   meta << "storage=streamed_upper_triangle" << std::endl;
//   meta << "n_bed=" << n_bed << std::endl;
//   meta << "n_used=" << n << std::endl;
//   meta << "n_samples_used=" << n << std::endl;
//   meta << "n_variants=" << m << std::endl;
//   meta << "nnz=" << static_cast<unsigned long long>(res.nnz) << std::endl;
//   meta << "triangle=upper" << std::endl;
//   meta << "diagonal=implicit_1" << std::endl;
//   meta << "row_ptr_file=" << res.row_file << std::endl;
//   meta << "col_idx_file=" << res.col_file << std::endl;
//   meta << "values_file=" << res.val_file << std::endl;
//   meta << "row_ptr_type=uint64" << std::endl;
//   meta << "col_idx_type=uint32" << std::endl;
//   meta << "values_type=float32" << std::endl;
//   meta << "index_base=0" << std::endl;
//   meta << "value=r" << std::endl;
//   meta << "r2_threshold=" << static_cast<double>(r2_threshold) << std::endl;
//   meta << "max_distance_bp=" << max_distance_bp << std::endl;
//   meta << "max_distance_variants=" << max_distance_variants << std::endl;
//   meta << "block_size=" << block_size << std::endl;
//   meta << "nthreads=" << nthreads << std::endl;
//   meta.close();
//
//  } catch (...) {
//   if (row_fs) std::fclose(row_fs);
//   if (col_fs) std::fclose(col_fs);
//   if (val_fs) std::fclose(val_fs);
//   throw;
//  }
//
//  return res;
// }
//
// // -----------------------------------------------------------------------------
// // Rcpp helpers
// // -----------------------------------------------------------------------------
//
// static Rcpp::NumericVector u64_to_numeric(const std::vector<uint64_t>& x) {
//  Rcpp::NumericVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<double>(x[i]);
//  return out;
// }
//
// static Rcpp::IntegerVector u32_to_integer_1based(const std::vector<uint32_t>& x) {
//  Rcpp::IntegerVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<int>(x[i]) + 1;
//  return out;
// }
//
// static Rcpp::NumericVector f32_to_numeric(const std::vector<float>& x) {
//  Rcpp::NumericVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<double>(x[i]);
//  return out;
// }
//
// static std::vector<std::string> copy_bed_files(Rcpp::CharacterVector bed_files) {
//  std::vector<std::string> out(bed_files.size());
//  for (int i = 0; i < bed_files.size(); ++i) out[i] = Rcpp::as<std::string>(bed_files[i]);
//  return out;
// }
//
// static std::vector<std::vector<int>> copy_int_list(Rcpp::List xlist) {
//  std::vector<std::vector<int>> out(xlist.size());
//  for (int f = 0; f < xlist.size(); ++f) {
//   Rcpp::IntegerVector x = xlist[f];
//   out[static_cast<std::size_t>(f)].resize(x.size());
//   for (int i = 0; i < x.size(); ++i) out[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)] = x[i];
//  }
//  return out;
// }
//
// static std::vector<double> flatten_af_list(Rcpp::Nullable<Rcpp::List> af) {
//  std::vector<double> out;
//  if (af.isNotNull()) {
//   Rcpp::List af_list = Rcpp::as<Rcpp::List>(af.get());
//   for (int f = 0; f < af_list.size(); ++f) {
//    Rcpp::NumericVector x = af_list[f];
//    for (int i = 0; i < x.size(); ++i) out.push_back(x[i]);
//   }
//  }
//  return out;
// }
//
// static std::vector<int> flatten_pos_list_or_empty(Rcpp::Nullable<Rcpp::List> pos_bp) {
//  std::vector<int> out;
//  if (pos_bp.isNotNull()) {
//   Rcpp::List pos_list = Rcpp::as<Rcpp::List>(pos_bp.get());
//   for (int f = 0; f < pos_list.size(); ++f) {
//    Rcpp::IntegerVector x = pos_list[f];
//    for (int i = 0; i < x.size(); ++i) out.push_back(x[i]);
//   }
//  }
//  return out;
// }
//
// static std::vector<int> copy_rows0_or_empty(Rcpp::Nullable<Rcpp::IntegerVector> rows, int n) {
//  std::vector<int> out;
//  if (rows.isNotNull()) {
//   Rcpp::IntegerVector r(rows);
//   out.resize(r.size());
//   for (int i = 0; i < r.size(); ++i) {
//    if (r[i] == NA_INTEGER) throw std::runtime_error("rows contains NA.");
//    out[static_cast<std::size_t>(i)] = r[i] - 1;
//   }
//
//   if (static_cast<int>(out.size()) == n) {
//    bool identity = true;
//    for (int i = 0; i < n; ++i) {
//     if (out[static_cast<std::size_t>(i)] != i) {
//      identity = false;
//      break;
//     }
//    }
//    if (identity) out.clear();
//   }
//  }
//  return out;
// }
//
// // -----------------------------------------------------------------------------
// // Rcpp exports
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_stream_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   std::string out_prefix,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  if (out_prefix.empty()) Rcpp::stop("out_prefix must not be empty.");
//
//  std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
//  std::vector<int> rows0 = copy_rows0_or_empty(rows, n);
//
//  const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
//  const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bed_files_cpp,
//   n,
//   rows0_ptr,
//   n_rows,
//   cls_by_file
//  );
//
//  std::vector<double> af_cpp = flatten_af_list(af);
//  const bool af_computed = af_cpp.empty();
//  if (af_computed) af_cpp = compute_af_from_packed(G);
//
//  std::vector<int> pos_cpp = flatten_pos_list_or_empty(pos_bp);
//  validate_af_and_pos(af_cpp, pos_cpp, G.m);
//  const int* pos_ptr = pos_cpp.empty() ? nullptr : pos_cpp.data();
//
//  StreamingCSRResult res = sparse_ld_packed_stream_to_disk_core(
//   G,
//   n,
//   af_cpp.data(),
//   pos_ptr,
//   max_distance_bp,
//   max_distance_variants,
//   static_cast<float>(r2_threshold),
//   block_size,
//   nthreads,
//   out_prefix
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr_file") = res.row_file,
//   Rcpp::Named("col_idx_file") = res.col_file,
//   Rcpp::Named("values_file") = res.val_file,
//   Rcpp::Named("meta_file") = res.meta_file,
//   Rcpp::Named("nrow") = res.m,
//   Rcpp::Named("ncol") = res.m,
//   Rcpp::Named("max_distance_bp") = max_distance_bp,
//   Rcpp::Named("max_distance_variants") = max_distance_variants,
//   Rcpp::Named("r2_threshold") = r2_threshold,
//   Rcpp::Named("block_size") = block_size,
//   Rcpp::Named("nthreads") = nthreads,
//   Rcpp::Named("nnz") = static_cast<double>(res.nnz),
//   Rcpp::Named("n_bed") = n,
//   Rcpp::Named("n_used") = G.n,
//   Rcpp::Named("af_computed") = af_computed,
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = 0,
//   Rcpp::Named("value") = "r"
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_to_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
//  std::vector<int> rows0 = copy_rows0_or_empty(rows, n);
//
//  const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
//  const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bed_files_cpp,
//   n,
//   rows0_ptr,
//   n_rows,
//   cls_by_file
//  );
//
//  std::vector<double> af_cpp = flatten_af_list(af);
//  const bool af_computed = af_cpp.empty();
//  if (af_computed) af_cpp = compute_af_from_packed(G);
//
//  std::vector<int> pos_cpp = flatten_pos_list_or_empty(pos_bp);
//  validate_af_and_pos(af_cpp, pos_cpp, G.m);
//  const int* pos_ptr = pos_cpp.empty() ? nullptr : pos_cpp.data();
//
//  SparseLDCSR csr = sparse_ld_packed_core(
//   G,
//   af_cpp.data(),
//   pos_ptr,
//   max_distance_bp,
//   max_distance_variants,
//   static_cast<float>(r2_threshold),
//   block_size,
//   nthreads
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr") = u64_to_numeric(csr.row_ptr),
//   Rcpp::Named("col_idx") = u32_to_integer_1based(csr.col_idx),
//   Rcpp::Named("values") = f32_to_numeric(csr.values),
//   Rcpp::Named("nrow") = G.m,
//   Rcpp::Named("ncol") = G.m,
//   Rcpp::Named("nnz") = static_cast<double>(csr.row_ptr.back()),
//   Rcpp::Named("n_bed") = n,
//   Rcpp::Named("n_used") = G.n,
//   Rcpp::Named("max_distance_bp") = max_distance_bp,
//   Rcpp::Named("max_distance_variants") = max_distance_variants,
//   Rcpp::Named("r2_threshold") = r2_threshold,
//   Rcpp::Named("block_size") = block_size,
//   Rcpp::Named("nthreads") = nthreads,
//   Rcpp::Named("af_computed") = af_computed,
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = 1,
//   Rcpp::Named("value") = "r"
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_write_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   std::string out_prefix,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  return sparseLD_stream_CSR(
//   bed_files,
//   n,
//   cls,
//   out_prefix,
//   rows,
//   af,
//   pos_bp,
//   max_distance_bp,
//   max_distance_variants,
//   r2_threshold,
//   block_size,
//   nthreads
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_read_CSR(
//   std::string prefix,
//   bool one_based = true
// ) {
//  const std::string meta_file = prefix + ".meta.txt";
//  const std::string row_file = prefix + ".row_ptr.u64.bin";
//  const std::string col_file = prefix + ".col_idx.u32.0based.bin";
//  const std::string val_file = prefix + ".values.f32.bin";
//
//  std::ifstream meta(meta_file.c_str());
//  if (!meta.is_open()) Rcpp::stop("Could not open metadata file: %s", meta_file);
//
//  int m = -1;
//  double nnz_d = -1.0;
//
//  int max_distance_bp = -1;
//  int max_distance_variants = -1;
//  double r2_threshold = R_NaN;
//  int block_size = -1;
//  int nthreads = -1;
//  int n_bed = -1;
//  int n_used = -1;
//
//  std::string line;
//  while (std::getline(meta, line)) {
//   const std::string key_m = "n_variants=";
//   const std::string key_nnz = "nnz=";
//   const std::string key_max_bp = "max_distance_bp=";
//   const std::string key_max_var = "max_distance_variants=";
//   const std::string key_r2 = "r2_threshold=";
//   const std::string key_block = "block_size=";
//   const std::string key_threads = "nthreads=";
//   const std::string key_n_bed = "n_bed=";
//   const std::string key_n_used = "n_used=";
//   const std::string key_n_samples_used = "n_samples_used=";
//
//   if (line.rfind(key_m, 0) == 0) {
//    m = std::stoi(line.substr(key_m.size()));
//   } else if (line.rfind(key_nnz, 0) == 0) {
//    nnz_d = std::stod(line.substr(key_nnz.size()));
//   } else if (line.rfind(key_max_bp, 0) == 0) {
//    max_distance_bp = std::stoi(line.substr(key_max_bp.size()));
//   } else if (line.rfind(key_max_var, 0) == 0) {
//    max_distance_variants = std::stoi(line.substr(key_max_var.size()));
//   } else if (line.rfind(key_r2, 0) == 0) {
//    r2_threshold = std::stod(line.substr(key_r2.size()));
//   } else if (line.rfind(key_block, 0) == 0) {
//    block_size = std::stoi(line.substr(key_block.size()));
//   } else if (line.rfind(key_threads, 0) == 0) {
//    nthreads = std::stoi(line.substr(key_threads.size()));
//   } else if (line.rfind(key_n_bed, 0) == 0) {
//    n_bed = std::stoi(line.substr(key_n_bed.size()));
//   } else if (line.rfind(key_n_used, 0) == 0) {
//    n_used = std::stoi(line.substr(key_n_used.size()));
//   } else if (line.rfind(key_n_samples_used, 0) == 0 && n_used < 0) {
//    n_used = std::stoi(line.substr(key_n_samples_used.size()));
//   }
//  }
//  meta.close();
//
//  if (m <= 0) Rcpp::stop("Could not read n_variants from metadata.");
//  if (nnz_d < 0.0) Rcpp::stop("Could not read nnz from metadata.");
//
//  const std::size_t nnz = static_cast<std::size_t>(nnz_d);
//
//  std::vector<uint64_t> row_ptr_u64(static_cast<std::size_t>(m) + 1);
//  std::vector<uint32_t> col_idx_u32(nnz);
//  std::vector<float> values_f32(nnz);
//
//  auto read_file = [](const std::string& path, void* data, std::size_t nbytes) {
//   FILE* fs = std::fopen(path.c_str(), "rb");
//   if (!fs) throw std::runtime_error("Could not open file: " + path);
//   const std::size_t got = std::fread(data, 1, nbytes, fs);
//   std::fclose(fs);
//   if (got != nbytes) throw std::runtime_error("Short read from file: " + path);
//  };
//
//  read_file(row_file, row_ptr_u64.data(), row_ptr_u64.size() * sizeof(uint64_t));
//  read_file(col_file, col_idx_u32.data(), col_idx_u32.size() * sizeof(uint32_t));
//  read_file(val_file, values_f32.data(), values_f32.size() * sizeof(float));
//
//  if (row_ptr_u64[0] != 0 || row_ptr_u64[static_cast<std::size_t>(m)] != static_cast<uint64_t>(nnz)) {
//   Rcpp::stop("Invalid row_ptr: expected 0-based pointer ending at nnz.");
//  }
//
//  for (int i = 0; i < m; ++i) {
//   if (row_ptr_u64[static_cast<std::size_t>(i + 1)] < row_ptr_u64[static_cast<std::size_t>(i)]) {
//    Rcpp::stop("Invalid row_ptr: row pointers are not nondecreasing.");
//   }
//  }
//
//  Rcpp::NumericVector row_ptr(row_ptr_u64.size());
//  for (std::size_t i = 0; i < row_ptr_u64.size(); ++i) {
//   row_ptr[i] = static_cast<double>(row_ptr_u64[i]);
//  }
//
//  Rcpp::IntegerVector col_idx(col_idx_u32.size());
//  for (std::size_t i = 0; i < col_idx_u32.size(); ++i) {
//   if (col_idx_u32[i] >= static_cast<uint32_t>(m)) {
//    Rcpp::stop("LD column index out of range.");
//   }
//   col_idx[i] = static_cast<int>(col_idx_u32[i]) + (one_based ? 1 : 0);
//  }
//
//  Rcpp::NumericVector values(values_f32.size());
//  for (std::size_t i = 0; i < values_f32.size(); ++i) {
//   if (!std::isfinite(values_f32[i])) {
//    Rcpp::stop("LD value contains NaN/Inf.");
//   }
//   values[i] = static_cast<double>(values_f32[i]);
//  }
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr") = row_ptr,
//   Rcpp::Named("col_idx") = col_idx,
//   Rcpp::Named("values") = values,
//   Rcpp::Named("nrow") = m,
//   Rcpp::Named("ncol") = m,
//   Rcpp::Named("nnz") = static_cast<double>(nnz),
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = one_based ? 1 : 0,
//   Rcpp::Named("value") = "r",
//   Rcpp::Named("max_distance_bp") =
//    max_distance_bp >= 0 ? Rcpp::wrap(max_distance_bp) : Rcpp::wrap(NA_INTEGER),
//     Rcpp::Named("max_distance_variants") =
//      max_distance_variants >= 0 ? Rcpp::wrap(max_distance_variants) : Rcpp::wrap(NA_INTEGER),
//       Rcpp::Named("r2_threshold") =
//        std::isfinite(r2_threshold) ? Rcpp::wrap(r2_threshold) : Rcpp::wrap(NA_REAL),
//        Rcpp::Named("block_size") =
//         block_size >= 0 ? Rcpp::wrap(block_size) : Rcpp::wrap(NA_INTEGER),
//          Rcpp::Named("nthreads") =
//           nthreads >= 0 ? Rcpp::wrap(nthreads) : Rcpp::wrap(NA_INTEGER),
//            Rcpp::Named("n_bed") =
//             n_bed >= 0 ? Rcpp::wrap(n_bed) : Rcpp::wrap(NA_INTEGER),
//              Rcpp::Named("n_used") =
//               n_used >= 0 ? Rcpp::wrap(n_used) : Rcpp::wrap(NA_INTEGER)
//  );
// }

// =============================================================================
// Block dependency graph + greedy coloring from sparse LD CSR
// =============================================================================
//
// Purpose:
//   Given marker-to-block assignments ("sets") and a sparse LD matrix, build a
//   graph where two blocks are connected if any marker pair across the blocks has
//   |r| >= ld_abs_threshold or r^2 >= ld_r2_threshold.
//
// Output:
//   block_color: length n_blocks, 1-based color per block
//   marker_block: length m, 1-based compact block index per marker
//   blocks_by_color: list of block indices per color, 1-based
//   n_colors, n_blocks, n_edges
//
// Notes:
//   - sets can use arbitrary labels, e.g. 1,1,1,2,2,10,10.
//   - output uses compact block IDs 1:n_blocks in first-occurrence order.
//   - ld_row_ptr may be 0-based or 1-based.
//   - ld_col_idx may be 0-based or 1-based, controlled by ld_col_idx_one_based.
//   - works with upper-triangle sparse LD, symmetric sparse LD, or mixed.
// =============================================================================

// [[Rcpp::export]]
Rcpp::List stblr_color_blocks_from_ld(
  std::vector<int> sets,
  std::vector<double> ld_row_ptr,
  std::vector<int> ld_col_idx,
  std::vector<double> ld_values,
  bool ld_col_idx_one_based = true,
  double ld_abs_threshold = 0.01,
  double ld_r2_threshold = 0.0
) {
 const int m = static_cast<int>(sets.size());

 if (m <= 0) {
  throw std::runtime_error("stblr_color_blocks_from_ld: sets must be non-empty.");
 }

 if (static_cast<int>(ld_row_ptr.size()) != m + 1) {
  throw std::runtime_error("stblr_color_blocks_from_ld: ld_row_ptr must have length m + 1.");
 }

 if (ld_col_idx.size() != ld_values.size()) {
  throw std::runtime_error("stblr_color_blocks_from_ld: ld_col_idx and ld_values must have same length.");
 }

 if (ld_abs_threshold < 0.0 || ld_r2_threshold < 0.0) {
  throw std::runtime_error("stblr_color_blocks_from_ld: LD thresholds must be non-negative.");
 }

 const std::size_t nnz = ld_values.size();

 int row_ptr_base = -1;

 if (ld_row_ptr[0] == 0.0 && ld_row_ptr[m] == static_cast<double>(nnz)) {
  row_ptr_base = 0;
 } else if (ld_row_ptr[0] == 1.0 && ld_row_ptr[m] == static_cast<double>(nnz + 1)) {
  row_ptr_base = 1;
 } else {
  throw std::runtime_error(
    "stblr_color_blocks_from_ld: row_ptr base not recognized. "
    "Expected row_ptr[0]=0,row_ptr[m]=nnz or row_ptr[0]=1,row_ptr[m]=nnz+1."
  );
 }

 // --------------------------------------------------------------------------
 // Compact block labels in first-occurrence order.
 // --------------------------------------------------------------------------

 std::vector<int> block_labels;
 block_labels.reserve(m);

 std::unordered_map<int, int> label_to_block;
 std::vector<int> marker_block0(m, -1);  // 0-based compact block index

 for (int i = 0; i < m; ++i) {
  const int lab = sets[static_cast<std::size_t>(i)];

  auto it = label_to_block.find(lab);

  if (it == label_to_block.end()) {
   const int b = static_cast<int>(block_labels.size());
   label_to_block[lab] = b;
   block_labels.push_back(lab);
   marker_block0[static_cast<std::size_t>(i)] = b;
  } else {
   marker_block0[static_cast<std::size_t>(i)] = it->second;
  }
 }

 const int n_blocks = static_cast<int>(block_labels.size());

 if (n_blocks <= 0) {
  throw std::runtime_error("stblr_color_blocks_from_ld: no blocks found.");
 }

 // --------------------------------------------------------------------------
 // Build block graph.
 // --------------------------------------------------------------------------

 std::vector<std::unordered_set<int>> graph(
   static_cast<std::size_t>(n_blocks)
 );

 std::size_t n_ld_edges_used = 0;

 for (int i = 0; i < m; ++i) {
  const long long start_ll =
   static_cast<long long>(ld_row_ptr[static_cast<std::size_t>(i)]) - row_ptr_base;

  const long long end_ll =
   static_cast<long long>(ld_row_ptr[static_cast<std::size_t>(i + 1)]) - row_ptr_base;

  if (start_ll < 0 || end_ll < start_ll ||
      static_cast<std::size_t>(end_ll) > nnz) {
   throw std::runtime_error(
     "stblr_color_blocks_from_ld: invalid row_ptr slice at row " +
      std::to_string(i)
   );
  }

  const int bi = marker_block0[static_cast<std::size_t>(i)];

  for (std::size_t kk = static_cast<std::size_t>(start_ll);
       kk < static_cast<std::size_t>(end_ll);
       ++kk) {
   int j = ld_col_idx[kk];

   if (ld_col_idx_one_based) {
    --j;
   }

   if (j < 0 || j >= m) {
    throw std::runtime_error(
      "stblr_color_blocks_from_ld: LD column index out of range at row " +
       std::to_string(i)
    );
   }

   if (j == i) {
    continue;
   }

   const double r = ld_values[kk];

   if (!std::isfinite(r)) {
    throw std::runtime_error(
      "stblr_color_blocks_from_ld: LD value contains NaN/Inf at row " +
       std::to_string(i)
    );
   }

   const double abs_r = std::abs(r);
   const double r2 = r * r;

   const bool linked =
    (ld_abs_threshold > 0.0 && abs_r >= ld_abs_threshold) ||
    (ld_r2_threshold > 0.0 && r2 >= ld_r2_threshold);

   if (!linked) {
    continue;
   }

   const int bj = marker_block0[static_cast<std::size_t>(j)];

   if (bi == bj) {
    continue;
   }

   graph[static_cast<std::size_t>(bi)].insert(bj);
   graph[static_cast<std::size_t>(bj)].insert(bi);

   ++n_ld_edges_used;
  }
 }

 // Count unique block edges.
 std::size_t n_block_edges = 0;

 for (int b = 0; b < n_blocks; ++b) {
  for (const int nb : graph[static_cast<std::size_t>(b)]) {
   if (nb > b) {
    ++n_block_edges;
   }
  }
 }

 // --------------------------------------------------------------------------
 // Greedy coloring.
 //
 // Heuristic: color high-degree blocks first, then map colors back to original
 // block order.
 // --------------------------------------------------------------------------

 std::vector<int> block_order(n_blocks);
 std::iota(block_order.begin(), block_order.end(), 0);

 std::sort(block_order.begin(), block_order.end(), [&](int a, int b) {
  const std::size_t da = graph[static_cast<std::size_t>(a)].size();
  const std::size_t db = graph[static_cast<std::size_t>(b)].size();

  if (da != db) {
   return da > db;
  }

  return a < b;
 });

 std::vector<int> color0(n_blocks, -1);  // 0-based color

 for (const int b : block_order) {
  std::unordered_set<int> used_colors;

  for (const int nb : graph[static_cast<std::size_t>(b)]) {
   const int c = color0[static_cast<std::size_t>(nb)];

   if (c >= 0) {
    used_colors.insert(c);
   }
  }

  int c = 0;

  while (used_colors.find(c) != used_colors.end()) {
   ++c;
  }

  color0[static_cast<std::size_t>(b)] = c;
 }

 int n_colors = 0;

 for (int b = 0; b < n_blocks; ++b) {
  n_colors = std::max(n_colors, color0[static_cast<std::size_t>(b)] + 1);
 }

 // --------------------------------------------------------------------------
 // Build R outputs.
 // --------------------------------------------------------------------------

 Rcpp::IntegerVector marker_block(m);

 for (int i = 0; i < m; ++i) {
  marker_block[i] = marker_block0[static_cast<std::size_t>(i)] + 1;
 }

 Rcpp::IntegerVector block_label_out(n_blocks);
 Rcpp::IntegerVector block_color(n_blocks);
 Rcpp::IntegerVector block_degree(n_blocks);

 for (int b = 0; b < n_blocks; ++b) {
  block_label_out[b] = block_labels[static_cast<std::size_t>(b)];
  block_color[b] = color0[static_cast<std::size_t>(b)] + 1;
  block_degree[b] = static_cast<int>(graph[static_cast<std::size_t>(b)].size());
 }

 Rcpp::List blocks_by_color(n_colors);

 for (int c = 0; c < n_colors; ++c) {
  std::vector<int> blocks_c;

  for (int b = 0; b < n_blocks; ++b) {
   if (color0[static_cast<std::size_t>(b)] == c) {
    blocks_c.push_back(b + 1);  // 1-based compact block index
   }
  }

  blocks_by_color[c] = Rcpp::wrap(blocks_c);
 }

 Rcpp::IntegerVector color_sizes(n_colors);

 for (int c = 0; c < n_colors; ++c) {
  Rcpp::IntegerVector bc = blocks_by_color[c];
  color_sizes[c] = bc.size();
 }

 return Rcpp::List::create(
  Rcpp::Named("marker_block") = marker_block,
  Rcpp::Named("block_labels") = block_label_out,
  Rcpp::Named("block_color") = block_color,
  Rcpp::Named("block_degree") = block_degree,
  Rcpp::Named("blocks_by_color") = blocks_by_color,
  Rcpp::Named("color_sizes") = color_sizes,
  Rcpp::Named("n_blocks") = n_blocks,
  Rcpp::Named("n_colors") = n_colors,
  Rcpp::Named("n_block_edges") = static_cast<double>(n_block_edges),
  Rcpp::Named("n_ld_edges_used") = static_cast<double>(n_ld_edges_used),
  Rcpp::Named("ld_abs_threshold") = ld_abs_threshold,
  Rcpp::Named("ld_r2_threshold") = ld_r2_threshold,
  Rcpp::Named("ld_col_idx_one_based") = ld_col_idx_one_based
 );
}


// // sparse_ld_bed_core.cpp
// //
// // Compute sparse LD from one or more PLINK BED-style 2-bit genotype files
// // and write CSR to disk.
// //
// // BED coding assumed:
// //   00 -> dosage 2
// //   01 -> missing
// //   10 -> dosage 1
// //   11 -> dosage 0
// //
// // Missing genotypes are mean-imputed after standardization, i.e. z = 0.
// // LD values stored are correlations r, not r^2.
// // Streaming CSR output stores upper triangle only, with implicit diagonal 1.
//
// #include <algorithm>
// #include <cassert>
// #include <cmath>
// #include <cstdint>
// #include <cstdio>
// #include <fstream>
// #include <stdexcept>
// #include <string>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// #include <Rcpp.h>
//
// #include "packed_bed.h"
//
// extern "C" {
// #include <cblas.h>
// }
//
// // -----------------------------------------------------------------------------
// // Small helpers
// // -----------------------------------------------------------------------------
//
//
// static inline void checked_fwrite(
//   const void* data,
//   std::size_t elem_size,
//   std::size_t n_elem,
//   FILE* fs,
//   const std::string& what
// ) {
//  if (n_elem == 0) return;
//  const std::size_t got = std::fwrite(data, elem_size, n_elem, fs);
//  if (got != n_elem) throw std::runtime_error("Short write while writing " + what);
// }
//
//
//
// // -----------------------------------------------------------------------------
// // Decode packed block to standardized float matrix
// // -----------------------------------------------------------------------------
//
// static void decode_packed_block_float(
//   const PackedBedMatrix& G,
//   int marker_start,
//   int marker_len,
//   const double* af,
//   float* Z
// ) {
//  const int n = G.n;
//  const std::size_t nbytes = G.nbytes;
//
// #pragma omp parallel for schedule(static)
//  for (int ii = 0; ii < marker_len; ++ii) {
//   const int global_i = marker_start + ii;
//   const uint8_t* packed = G.row(global_i);
//   float* z = Z + static_cast<std::size_t>(ii) * n;
//
//   const double p = af[global_i];
//   const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//   float lut[4];
//   if (denom <= 0.0 || !std::isfinite(denom)) {
//    lut[0] = lut[1] = lut[2] = lut[3] = 0.0f;
//   } else {
//    lut[0] = static_cast<float>((2.0 - 2.0 * p) / denom);
//    lut[1] = 0.0f;
//    lut[2] = static_cast<float>((1.0 - 2.0 * p) / denom);
//    lut[3] = static_cast<float>((0.0 - 2.0 * p) / denom);
//   }
//
//   for (std::size_t kb = 0; kb < nbytes; ++kb) {
//    const uint8_t x = packed[kb];
//    const int jbase = static_cast<int>(kb << 2);
//
//    if (jbase + 0 < n) z[jbase + 0] = lut[(x >> 0) & 3u];
//    if (jbase + 1 < n) z[jbase + 1] = lut[(x >> 2) & 3u];
//    if (jbase + 2 < n) z[jbase + 2] = lut[(x >> 4) & 3u];
//    if (jbase + 3 < n) z[jbase + 3] = lut[(x >> 6) & 3u];
//   }
//  }
// }
//
// // -----------------------------------------------------------------------------
// // CSR objects and LD helpers
// // -----------------------------------------------------------------------------
//
// struct SparseLDCSR {
//  int m = 0;
//  std::vector<uint64_t> row_ptr;
//  std::vector<uint32_t> col_idx;
//  std::vector<float> values;
// };
//
// struct StreamingCSRResult {
//  int m = 0;
//  uint64_t nnz = 0;
//  std::string row_file;
//  std::string col_file;
//  std::string val_file;
//  std::string meta_file;
// };
//
// static void compute_ld_block_sgemm(
//   const float* ZA,
//   const float* ZB,
//   int n,
//   int ma,
//   int mb,
//   float* R
// ) {
//  const float alpha = 1.0f / static_cast<float>(n - 1);
//  const float beta = 0.0f;
//
//  cblas_sgemm(
//   CblasRowMajor,
//   CblasNoTrans,
//   CblasTrans,
//   ma,
//   mb,
//   n,
//   alpha,
//   ZA,
//   n,
//   ZB,
//   n,
//   beta,
//   R,
//   mb
//  );
// }
//
// static inline bool within_window(
//   int i,
//   int j,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants
// ) {
//  if (j <= i) return false;
//  if (max_distance_variants > 0 && j - i > max_distance_variants) return false;
//  if (pos_bp && max_distance_bp > 0 && pos_bp[j] - pos_bp[i] > max_distance_bp) return false;
//  return true;
// }
//
// static inline bool block_may_overlap_window(
//   int a0,
//   int ma,
//   int b0,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants
// ) {
//  const int last_a = a0 + ma - 1;
//  if (max_distance_variants > 0 && b0 - last_a > max_distance_variants) return false;
//  if (pos_bp && max_distance_bp > 0 && pos_bp[b0] - pos_bp[last_a] > max_distance_bp) return false;
//  return true;
// }
//
// // -----------------------------------------------------------------------------
// // In-memory CSR core
// // -----------------------------------------------------------------------------
//
// static SparseLDCSR sparse_ld_packed_core(
//   const PackedBedMatrix& G,
//   const double* af,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants,
//   float r2_threshold,
//   int block_size,
//   int nthreads
// ) {
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  const int n = G.n;
//  const int m = G.m;
//  const int nb = (m + block_size - 1) / block_size;
//
//  std::vector<uint64_t> row_counts(static_cast<std::size_t>(m), 0);
//
//  auto count_or_fill = [&](SparseLDCSR* csr_ptr) {
//   std::vector<uint64_t> write_pos;
//   if (csr_ptr) write_pos = csr_ptr->row_ptr;
//
//   std::vector<float> ZA(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> ZB(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> R(static_cast<std::size_t>(block_size) * block_size);
//
//   for (int a_block = 0; a_block < nb; ++a_block) {
//    const int a0 = a_block * block_size;
//    const int ma = std::min(block_size, m - a0);
//    decode_packed_block_float(G, a0, ma, af, ZA.data());
//
//    for (int b_block = a_block; b_block < nb; ++b_block) {
//     const int b0 = b_block * block_size;
//     const int mb = std::min(block_size, m - b0);
//
//     if (!block_may_overlap_window(a0, ma, b0, pos_bp, max_distance_bp, max_distance_variants)) break;
//
//     if (b_block == a_block) {
//      compute_ld_block_sgemm(ZA.data(), ZA.data(), n, ma, ma, R.data());
// #pragma omp parallel for schedule(static)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       uint64_t cnt = 0;
//       uint64_t pos = csr_ptr ? write_pos[i] : 0;
//
//       for (int ib = ia + 1; ib < ma; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * ma + ib];
//        if (r * r >= r2_threshold) {
//         if (csr_ptr) {
//          csr_ptr->col_idx[static_cast<std::size_t>(pos)] = static_cast<uint32_t>(j);
//          csr_ptr->values[static_cast<std::size_t>(pos)] = r;
//          ++pos;
//         } else {
//          ++cnt;
//         }
//        }
//       }
//
//       if (csr_ptr) write_pos[i] = pos;
//       else row_counts[i] += cnt;
//      }
//     } else {
//      decode_packed_block_float(G, b0, mb, af, ZB.data());
//      compute_ld_block_sgemm(ZA.data(), ZB.data(), n, ma, mb, R.data());
// #pragma omp parallel for schedule(static)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       uint64_t cnt = 0;
//       uint64_t pos = csr_ptr ? write_pos[i] : 0;
//
//       for (int ib = 0; ib < mb; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * mb + ib];
//        if (r * r >= r2_threshold) {
//         if (csr_ptr) {
//          csr_ptr->col_idx[static_cast<std::size_t>(pos)] = static_cast<uint32_t>(j);
//          csr_ptr->values[static_cast<std::size_t>(pos)] = r;
//          ++pos;
//         } else {
//          ++cnt;
//         }
//        }
//       }
//
//       if (csr_ptr) write_pos[i] = pos;
//       else row_counts[i] += cnt;
//      }
//     }
//    }
//   }
//
// #ifndef NDEBUG
//   if (csr_ptr) {
//    for (int i = 0; i < m; ++i) {
//     assert(write_pos[i] == csr_ptr->row_ptr[static_cast<std::size_t>(i) + 1]);
//    }
//   }
// #endif
//  };
//
//  count_or_fill(nullptr);
//
//  SparseLDCSR csr;
//  csr.m = m;
//  csr.row_ptr.resize(static_cast<std::size_t>(m) + 1);
//  csr.row_ptr[0] = 0;
//  for (int i = 0; i < m; ++i) {
//   csr.row_ptr[static_cast<std::size_t>(i) + 1] = csr.row_ptr[i] + row_counts[i];
//  }
//
//  const uint64_t nnz = csr.row_ptr[m];
//  csr.col_idx.resize(static_cast<std::size_t>(nnz));
//  csr.values.resize(static_cast<std::size_t>(nnz));
//
//  count_or_fill(&csr);
//  return csr;
// }
//
// // -----------------------------------------------------------------------------
// // Streaming CSR core
// // -----------------------------------------------------------------------------
//
// static StreamingCSRResult sparse_ld_packed_stream_to_disk_core(
//   const PackedBedMatrix& G,
//   const double* af,
//   const int* pos_bp,
//   int max_distance_bp,
//   int max_distance_variants,
//   float r2_threshold,
//   int block_size,
//   int nthreads,
//   const std::string& out_prefix
// ) {
// #ifdef _OPENMP
//  if (nthreads > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//  }
// #endif
//
//  const int n = G.n;
//  const int m = G.m;
//  const int nb = (m + block_size - 1) / block_size;
//
//  StreamingCSRResult res;
//  res.m = m;
//  res.nnz = 0;
//  res.row_file = out_prefix + ".row_ptr.u64.bin";
//  res.col_file = out_prefix + ".col_idx.u32.0based.bin";
//  res.val_file = out_prefix + ".values.f32.bin";
//  res.meta_file = out_prefix + ".meta.txt";
//
//  FILE* row_fs = std::fopen(res.row_file.c_str(), "wb");
//  FILE* col_fs = std::fopen(res.col_file.c_str(), "wb");
//  FILE* val_fs = std::fopen(res.val_file.c_str(), "wb");
//
//  if (!row_fs || !col_fs || !val_fs) {
//   if (row_fs) std::fclose(row_fs);
//   if (col_fs) std::fclose(col_fs);
//   if (val_fs) std::fclose(val_fs);
//   throw std::runtime_error("Could not open one or more CSR output files.");
//  }
//
//  try {
//   std::vector<float> ZA(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> ZB(static_cast<std::size_t>(block_size) * n);
//   std::vector<float> R(static_cast<std::size_t>(block_size) * block_size);
//
//   uint64_t cumulative_nnz = 0;
//   checked_fwrite(&cumulative_nnz, sizeof(uint64_t), 1, row_fs, "row_ptr[0]");
//
//   for (int a_block = 0; a_block < nb; ++a_block) {
//    const int a0 = a_block * block_size;
//    const int ma = std::min(block_size, m - a0);
//
//    decode_packed_block_float(G, a0, ma, af, ZA.data());
//
//    std::vector<std::vector<uint32_t>> block_cols(static_cast<std::size_t>(ma));
//    std::vector<std::vector<float>> block_vals(static_cast<std::size_t>(ma));
//    for (int ia = 0; ia < ma; ++ia) {
//     block_cols[static_cast<std::size_t>(ia)].reserve(128);
//     block_vals[static_cast<std::size_t>(ia)].reserve(128);
//    }
//
//    for (int b_block = a_block; b_block < nb; ++b_block) {
//     const int b0 = b_block * block_size;
//     const int mb = std::min(block_size, m - b0);
//
//     if (!block_may_overlap_window(a0, ma, b0, pos_bp, max_distance_bp, max_distance_variants)) break;
//
//     if (b_block == a_block) {
//      compute_ld_block_sgemm(ZA.data(), ZA.data(), n, ma, ma, R.data());
// #pragma omp parallel for schedule(static)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//       std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//       for (int ib = ia + 1; ib < ma; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * ma + ib];
//        if (r * r >= r2_threshold) {
//         cols.push_back(static_cast<uint32_t>(j));
//         vals.push_back(r);
//        }
//       }
//      }
//     } else {
//      decode_packed_block_float(G, b0, mb, af, ZB.data());
//      compute_ld_block_sgemm(ZA.data(), ZB.data(), n, ma, mb, R.data());
// #pragma omp parallel for schedule(static)
//      for (int ia = 0; ia < ma; ++ia) {
//       const int i = a0 + ia;
//       std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//       std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//       for (int ib = 0; ib < mb; ++ib) {
//        const int j = b0 + ib;
//        if (!within_window(i, j, pos_bp, max_distance_bp, max_distance_variants)) continue;
//        const float r = R[static_cast<std::size_t>(ia) * mb + ib];
//        if (r * r >= r2_threshold) {
//         cols.push_back(static_cast<uint32_t>(j));
//         vals.push_back(r);
//        }
//       }
//      }
//     }
//    }
//
//    for (int ia = 0; ia < ma; ++ia) {
//     const std::vector<uint32_t>& cols = block_cols[static_cast<std::size_t>(ia)];
//     const std::vector<float>& vals = block_vals[static_cast<std::size_t>(ia)];
//
//     checked_fwrite(cols.data(), sizeof(uint32_t), cols.size(), col_fs, "col_idx");
//     checked_fwrite(vals.data(), sizeof(float), vals.size(), val_fs, "values");
//
//     cumulative_nnz += static_cast<uint64_t>(cols.size());
//     checked_fwrite(&cumulative_nnz, sizeof(uint64_t), 1, row_fs, "row_ptr");
//    }
//   }
//
//   res.nnz = cumulative_nnz;
//
//   std::fclose(row_fs);
//   std::fclose(col_fs);
//   std::fclose(val_fs);
//   row_fs = nullptr;
//   col_fs = nullptr;
//   val_fs = nullptr;
//
//   std::ofstream meta(res.meta_file.c_str());
//   if (!meta.is_open()) throw std::runtime_error("Could not open metadata output file.");
//
//   meta << "format=sparse_ld_csr" << std::endl;
//   meta << "storage=streamed_upper_triangle" << std::endl;
//   meta << "n_samples_used=" << n << std::endl;
//   meta << "n_variants=" << m << std::endl;
//   meta << "nnz=" << static_cast<unsigned long long>(res.nnz) << std::endl;
//   meta << "triangle=upper" << std::endl;
//   meta << "diagonal=implicit_1" << std::endl;
//   meta << "row_ptr_file=" << res.row_file << std::endl;
//   meta << "col_idx_file=" << res.col_file << std::endl;
//   meta << "values_file=" << res.val_file << std::endl;
//   meta << "row_ptr_type=uint64" << std::endl;
//   meta << "col_idx_type=uint32" << std::endl;
//   meta << "values_type=float32" << std::endl;
//   meta << "index_base=0" << std::endl;
//   meta << "value=r" << std::endl;
//   meta << "r2_threshold=" << static_cast<double>(r2_threshold) << std::endl;
//   meta << "max_distance_bp=" << max_distance_bp << std::endl;
//   meta << "max_distance_variants=" << max_distance_variants << std::endl;
//   meta << "block_size=" << block_size << std::endl;
//   meta.close();
//
//  } catch (...) {
//   if (row_fs) std::fclose(row_fs);
//   if (col_fs) std::fclose(col_fs);
//   if (val_fs) std::fclose(val_fs);
//   throw;
//  }
//
//  return res;
// }
//
// // -----------------------------------------------------------------------------
// // Rcpp helpers
// // -----------------------------------------------------------------------------
//
// static Rcpp::NumericVector u64_to_numeric(const std::vector<uint64_t>& x) {
//  Rcpp::NumericVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<double>(x[i]);
//  return out;
// }
//
// static Rcpp::IntegerVector u32_to_integer_1based(const std::vector<uint32_t>& x) {
//  Rcpp::IntegerVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<int>(x[i]) + 1;
//  return out;
// }
//
// static Rcpp::NumericVector f32_to_numeric(const std::vector<float>& x) {
//  Rcpp::NumericVector out(x.size());
//  for (std::size_t i = 0; i < x.size(); ++i) out[i] = static_cast<double>(x[i]);
//  return out;
// }
//
// static std::vector<std::string> copy_bed_files(Rcpp::CharacterVector bed_files) {
//  std::vector<std::string> out(bed_files.size());
//  for (int i = 0; i < bed_files.size(); ++i) out[i] = Rcpp::as<std::string>(bed_files[i]);
//  return out;
// }
//
// static std::vector<std::vector<int>> copy_int_list(Rcpp::List xlist) {
//  std::vector<std::vector<int>> out(xlist.size());
//  for (int f = 0; f < xlist.size(); ++f) {
//   Rcpp::IntegerVector x = xlist[f];
//   out[static_cast<std::size_t>(f)].resize(x.size());
//   for (int i = 0; i < x.size(); ++i) out[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)] = x[i];
//  }
//  return out;
// }
//
// static std::vector<double> flatten_af_list(Rcpp::Nullable<Rcpp::List> af) {
//  std::vector<double> out;
//  if (af.isNotNull()) {
//   Rcpp::List af_list = Rcpp::as<Rcpp::List>(af.get());
//   for (int f = 0; f < af_list.size(); ++f) {
//    Rcpp::NumericVector x = af_list[f];
//    for (int i = 0; i < x.size(); ++i) out.push_back(x[i]);
//   }
//  }
//  return out;
// }
//
// static std::vector<int> flatten_pos_list_or_empty(Rcpp::Nullable<Rcpp::List> pos_bp) {
//  std::vector<int> out;
//  if (pos_bp.isNotNull()) {
//   Rcpp::List pos_list = Rcpp::as<Rcpp::List>(pos_bp.get());
//   for (int f = 0; f < pos_list.size(); ++f) {
//    Rcpp::IntegerVector x = pos_list[f];
//    for (int i = 0; i < x.size(); ++i) out.push_back(x[i]);
//   }
//  }
//  return out;
// }
//
// static std::vector<int> copy_rows0_or_empty(Rcpp::Nullable<Rcpp::IntegerVector> rows, int n) {
//  std::vector<int> out;
//  if (rows.isNotNull()) {
//   Rcpp::IntegerVector r(rows);
//   out.resize(r.size());
//   for (int i = 0; i < r.size(); ++i) out[static_cast<std::size_t>(i)] = r[i] - 1;
//
//   if (static_cast<int>(out.size()) == n) {
//    bool identity = true;
//    for (int i = 0; i < n; ++i) {
//     if (out[static_cast<std::size_t>(i)] != i) {
//      identity = false;
//      break;
//     }
//    }
//    if (identity) out.clear();
//   }
//  }
//  return out;
// }
//
//
// // -----------------------------------------------------------------------------
// // Rcpp exports
// // -----------------------------------------------------------------------------
//
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_stream_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   std::string out_prefix,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  if (out_prefix.empty()) Rcpp::stop("out_prefix must not be empty.");
//
//  std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
//  std::vector<int> rows0 = copy_rows0_or_empty(rows, n);
//
//  const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
//  const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bed_files_cpp,
//   n,
//   rows0_ptr,
//   n_rows,
//   cls_by_file
//  );
//
//  std::vector<double> af_cpp = flatten_af_list(af);
//  const bool af_computed = af_cpp.empty();
//  if (af_computed) af_cpp = compute_af_from_packed(G);
//
//  std::vector<int> pos_cpp = flatten_pos_list_or_empty(pos_bp);
//  const int* pos_ptr = pos_cpp.empty() ? nullptr : pos_cpp.data();
//
//  StreamingCSRResult res = sparse_ld_packed_stream_to_disk_core(
//   G,
//   af_cpp.data(),
//   pos_ptr,
//   max_distance_bp,
//   max_distance_variants,
//   static_cast<float>(r2_threshold),
//   block_size,
//   nthreads,
//   out_prefix
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr_file") = res.row_file,
//   Rcpp::Named("col_idx_file") = res.col_file,
//   Rcpp::Named("values_file") = res.val_file,
//   Rcpp::Named("meta_file") = res.meta_file,
//   Rcpp::Named("nrow") = res.m,
//   Rcpp::Named("ncol") = res.m,
//   Rcpp::Named("max_distance_variants") = max_distance_variants,
//   Rcpp::Named("nnz") = static_cast<double>(res.nnz),
//   Rcpp::Named("n_bed") = n,
//   Rcpp::Named("n_used") = G.n,
//   Rcpp::Named("af_computed") = af_computed,
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = 0,
//   Rcpp::Named("value") = "r"
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_to_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  std::vector<std::string> bed_files_cpp = copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = copy_int_list(cls);
//  std::vector<int> rows0 = copy_rows0_or_empty(rows, n);
//
//  const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
//  const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());
//
//  PackedBedMatrix G = read_bedfiles_to_packed_matrix(
//   bed_files_cpp,
//   n,
//   rows0_ptr,
//   n_rows,
//   cls_by_file
//  );
//
//  std::vector<double> af_cpp = flatten_af_list(af);
//  const bool af_computed = af_cpp.empty();
//  if (af_computed) af_cpp = compute_af_from_packed(G);
//
//  std::vector<int> pos_cpp = flatten_pos_list_or_empty(pos_bp);
//  const int* pos_ptr = pos_cpp.empty() ? nullptr : pos_cpp.data();
//
//  SparseLDCSR csr = sparse_ld_packed_core(
//   G,
//   af_cpp.data(),
//   pos_ptr,
//   max_distance_bp,
//   max_distance_variants,
//   static_cast<float>(r2_threshold),
//   block_size,
//   nthreads
//  );
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr") = u64_to_numeric(csr.row_ptr),
//   Rcpp::Named("col_idx") = u32_to_integer_1based(csr.col_idx),
//   Rcpp::Named("values") = f32_to_numeric(csr.values),
//   Rcpp::Named("nrow") = G.m,
//   Rcpp::Named("ncol") = G.m,
//   Rcpp::Named("nnz") = static_cast<double>(csr.row_ptr.back()),
//   Rcpp::Named("n_bed") = n,
//   Rcpp::Named("n_used") = G.n,
//   Rcpp::Named("af_computed") = af_computed,
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = 1,
//   Rcpp::Named("value") = "r"
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_write_CSR(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   std::string out_prefix,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> af = R_NilValue,
//   Rcpp::Nullable<Rcpp::List> pos_bp = R_NilValue,
//   int max_distance_bp = 1000000,
//   int max_distance_variants = 1000,
//   double r2_threshold = 0.01,
//   int block_size = 1024,
//   int nthreads = 1
// ) {
//  return sparseLD_stream_CSR(
//   bed_files,
//   n,
//   cls,
//   out_prefix,
//   rows,
//   af,
//   pos_bp,
//   max_distance_bp,
//   max_distance_variants,
//   r2_threshold,
//   block_size,
//   nthreads
//  );
// }
//
// // [[Rcpp::export]]
// Rcpp::List sparseLD_read_CSR(
//   std::string prefix,
//   bool one_based = true
// ) {
//  const std::string meta_file = prefix + ".meta.txt";
//  const std::string row_file = prefix + ".row_ptr.u64.bin";
//  const std::string col_file = prefix + ".col_idx.u32.0based.bin";
//  const std::string val_file = prefix + ".values.f32.bin";
//
//  std::ifstream meta(meta_file.c_str());
//  if (!meta.is_open()) Rcpp::stop("Could not open metadata file: %s", meta_file);
//
//  int m = -1;
//  double nnz_d = -1.0;
//
//  int max_distance_bp = -1;
//  int max_distance_variants = -1;
//  double r2_threshold = R_NaN;
//  int block_size = -1;
//  int nthreads = -1;
//  int n_bed = -1;
//  int n_used = -1;
//
//  std::string line;
//  while (std::getline(meta, line)) {
//   const std::string key_m = "n_variants=";
//   const std::string key_nnz = "nnz=";
//   const std::string key_max_bp = "max_distance_bp=";
//   const std::string key_max_var = "max_distance_variants=";
//   const std::string key_r2 = "r2_threshold=";
//   const std::string key_block = "block_size=";
//   const std::string key_threads = "nthreads=";
//   const std::string key_n_bed = "n_bed=";
//   const std::string key_n_used = "n_used=";
//
//   if (line.rfind(key_m, 0) == 0) {
//    m = std::stoi(line.substr(key_m.size()));
//   } else if (line.rfind(key_nnz, 0) == 0) {
//    nnz_d = std::stod(line.substr(key_nnz.size()));
//   } else if (line.rfind(key_max_bp, 0) == 0) {
//    max_distance_bp = std::stoi(line.substr(key_max_bp.size()));
//   } else if (line.rfind(key_max_var, 0) == 0) {
//    max_distance_variants = std::stoi(line.substr(key_max_var.size()));
//   } else if (line.rfind(key_r2, 0) == 0) {
//    r2_threshold = std::stod(line.substr(key_r2.size()));
//   } else if (line.rfind(key_block, 0) == 0) {
//    block_size = std::stoi(line.substr(key_block.size()));
//   } else if (line.rfind(key_threads, 0) == 0) {
//    nthreads = std::stoi(line.substr(key_threads.size()));
//   } else if (line.rfind(key_n_bed, 0) == 0) {
//    n_bed = std::stoi(line.substr(key_n_bed.size()));
//   } else if (line.rfind(key_n_used, 0) == 0) {
//    n_used = std::stoi(line.substr(key_n_used.size()));
//   }
//  }
//  meta.close();
//
//  if (m <= 0) Rcpp::stop("Could not read n_variants from metadata.");
//  if (nnz_d < 0.0) Rcpp::stop("Could not read nnz from metadata.");
//
//  const std::size_t nnz = static_cast<std::size_t>(nnz_d);
//
//  std::vector<uint64_t> row_ptr_u64(static_cast<std::size_t>(m) + 1);
//  std::vector<uint32_t> col_idx_u32(nnz);
//  std::vector<float> values_f32(nnz);
//
//  auto read_file = [](const std::string& path, void* data, std::size_t nbytes) {
//   FILE* fs = std::fopen(path.c_str(), "rb");
//   if (!fs) throw std::runtime_error("Could not open file: " + path);
//   const std::size_t got = std::fread(data, 1, nbytes, fs);
//   std::fclose(fs);
//   if (got != nbytes) throw std::runtime_error("Short read from file: " + path);
//  };
//
//  read_file(row_file, row_ptr_u64.data(), row_ptr_u64.size() * sizeof(uint64_t));
//  read_file(col_file, col_idx_u32.data(), col_idx_u32.size() * sizeof(uint32_t));
//  read_file(val_file, values_f32.data(), values_f32.size() * sizeof(float));
//
//  Rcpp::NumericVector row_ptr(row_ptr_u64.size());
//  for (std::size_t i = 0; i < row_ptr_u64.size(); ++i) {
//   row_ptr[i] = static_cast<double>(row_ptr_u64[i]);
//  }
//
//  Rcpp::IntegerVector col_idx(col_idx_u32.size());
//  for (std::size_t i = 0; i < col_idx_u32.size(); ++i) {
//   col_idx[i] = static_cast<int>(col_idx_u32[i]) + (one_based ? 1 : 0);
//  }
//
//  Rcpp::NumericVector values(values_f32.size());
//  for (std::size_t i = 0; i < values_f32.size(); ++i) {
//   values[i] = static_cast<double>(values_f32[i]);
//  }
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr") = row_ptr,
//   Rcpp::Named("col_idx") = col_idx,
//   Rcpp::Named("values") = values,
//
//   Rcpp::Named("nrow") = m,
//   Rcpp::Named("ncol") = m,
//   Rcpp::Named("nnz") = static_cast<double>(nnz),
//
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = one_based ? 1 : 0,
//   Rcpp::Named("value") = "r",
//
//   Rcpp::Named("max_distance_bp") =
//    max_distance_bp >= 0 ? Rcpp::wrap(max_distance_bp) : Rcpp::wrap(NA_INTEGER),
//
//     Rcpp::Named("max_distance_variants") =
//      max_distance_variants >= 0 ? Rcpp::wrap(max_distance_variants) : Rcpp::wrap(NA_INTEGER),
//
//       Rcpp::Named("r2_threshold") =
//        std::isfinite(r2_threshold) ? Rcpp::wrap(r2_threshold) : Rcpp::wrap(NA_REAL),
//
//        Rcpp::Named("block_size") =
//         block_size >= 0 ? Rcpp::wrap(block_size) : Rcpp::wrap(NA_INTEGER),
//
//          Rcpp::Named("nthreads") =
//           nthreads >= 0 ? Rcpp::wrap(nthreads) : Rcpp::wrap(NA_INTEGER),
//
//            Rcpp::Named("n_bed") =
//             n_bed >= 0 ? Rcpp::wrap(n_bed) : Rcpp::wrap(NA_INTEGER),
//
//              Rcpp::Named("n_used") =
//               n_used >= 0 ? Rcpp::wrap(n_used) : Rcpp::wrap(NA_INTEGER)
//  );
// }



// // [[Rcpp::export]]
// Rcpp::List sparseLD_read_CSR(
//   std::string prefix,
//   bool one_based = true
// ) {
//  const std::string meta_file = prefix + ".meta.txt";
//  const std::string row_file = prefix + ".row_ptr.u64.bin";
//  const std::string col_file = prefix + ".col_idx.u32.0based.bin";
//  const std::string val_file = prefix + ".values.f32.bin";
//
//  std::ifstream meta(meta_file.c_str());
//  if (!meta.is_open()) Rcpp::stop("Could not open metadata file: %s", meta_file);
//
//  int m = -1;
//  double nnz_d = -1.0;
//  std::string line;
//  while (std::getline(meta, line)) {
//   const std::string key_m = "n_variants=";
//   const std::string key_nnz = "nnz=";
//   if (line.rfind(key_m, 0) == 0) {
//    m = std::stoi(line.substr(key_m.size()));
//   } else if (line.rfind(key_nnz, 0) == 0) {
//    nnz_d = std::stod(line.substr(key_nnz.size()));
//   }
//  }
//  meta.close();
//
//  if (m <= 0) Rcpp::stop("Could not read n_variants from metadata.");
//  if (nnz_d < 0.0) Rcpp::stop("Could not read nnz from metadata.");
//
//  const std::size_t nnz = static_cast<std::size_t>(nnz_d);
//
//  std::vector<uint64_t> row_ptr_u64(static_cast<std::size_t>(m) + 1);
//  std::vector<uint32_t> col_idx_u32(nnz);
//  std::vector<float> values_f32(nnz);
//
//  auto read_file = [](const std::string& path, void* data, std::size_t nbytes) {
//   FILE* fs = std::fopen(path.c_str(), "rb");
//   if (!fs) throw std::runtime_error("Could not open file: " + path);
//   const std::size_t got = std::fread(data, 1, nbytes, fs);
//   std::fclose(fs);
//   if (got != nbytes) throw std::runtime_error("Short read from file: " + path);
//  };
//
//  read_file(row_file, row_ptr_u64.data(), row_ptr_u64.size() * sizeof(uint64_t));
//  read_file(col_file, col_idx_u32.data(), col_idx_u32.size() * sizeof(uint32_t));
//  read_file(val_file, values_f32.data(), values_f32.size() * sizeof(float));
//
//  Rcpp::NumericVector row_ptr(row_ptr_u64.size());
//  for (std::size_t i = 0; i < row_ptr_u64.size(); ++i) {
//   row_ptr[i] = static_cast<double>(row_ptr_u64[i]);
//  }
//
//  Rcpp::IntegerVector col_idx(col_idx_u32.size());
//  for (std::size_t i = 0; i < col_idx_u32.size(); ++i) {
//   col_idx[i] = static_cast<int>(col_idx_u32[i]) + (one_based ? 1 : 0);
//  }
//
//  Rcpp::NumericVector values(values_f32.size());
//  for (std::size_t i = 0; i < values_f32.size(); ++i) {
//   values[i] = static_cast<double>(values_f32[i]);
//  }
//
//  return Rcpp::List::create(
//   Rcpp::Named("row_ptr") = row_ptr,
//   Rcpp::Named("col_idx") = col_idx,
//   Rcpp::Named("values") = values,
//   Rcpp::Named("nrow") = m,
//   Rcpp::Named("ncol") = m,
//   Rcpp::Named("nnz") = static_cast<double>(nnz),
//   Rcpp::Named("upper_triangle") = true,
//   Rcpp::Named("diag") = "implicit_1",
//   Rcpp::Named("index_base") = one_based ? 1 : 0,
//   Rcpp::Named("value") = "r"
//  );
// }
