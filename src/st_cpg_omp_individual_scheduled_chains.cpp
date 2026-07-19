// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "packed_bed.h"
#include "blr_bed_scheduled_bayesc_rng.h"
#include "blr_bed_scheduled_bayesc_types.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

// =============================================================================
// Packed BED marker-wise ST BayesC sampler with scheduled null-marker updates
// and multiple independent chains.
//
// Main difference from stblr_cpg_omp_bed_marker_scheduled():
//   - BED is read/packed once.
//   - Marker maps are built once.
//   - Independent chains are run in parallel over chain x trait jobs.
//   - Posterior summaries are averaged across chains before returning.
//   - BED setup uses a blocked reader, modeled after mtgrsbed_core(), to speed
//     up raw BED -> selected-row packed matrix construction.
//
// Return layout uses the same 23-slot base structure as the single-chain
// scheduled marker function, with optional CSR-compatible chain summaries:
//   result[[1]]  / result[0]  = posterior mean effects, averaged over chains
//   result[[2]]  / result[1]  = posterior inclusion probabilities, averaged over chains
//   result[[5]]  / result[4]  = final effects, averaged over chains
//   result[[6]]  / result[5]  = final indicators, averaged over chains
//   result[[8:10]] / result[7:9] = variance traces, averaged over chains
//   result[[17:18]] / result[16:17] = final and posterior mean pi, averaged over chains
//   result[[19]] / result[18] = diagnostics: mean total log-CPO, mean log-CPO per individual,
//                                  mean seconds, max seconds
//   result[[20]] / result[19] = diagnostics: mean nsamples, n_used
//   result[[24:29]] / result[23:28] = bm_sd, bm_min, bm_max,
//                                      dm_sd, dm_min, dm_max
// =============================================================================

using MarkerMapSTScheduledChains=sblr::core::BedScheduledBayesCMarkerMap;

struct FastPackedBedMatrix {
 int n = 0;
 int m = 0;
 std::size_t nbytes = 0;
 std::size_t stride = 0;
 std::vector<uint8_t> data;

 inline uint8_t* row(int marker) {
  return data.data() + static_cast<std::size_t>(marker) * stride;
 }

 inline const uint8_t* row(int marker) const {
  return data.data() + static_cast<std::size_t>(marker) * stride;
 }
};

static inline unsigned int fast_get_bed_code_from_row(
  const uint8_t* packed,
  int sample
) {
 return (packed[static_cast<std::size_t>(sample >> 2)] >> (2 * (sample & 3))) & 3u;
}

static FastPackedBedMatrix read_bedfiles_to_fast_packed_matrix_blocked(
  const std::vector<std::string>& bed_files,
  int n,
  const int* rows0,
  int n_rows,
  const std::vector<std::vector<int>>& cls_by_file,
  int read_block_size,
  int nthreads
) {
 if (n <= 0) {
  throw std::runtime_error("read_bedfiles_to_fast_packed_matrix_blocked: n must be positive.");
 }
 if (read_block_size <= 0) {
  throw std::runtime_error("read_bedfiles_to_fast_packed_matrix_blocked: read_block_size must be positive.");
 }
 if (bed_files.size() != cls_by_file.size()) {
  throw std::runtime_error("read_bedfiles_to_fast_packed_matrix_blocked: bed_files and cls lengths differ.");
 }

 const int n_used = n_rows > 0 ? n_rows : n;
 const std::size_t raw_nbytes = (static_cast<std::size_t>(n) + 3u) / 4u;
 const std::size_t out_nbytes = (static_cast<std::size_t>(n_used) + 3u) / 4u;

 int m_total = 0;
 for (std::size_t f = 0; f < cls_by_file.size(); ++f) {
  m_total += static_cast<int>(cls_by_file[f].size());
 }
 if (m_total <= 0) {
  throw std::runtime_error("read_bedfiles_to_fast_packed_matrix_blocked: no markers selected.");
 }

 FastPackedBedMatrix G;
 G.n = n_used;
 G.m = m_total;
 G.nbytes = out_nbytes;
 G.stride = out_nbytes;
 G.data.assign(static_cast<std::size_t>(m_total) * out_nbytes, static_cast<uint8_t>(0));

 std::vector<int> src_byte(static_cast<std::size_t>(n_used));
 std::vector<int> src_shift(static_cast<std::size_t>(n_used));
 std::vector<int> dst_byte(static_cast<std::size_t>(n_used));
 std::vector<int> dst_shift(static_cast<std::size_t>(n_used));

 if (n_rows > 0) {
  for (int k = 0; k < n_used; ++k) {
   const int r = rows0[k];
   if (r < 0 || r >= n) {
    throw std::runtime_error("read_bedfiles_to_fast_packed_matrix_blocked: rows0 contains index outside [0, n).");
   }
   src_byte[static_cast<std::size_t>(k)] = r >> 2;
   src_shift[static_cast<std::size_t>(k)] = 2 * (r & 3);
   dst_byte[static_cast<std::size_t>(k)] = k >> 2;
   dst_shift[static_cast<std::size_t>(k)] = 2 * (k & 3);
  }
 } else {
  for (int k = 0; k < n_used; ++k) {
   src_byte[static_cast<std::size_t>(k)] = k >> 2;
   src_shift[static_cast<std::size_t>(k)] = 2 * (k & 3);
   dst_byte[static_cast<std::size_t>(k)] = k >> 2;
   dst_shift[static_cast<std::size_t>(k)] = 2 * (k & 3);
  }
 }

 int global_marker = 0;

 for (std::size_t f = 0; f < bed_files.size(); ++f) {
  FILE* fs = std::fopen(bed_files[f].c_str(), "rb");
  if (!fs) {
   throw std::runtime_error("Could not open BED file: " + bed_files[f]);
  }

  unsigned char header[3];
  const std::size_t nhead = std::fread(header, sizeof(unsigned char), 3, fs);
  if (nhead != 3 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
   std::fclose(fs);
   throw std::runtime_error("Invalid or unsupported PLINK BED header in file: " + bed_files[f]);
  }

  const std::vector<int>& cls_f = cls_by_file[f];
  const int mf = static_cast<int>(cls_f.size());
  std::vector<uint8_t> block_buffer(
    static_cast<std::size_t>(read_block_size) * raw_nbytes
  );

  for (int i0 = 0; i0 < mf; i0 += read_block_size) {
   const int imax = std::min(i0 + read_block_size, mf);
   const int mlen = imax - i0;

   for (int ii = 0; ii < mlen; ++ii) {
    const int cls_index_1based = cls_f[static_cast<std::size_t>(i0 + ii)];
    if (cls_index_1based <= 0) {
     std::fclose(fs);
     throw std::runtime_error("read_bedfiles_to_fast_packed_matrix_blocked: cls contains non-positive marker index.");
    }

    const long long offset =
     3LL + static_cast<long long>(cls_index_1based - 1) * static_cast<long long>(raw_nbytes);

    if (std::fseek(fs, offset, SEEK_SET) != 0) {
     std::fclose(fs);
     throw std::runtime_error("fseek failed while reading BED file: " + bed_files[f]);
    }

    uint8_t* raw = block_buffer.data() + static_cast<std::size_t>(ii) * raw_nbytes;
    const std::size_t got = std::fread(raw, sizeof(uint8_t), raw_nbytes, fs);
    if (got != raw_nbytes) {
     std::fclose(fs);
     throw std::runtime_error("Short read while reading BED file: " + bed_files[f]);
    }
   }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
   for (int ii = 0; ii < mlen; ++ii) {
    const uint8_t* raw = block_buffer.data() + static_cast<std::size_t>(ii) * raw_nbytes;
    uint8_t* dst = G.row(global_marker + ii);
    std::memset(dst, 0, out_nbytes);

    if (n_rows == 0 && n_used == n && out_nbytes == raw_nbytes) {
     std::memcpy(dst, raw, raw_nbytes);
    } else {
     for (int k = 0; k < n_used; ++k) {
      const std::size_t ku = static_cast<std::size_t>(k);
      const uint8_t code =
       static_cast<uint8_t>((raw[static_cast<std::size_t>(src_byte[ku])] >> src_shift[ku]) & 3u);
      dst[static_cast<std::size_t>(dst_byte[ku])] |=
       static_cast<uint8_t>(code << dst_shift[ku]);
     }
    }
   }

   global_marker += mlen;
  }

  std::fclose(fs);
 }

 if (global_marker != m_total) {
  throw std::runtime_error("read_bedfiles_to_fast_packed_matrix_blocked: internal marker count mismatch.");
 }

 return G;
}

static std::vector<double> compute_af_from_fast_packed(
  const FastPackedBedMatrix& G
) {
 std::vector<double> af(static_cast<std::size_t>(G.m), 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
 for (int j = 0; j < G.m; ++j) {
  const uint8_t* row = G.row(j);
  double dosage_sum = 0.0;
  double allele_count = 0.0;

  for (int i = 0; i < G.n; ++i) {
   const unsigned int code = fast_get_bed_code_from_row(row, i);
   if (code == 1u) {
    continue;
   }

   if (code == 0u) dosage_sum += 2.0;
   else if (code == 2u) dosage_sum += 1.0;
   else dosage_sum += 0.0;

   allele_count += 2.0;
  }

  af[static_cast<std::size_t>(j)] =
   allele_count > 0.0 ? dosage_sum / allele_count : 0.0;
 }

 return af;
}

static std::vector<std::string> copy_bed_files_scheduled_chains(Rcpp::CharacterVector bed_files) {
 std::vector<std::string> out(static_cast<std::size_t>(bed_files.size()));
 for (int i = 0; i < bed_files.size(); ++i) {
  out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(bed_files[i]);
 }
 return out;
}

static std::vector<std::vector<int>> copy_int_list_scheduled_chains(Rcpp::List xlist) {
 std::vector<std::vector<int>> out(static_cast<std::size_t>(xlist.size()));
 for (int f = 0; f < xlist.size(); ++f) {
  Rcpp::IntegerVector x = xlist[f];
  out[static_cast<std::size_t>(f)].resize(static_cast<std::size_t>(x.size()));
  for (int i = 0; i < x.size(); ++i) {
   out[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)] = x[i];
  }
 }
 return out;
}

static std::vector<int> copy_rows0_or_empty_scheduled_chains(
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  int n
) {
 std::vector<int> out;
 if (rows.isNotNull()) {
  Rcpp::IntegerVector r(rows);
  out.resize(static_cast<std::size_t>(r.size()));

  for (int i = 0; i < r.size(); ++i) {
   if (r[i] == NA_INTEGER) throw std::runtime_error("rows contains NA.");
   if (r[i] < 1 || r[i] > n) throw std::runtime_error("rows contains index outside [1, n].");
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

static std::vector<double> flatten_af_list_or_empty_scheduled_chains(Rcpp::Nullable<Rcpp::List> af) {
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

static inline void bed_marker_map_values_scheduled_chains(
  double p,
  bool scale,
  MarkerMapSTScheduledChains& map
) {
 if (scale) {
  const double denom = std::sqrt(2.0 * p * (1.0 - p));
  if (denom <= 0.0 || !std::isfinite(denom)) {
   map.val[0] = map.val[1] = map.val[2] = map.val[3] = 0.0;
  } else {
   // PLINK 2-bit BED coding used by packed_bed.h:
   // 0 -> genotype 2, 1 -> missing, 2 -> genotype 1, 3 -> genotype 0
   map.val[0] = (2.0 - 2.0 * p) / denom;
   map.val[1] = 0.0;
   map.val[2] = (1.0 - 2.0 * p) / denom;
   map.val[3] = (0.0 - 2.0 * p) / denom;
  }
 } else {
  map.val[0] = 2.0;
  map.val[1] = 2.0 * p;
  map.val[2] = 1.0;
  map.val[3] = 0.0;
 }
}

static inline double marker_xx_from_packed_scheduled_chains(
  const FastPackedBedMatrix& G,
  int marker,
  const MarkerMapSTScheduledChains& map
) {
 const uint8_t* packed = G.row(marker);
 const int n = G.n;
 const std::size_t nbytes = G.nbytes;
 double xx = 0.0;

 for (std::size_t kb = 0; kb < nbytes; ++kb) {
  const unsigned char byte = packed[kb];
  const int jbase = static_cast<int>(kb << 2);

  if (jbase + 0 < n) {
   const double x = map.val[(byte >> 0) & 3u];
   xx += x * x;
  }
  if (jbase + 1 < n) {
   const double x = map.val[(byte >> 2) & 3u];
   xx += x * x;
  }
  if (jbase + 2 < n) {
   const double x = map.val[(byte >> 4) & 3u];
   xx += x * x;
  }
  if (jbase + 3 < n) {
   const double x = map.val[(byte >> 6) & 3u];
   xx += x * x;
  }
 }

 return xx;
}

static std::vector<MarkerMapSTScheduledChains> build_marker_maps_scheduled_chains(
  const FastPackedBedMatrix& G,
  const std::vector<double>& af,
  bool scale,
  int nthreads
) {
 const int m = G.m;
 if (static_cast<int>(af.size()) != m) {
  throw std::runtime_error("build_marker_maps_scheduled_chains: af length must equal number of markers.");
 }

 std::vector<MarkerMapSTScheduledChains> maps(static_cast<std::size_t>(m));

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
 for (int j = 0; j < m; ++j) {
  MarkerMapSTScheduledChains map;
  bed_marker_map_values_scheduled_chains(af[static_cast<std::size_t>(j)], scale, map);
  map.xx = marker_xx_from_packed_scheduled_chains(G, j, map);
  maps[static_cast<std::size_t>(j)] = map;
 }

 for (int j = 0; j < m; ++j) {
  if (!std::isfinite(maps[static_cast<std::size_t>(j)].xx) ||
      maps[static_cast<std::size_t>(j)].xx <= 0.0) {
   throw std::runtime_error(
     "build_marker_maps_scheduled_chains: invalid x'x for marker " + std::to_string(j) +
      ". Check allele frequencies and monomorphic markers."
   );
  }
 }

 return maps;
}

static inline double bed_marker_dot_residual_scheduled_chains(
  const FastPackedBedMatrix& G,
  int marker,
  const MarkerMapSTScheduledChains& map,
  const double* e
) {
 const uint8_t* packed = G.row(marker);
 const int n = G.n;
 const std::size_t nbytes = G.nbytes;
 double out = 0.0;

 for (std::size_t kb = 0; kb < nbytes; ++kb) {
  const unsigned char byte = packed[kb];
  const int jbase = static_cast<int>(kb << 2);

  if (jbase + 0 < n) out += map.val[(byte >> 0) & 3u] * e[jbase + 0];
  if (jbase + 1 < n) out += map.val[(byte >> 2) & 3u] * e[jbase + 1];
  if (jbase + 2 < n) out += map.val[(byte >> 4) & 3u] * e[jbase + 2];
  if (jbase + 3 < n) out += map.val[(byte >> 6) & 3u] * e[jbase + 3];
 }

 return out;
}

static inline void bed_marker_update_residual_scheduled_chains(
  const FastPackedBedMatrix& G,
  int marker,
  const MarkerMapSTScheduledChains& map,
  double* e,
  double diff
) {
 const uint8_t* packed = G.row(marker);
 const int n = G.n;
 const std::size_t nbytes = G.nbytes;

 for (std::size_t kb = 0; kb < nbytes; ++kb) {
  const unsigned char byte = packed[kb];
  const int jbase = static_cast<int>(kb << 2);

  if (jbase + 0 < n) e[jbase + 0] -= map.val[(byte >> 0) & 3u] * diff;
  if (jbase + 1 < n) e[jbase + 1] -= map.val[(byte >> 2) & 3u] * diff;
  if (jbase + 2 < n) e[jbase + 2] -= map.val[(byte >> 4) & 3u] * diff;
  if (jbase + 3 < n) e[jbase + 3] -= map.val[(byte >> 6) & 3u] * diff;
 }
}

static std::vector<int> make_marker_order_from_sets_scheduled_chains(
  const std::vector<int>& sets,
  int m
) {
 if (static_cast<int>(sets.size()) != m) {
  throw std::runtime_error("make_marker_order_from_sets_scheduled_chains: sets length must equal m.");
 }

 std::vector<int> labels;
 labels.reserve(static_cast<std::size_t>(m));
 std::unordered_map<int, int> label_to_block;

 for (int j = 0; j < m; ++j) {
  const int lab = sets[static_cast<std::size_t>(j)];
  if (label_to_block.find(lab) == label_to_block.end()) {
   const int block_id = static_cast<int>(labels.size());
   label_to_block[lab] = block_id;
   labels.push_back(lab);
  }
 }

 std::vector<std::vector<int>> block_markers(labels.size());
 for (int j = 0; j < m; ++j) {
  const int block_id = label_to_block[sets[static_cast<std::size_t>(j)]];
  block_markers[static_cast<std::size_t>(block_id)].push_back(j);
 }

 std::vector<int> order;
 order.reserve(static_cast<std::size_t>(m));
 for (std::size_t b = 0; b < block_markers.size(); ++b) {
  for (int j : block_markers[b]) order.push_back(j);
 }
 return order;
}

static inline double sample_marker_scheduled_chains(
  const FastPackedBedMatrix& G,
  int marker,
  const MarkerMapSTScheduledChains& map,
  const std::vector<double>& pi,
  double vb,
  double vei,
  arma::vec& e,
  double& b_j,
  int& d_j,
  sblr::core::BedScheduledBayesCChainRng& rng
) {
 const double xj2 = map.xx;
 const double vei_safe = std::max(vei, 1e-300);
 const double pi0 = std::max(pi[0], 1e-300);
 const double pi1 = std::max(pi[1], 1e-300);

 const double xte = bed_marker_dot_residual_scheduled_chains(G, marker, map, e.memptr());
 const double score = xte + xj2 * b_j;
 const double denom = std::max(vei_safe + xj2 * vb, 1e-300);

 const double logBF =
  0.5 * std::log(vei_safe / denom) +
  0.5 * score * score * vb / (vei_safe * denom);

 const double logp1 = std::log(pi1) + logBF;
 const double logp0 = std::log(pi0);
 const double delta_log = logp0 - logp1;

 double p1 = 0.0;
 if (delta_log > 35.0) p1 = 0.0;
 else if (delta_log < -35.0) p1 = 1.0;
 else p1 = 1.0 / (1.0 + std::exp(delta_log));

 const int d_new = (rng.uniform(rng.engine) < p1) ? 1 : 0;
 double b_new = 0.0;

 if (d_new == 1) {
  const double lhs = xj2 + vei_safe / vb;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * rng.normal(rng.engine);
 }

 const double diff = b_new - b_j;
 if (diff != 0.0) {
  bed_marker_update_residual_scheduled_chains(G, marker, map, e.memptr(), diff);
 }

 b_j = b_new;
 d_j = d_new;
 return p1;
}

static inline void sampleB_sparse_scheduled_chains(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (d(iu) > 0) {
   ssb += b(iu) * b(iu);
   dfb += 1.0;
  }
 }

 const double scale = ssb + nub * ssb_prior;
 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 vb = std::max(scale / chi2, 1e-12);
}

static inline void sampleE_sparse_scheduled_chains(
  double nue,
  double& ve,
  const arma::vec& e,
  double sse_prior,
  std::mt19937& gen
) {
 const double sse = arma::dot(e, e);
 const double scale = sse + nue * sse_prior;

 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleE_sparse_scheduled_chains: invalid residual scale.");
 }

 std::chi_squared_distribution<double> rchisq(static_cast<double>(e.n_elem) + nue);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 const double ve_new = scale / chi2;

 if (!std::isfinite(ve_new) || ve_new <= 0.0) {
  throw std::runtime_error("sampleE_sparse_scheduled_chains: sampled ve is invalid.");
 }

 ve = std::max(ve_new, 1e-12);
}

static inline void update_log_inv_cpo_gaussian_scheduled_chains(
  const arma::vec& e,
  double ve,
  arma::vec& log_inv_cpo
) {
 const double ve_safe = std::max(ve, 1e-300);
 const double log_norm = -0.5 * std::log(2.0 * std::acos(-1.0) * ve_safe);
 const arma::uword n = e.n_elem;

 if (log_inv_cpo.n_elem != n) {
  throw std::runtime_error("update_log_inv_cpo_gaussian_scheduled_chains: dimension mismatch.");
 }

 for (arma::uword i = 0; i < n; ++i) {
  const double ei = e(i);
  const double loglik = log_norm - 0.5 * ei * ei / ve_safe;
  const double z = -loglik;

  if (!std::isfinite(log_inv_cpo(i))) {
   log_inv_cpo(i) = z;
  } else {
   const double a = log_inv_cpo(i);
   const double mx = std::max(a, z);
   log_inv_cpo(i) = mx + std::log(std::exp(a - mx) + std::exp(z - mx));
  }
 }
}

static inline double compute_total_log_cpo_scheduled_chains(
  const arma::vec& log_inv_cpo,
  double nsamples
) {
 if (nsamples <= 0.0) return NA_REAL;

 const double log_nsamples = std::log(nsamples);
 double total = 0.0;

 for (arma::uword i = 0; i < log_inv_cpo.n_elem; ++i) {
  if (std::isfinite(log_inv_cpo(i))) {
   total += log_nsamples - log_inv_cpo(i);
  }
 }

 return total;
}

static inline double computeG_sparse_scheduled_chains(
  const arma::vec& y,
  const arma::vec& e
) {
 double ss = 0.0;
 const int n = static_cast<int>(y.n_elem);
 for (int i = 0; i < n; ++i) {
  const double g = y(static_cast<arma::uword>(i)) - e(static_cast<arma::uword>(i));
  ss += g * g;
 }
 return ss / static_cast<double>(n);
}

static inline double computeLE_sparse_scheduled_chains(
  const arma::rowvec& b,
  const std::vector<MarkerMapSTScheduledChains>& maps,
  int n
) {
 // Linkage-equilibrium contribution under the empirical genotype scaling used
 // by the packed matrix: sum_j Var(x_j) b_j^2, with Var(x_j) = x_j'x_j / n.
 // For scale=TRUE and little missingness this is close to sum_j b_j^2.
 double vle = 0.0;
 const int m = static_cast<int>(b.n_elem);
 const double inv_n = 1.0 / static_cast<double>(n);

 for (int j = 0; j < m; ++j) {
  const double bj = b(static_cast<arma::uword>(j));
  if (bj != 0.0) {
   vle += maps[static_cast<std::size_t>(j)].xx * inv_n * bj * bj;
  }
 }

 return vle;
}

static inline void samplePi_sparse_scheduled_chains(
  const arma::Row<int>& d,
  std::vector<double>& pi,
  double pi_prior_a,
  double pi_prior_b,
  std::mt19937& gen
) {
 // pi[1] is the inclusion probability and pi[0] is the exclusion probability.
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

static arma::vec bed_xb_from_b_scheduled_chains(
  const FastPackedBedMatrix& G,
  const std::vector<MarkerMapSTScheduledChains>& maps,
  const std::vector<int>& marker_order,
  const arma::rowvec& b
) {
 arma::vec xb(G.n, arma::fill::zeros);

 for (int marker : marker_order) {
  const double bj = b(static_cast<arma::uword>(marker));
  if (bj != 0.0) {
   // xb += x_j * bj, using residual-update routine with negative diff.
   bed_marker_update_residual_scheduled_chains(
    G,
    marker,
    maps[static_cast<std::size_t>(marker)],
        xb.memptr(),
        -bj
   );
  }
 }
 return xb;
}

static inline int adaptive_skip_length_scheduled_chains(
  double p1,
  int null_skip_base,
  int null_skip_max
) {
 if (null_skip_base <= 1) return 1;

 int skip = null_skip_base;

 if (p1 < 1e-6) skip = 4 * null_skip_base;
 else if (p1 < 1e-5) skip = 2 * null_skip_base;
 else if (p1 < 1e-4) skip = null_skip_base;
 else if (p1 < 1e-3) skip = std::max(1, null_skip_base / 2);
 else skip = 1;

 if (null_skip_max > 0) skip = std::min(skip, null_skip_max);
 return std::max(1, skip);
}

#include "blr_bed_scheduled_bayesc_core_impl.h"
#include "blr_bed_scheduled_bayesc_aggregate_impl.h"

struct BedScheduledBayesCBindingMetadata {
 int nit;
 int nburn;
 int nthin;
 int n_used;
};

static Rcpp::List stblr_bed_scheduled_bayesc_result_to_raw(
 const sblr::core::BedScheduledBayesCExecutionResult& result,
 const BedScheduledBayesCBindingMetadata& metadata
) {
 const int m=result.marker_count;
 const int nt=result.trait_count;
 const int n_trace=result.trace_length;
 auto marker_matrix=[&](const arma::mat& x) {
  Rcpp::NumericMatrix out(m,nt);
  for (int t=0; t<nt; ++t) for (int i=0; i<m; ++i)
   out(i,t)=x(static_cast<arma::uword>(t),static_cast<arma::uword>(i));
  return out;
 };
 auto trace_matrix=[&](const arma::mat& x) {
  Rcpp::NumericMatrix out(n_trace,nt);
  for (int t=0; t<nt; ++t) for (int it=0; it<n_trace; ++it)
   out(it,t)=x(static_cast<arma::uword>(t),static_cast<arma::uword>(it));
  return out;
 };
 auto diagonal_matrix=[&](const arma::vec& x) {
  Rcpp::NumericMatrix out(nt,nt);
  for (int t=0; t<nt; ++t) out(t,t)=x(static_cast<arma::uword>(t));
  return out;
 };

 Rcpp::NumericMatrix pi_final(nt,2),pi_mean(nt,2);
 Rcpp::NumericVector nsamples_out(nt),log_cpo_out(nt),mean_log_cpo_out(nt),
  seconds_mean_out(nt),seconds_max_out(nt);
 Rcpp::IntegerVector n_used_out(nt);
 for (int t=0; t<nt; ++t) {
  const arma::uword tu=static_cast<arma::uword>(t);
  pi_final(t,0)=1.0-result.final_inclusion_probability(tu);
  pi_final(t,1)=result.final_inclusion_probability(tu);
  pi_mean(t,0)=1.0-result.mean_inclusion_probability(tu);
  pi_mean(t,1)=result.mean_inclusion_probability(tu);
  nsamples_out[t]=result.mean_retained_samples(tu);
  n_used_out[t]=metadata.n_used;
  log_cpo_out[t]=result.mean_total_log_cpo(tu);
  mean_log_cpo_out[t]=result.mean_log_cpo(tu);
  seconds_mean_out[t]=result.mean_seconds(tu);
  seconds_max_out[t]=result.max_seconds(tu);
 }

 Rcpp::List marker=Rcpp::List::create(
  Rcpp::Named("bm")=marker_matrix(result.bm),
  Rcpp::Named("dm")=marker_matrix(result.dm),
  Rcpp::Named("wy")=marker_matrix(result.wy),
  Rcpp::Named("r")=marker_matrix(result.residual_scores),
  Rcpp::Named("b")=marker_matrix(result.final_effects),
  Rcpp::Named("state")=marker_matrix(result.final_states),
  Rcpp::Named("bm_sd")=marker_matrix(result.bm_sd),
  Rcpp::Named("bm_min")=marker_matrix(result.bm_min),
  Rcpp::Named("bm_max")=marker_matrix(result.bm_max),
  Rcpp::Named("dm_sd")=marker_matrix(result.dm_sd),
  Rcpp::Named("dm_min")=marker_matrix(result.dm_min),
  Rcpp::Named("dm_max")=marker_matrix(result.dm_max));
 Rcpp::List trace=Rcpp::List::create(
  Rcpp::Named("vbs")=trace_matrix(result.marker_variance_trace),
  Rcpp::Named("vgs")=trace_matrix(result.genetic_variance_trace),
  Rcpp::Named("ves")=trace_matrix(result.residual_variance_trace),
  Rcpp::Named("vle")=trace_matrix(result.vle_trace),
  Rcpp::Named("vld")=trace_matrix(result.vld_trace),
  Rcpp::Named("pis")=trace_matrix(result.inclusion_trace));
 Rcpp::List variance=Rcpp::List::create(
  Rcpp::Named("covb")=diagonal_matrix(result.final_marker_variance),
  Rcpp::Named("covg")=diagonal_matrix(result.final_genetic_variance),
  Rcpp::Named("cove")=diagonal_matrix(result.final_residual_variance),
  Rcpp::Named("vb")=diagonal_matrix(result.final_marker_variance),
  Rcpp::Named("vg")=diagonal_matrix(result.final_genetic_variance),
  Rcpp::Named("ve")=diagonal_matrix(result.final_residual_variance));
 Rcpp::List diagnostics=Rcpp::List::create(
  Rcpp::Named("nsamples")=nsamples_out,Rcpp::Named("n_used")=n_used_out,
  Rcpp::Named("log_cpo")=log_cpo_out,
  Rcpp::Named("mean_log_cpo")=mean_log_cpo_out,
  Rcpp::Named("seconds_mean")=seconds_mean_out,
  Rcpp::Named("seconds_max")=seconds_max_out,Rcpp::Named("ld_swap")=R_NilValue);
 Rcpp::List selection=Rcpp::List::create(
  Rcpp::Named("enabled")=false,Rcpp::Named("fixed")=false,
  Rcpp::Named("scale")="standardized_genotype_effect",
  Rcpp::Named("trace")=R_NilValue,Rcpp::Named("mean")=R_NilValue,
  Rcpp::Named("sd")=R_NilValue,Rcpp::Named("min")=R_NilValue,
  Rcpp::Named("max")=R_NilValue,Rcpp::Named("acceptance")=R_NilValue);
 Rcpp::List raw=Rcpp::List::create(
  Rcpp::Named("schema")=Rcpp::List::create(
   Rcpp::Named("class")="stblr_raw",Rcpp::Named("version")=1),
  Rcpp::Named("meta")=Rcpp::List::create(
   Rcpp::Named("model")="bayesc",
   Rcpp::Named("backend")="bed_scheduled_chains_bayesc",
   Rcpp::Named("data_level")="individual",Rcpp::Named("prior_type")="global",
   Rcpp::Named("m")=m,Rcpp::Named("nt")=nt,Rcpp::Named("n_trace")=n_trace,
   Rcpp::Named("nit")=metadata.nit,Rcpp::Named("nburn")=metadata.nburn,
   Rcpp::Named("nthin")=metadata.nthin,Rcpp::Named("nchains")=result.chain_count,
   Rcpp::Named("keep_chains")=false,Rcpp::Named("n_components")=2,
   Rcpp::Named("n_annotations")=0,Rcpp::Named("n_groups")=0),
  Rcpp::Named("marker")=marker,Rcpp::Named("trace")=trace,
  Rcpp::Named("variance")=variance,
  Rcpp::Named("pi")=Rcpp::List::create(
   Rcpp::Named("final")=pi_final,Rcpp::Named("mean")=pi_mean,
   Rcpp::Named("names")=Rcpp::CharacterVector::create("pi0","pi1")),
  Rcpp::Named("diagnostics")=diagnostics,Rcpp::Named("chains")=R_NilValue,
  Rcpp::Named("prior")=Rcpp::List::create(),Rcpp::Named("group")=Rcpp::List::create(),
  Rcpp::Named("annotation")=Rcpp::List::create(),
  Rcpp::Named("component")=Rcpp::List::create(),Rcpp::Named("selection")=selection);
 raw.attr("class")=Rcpp::CharacterVector::create("stblr_raw_v1","stblr_raw","list");
 return raw;
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_bed_marker_scheduled_chains(
  Rcpp::CharacterVector bed_files,
  int n,
  Rcpp::List cls,
  Rcpp::NumericMatrix y,
  std::vector<std::vector<double>> b_init,
  std::vector<int> sets,
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::Nullable<Rcpp::List> af,
  bool scale,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<double> pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  double adjE,
  int nit,
  int nburn,
  int nthin,
  int rebuild_every,
  int full_sweep_every,
  int null_skip_base,
  int null_skip_max,
  double candidate_threshold,
  int candidate_lifetime,
  bool skip_nulls_burnin_only,
  bool return_wy,
  bool return_r,
  int read_block_size,
  int progress_every,
  double pi_prior_a,
  double pi_prior_b,
  int nchains,
  int ncores,
  int seed
) {
 if (nit <= 0 || nburn < 0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_scheduled_chains: nit must be positive and nburn non-negative.");
 }
 if (nthin <= 0) throw std::runtime_error("nthin must be positive.");
 if (rebuild_every < 0) throw std::runtime_error("rebuild_every must be >= 0.");
 if (full_sweep_every < 0) throw std::runtime_error("full_sweep_every must be >= 0.");
 if (null_skip_base <= 0) throw std::runtime_error("null_skip_base must be positive.");
 if (null_skip_max < 0) throw std::runtime_error("null_skip_max must be >= 0.");
 if (nchains <= 0) throw std::runtime_error("nchains must be positive.");
 if (read_block_size <= 0) throw std::runtime_error("read_block_size must be positive.");
 if (progress_every < 0) throw std::runtime_error("progress_every must be >= 0.");
 if (!std::isfinite(pi_prior_a) || pi_prior_a <= 0.0) {
  throw std::runtime_error("pi_prior_a must be finite and positive.");
 }
 if (!std::isfinite(pi_prior_b) || pi_prior_b <= 0.0) {
  throw std::runtime_error("pi_prior_b must be finite and positive.");
 }
 if (!std::isfinite(candidate_threshold) || candidate_threshold < 0.0 || candidate_threshold > 1.0) {
  throw std::runtime_error("candidate_threshold must be in [0,1].");
 }
 if (candidate_lifetime < 0) throw std::runtime_error("candidate_lifetime must be >= 0.");

 std::vector<std::string> bed_files_cpp = copy_bed_files_scheduled_chains(bed_files);
 std::vector<std::vector<int>> cls_by_file = copy_int_list_scheduled_chains(cls);
 std::vector<int> rows0 = copy_rows0_or_empty_scheduled_chains(rows, n);

 const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
 const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());

 int read_nthreads = 1;
#ifdef _OPENMP
 read_nthreads = ncores > 0 ? std::max(1, ncores) : std::max(1, omp_get_max_threads());
#endif

 Rcpp::Rcout
 << "Reading/packing BED with blocked reader: n_bed=" << n
 << ", n_rows=" << n_rows
 << ", read_block_size=" << read_block_size
 << ", read_nthreads=" << read_nthreads
 << "\n";

 FastPackedBedMatrix G = read_bedfiles_to_fast_packed_matrix_blocked(
  bed_files_cpp,
  n,
  rows0_ptr,
  n_rows,
  cls_by_file,
  read_block_size,
  read_nthreads
 );

 const int n_used = G.n;
 const int m = G.m;
 const int nt = y.ncol();

 if (nt <= 0) throw std::runtime_error("y must have at least one trait column.");
 if (y.nrow() != n_used) {
  throw std::runtime_error("y rows must equal the number of samples used after rows filtering.");
 }
 if (static_cast<int>(sets.size()) != m) throw std::runtime_error("sets length must equal number of markers.");
 if (static_cast<int>(b_init.size()) != nt) throw std::runtime_error("b_init must have length nt.");
 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(b_init[static_cast<std::size_t>(t)].size()) != m) {
   throw std::runtime_error("each b_init[t] must have length m.");
  }
 }
 if (pi.size() != 2) throw std::runtime_error("pi must have length 2.");
 if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt) throw std::runtime_error("B must be nt x nt.");
 if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt) throw std::runtime_error("E must be nt x nt.");
 if (static_cast<int>(ssb_prior.size()) != nt || static_cast<int>(sse_prior.size()) != nt) {
  throw std::runtime_error("prior lists must have length nt.");
 }

 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(ssb_prior[t].size()) != nt || static_cast<int>(sse_prior[t].size()) != nt) {
   throw std::runtime_error("priors must be nt x nt.");
  }
  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(t, t2) = ssb_prior[t][t2];
   sse_prior_mat(t, t2) = sse_prior[t][t2];
  }
 }

 std::vector<double> af_cpp = flatten_af_list_or_empty_scheduled_chains(af);
 const bool af_computed = af_cpp.empty();
 if (af_computed) af_cpp = compute_af_from_fast_packed(G);
 if (static_cast<int>(af_cpp.size()) != m) throw std::runtime_error("af must have one value per marker.");

 const int njobs = nchains * nt;
 int nthreads = 1;
#ifdef _OPENMP
 if (ncores > 0) {
  omp_set_dynamic(0);
  omp_set_num_threads(ncores);
 }
 nthreads = ncores > 0 ? std::max(1, std::min(ncores, njobs)) : std::min(omp_get_max_threads(), njobs);
 nthreads = std::max(1, nthreads);
 omp_set_num_threads(nthreads);
#endif

 Rcpp::Rcout
 << "Building scheduled packed BED STBLR chains: n=" << n_used
 << ", m=" << m
 << ", nt=" << nt
 << ", nchains=" << nchains
 << ", jobs=" << njobs
 << ", scale=" << scale
 << ", af_computed=" << af_computed
 << ", read_block_size=" << read_block_size
 << ", progress_every=" << progress_every
 << ", pi_prior_a=" << pi_prior_a
 << ", pi_prior_b=" << pi_prior_b
 << "\n";

 std::vector<MarkerMapSTScheduledChains> marker_maps = build_marker_maps_scheduled_chains(
  G,
  af_cpp,
  scale,
  nthreads
 );
 std::vector<int> marker_order = make_marker_order_from_sets_scheduled_chains(sets, m);

 Rcpp::Rcout
 << "Scheduled chains sampler: full_sweep_every=" << full_sweep_every
 << ", null_skip_base=" << null_skip_base
 << ", null_skip_max=" << null_skip_max
 << ", candidate_threshold=" << candidate_threshold
 << ", candidate_lifetime=" << candidate_lifetime
 << ", skip_nulls_burnin_only=" << skip_nulls_burnin_only
 << ", return_wy=" << return_wy
 << ", return_r=" << return_r
 << ", pi_prior_a=" << pi_prior_a
 << ", pi_prior_b=" << pi_prior_b
 << "\n";

#ifdef _OPENMP
 Rcpp::Rcout
 << "STBLR scheduled packed BED chains OpenMP threads = " << nthreads
 << ", max threads = " << omp_get_max_threads()
 << ", num procs = " << omp_get_num_procs()
 << "\n";
#else
 Rcpp::Rcout << "STBLR scheduled packed BED chains compiled without OpenMP; using one thread.\n";
#endif

 arma::mat y_mat(n_used, nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < n_used; ++i) y_mat(i, t) = y(i, t);
 }

 std::vector<sblr::core::BedScheduledBayesCChainExecutionResult> job_results(
  static_cast<std::size_t>(njobs));

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int job = 0; job < njobs; ++job) {
  const auto task=sblr::core::make_bed_family_task_index(job,nt);
  const int chain=task.chain;
  const int t=task.trait;
  const std::uint64_t chain_seed=
   sblr::core::resolve_bed_family_logical_chain_seed(seed,t,chain);
  const sblr::core::BedPackedGenotypeView<FastPackedBedMatrix> genotype{
   G,G.data.data(),G.data.size(),static_cast<std::size_t>(G.m),
   static_cast<std::size_t>(G.n),G.nbytes,G.stride};
  const sblr::core::ScheduledSweepControl sweep{full_sweep_every,true};
  const sblr::core::NullSkipControl skip{
   null_skip_base,null_skip_max,skip_nulls_burnin_only,"probability_adaptive"};
  const sblr::core::CandidateControl candidate{
   candidate_threshold,candidate_lifetime};
  const sblr::core::BedScheduledBayesCChainExecutionContext<FastPackedBedMatrix> context{
   genotype,marker_maps,marker_order,y_mat,b_init,B,E,ssb_prior_mat,sse_prior_mat,
   pi,sweep,skip,candidate,nub,nue,adjE,pi_prior_a,pi_prior_b,nit,nburn,nthin,
   rebuild_every,progress_every,chain_seed,t,chain,updateB,updateE,updatePi};

  job_results[static_cast<std::size_t>(job)]=
   sblr::core::run_bed_scheduled_bayesc_chain(context);
 }

 int failed_total = 0;
 for (int job = 0; job < njobs; ++job) {
  const auto& r = job_results[static_cast<std::size_t>(job)];
  const auto task=sblr::core::make_bed_family_task_index(job,nt);
  const int chain=task.chain;
  const int t=task.trait;
  Rcpp::Rcout
  << "chain " << chain
  << ", trait " << t
  << ", failed=" << r.failed
  << ", seconds=" << r.seconds
  << ", nsamples=" << r.nsamples
  << ", final_pi=" << r.final_pi
  << ", log_cpo=" << r.log_cpo
  << ", mean_log_cpo=" << r.mean_log_cpo
  << "\n";

  if (r.failed) {
   ++failed_total;
   Rcpp::Rcout
   << "  error: " << r.error << "\n";
  }
 }

 if (failed_total > 0) {
  for (int job = 0; job < njobs; ++job) {
   const auto& r = job_results[static_cast<std::size_t>(job)];
   if (r.failed) {
    const auto task=sblr::core::make_bed_family_task_index(job,nt);
    const int chain=task.chain;
    const int t=task.trait;
    throw std::runtime_error(
      "stblr_cpg_omp_bed_marker_scheduled_chains failed for chain " +
       std::to_string(chain) +
       ", trait " +
       std::to_string(t) +
       ": " +
       r.error
    );
   }
  }
 }

 const sblr::core::BedScheduledBayesCAggregationContext<FastPackedBedMatrix>
  aggregation_context{
   {G,G.data.data(),G.data.size(),static_cast<std::size_t>(G.m),
    static_cast<std::size_t>(G.n),G.nbytes,G.stride},
   marker_maps,marker_order,y_mat,nt,nchains,nit+nburn,return_wy,return_r};
 const auto execution_result=sblr::core::aggregate_bed_scheduled_bayesc_results(
  job_results,aggregation_context);
 const BedScheduledBayesCBindingMetadata binding_metadata{nit,nburn,nthin,n_used};
 return stblr_bed_scheduled_bayesc_result_to_raw(execution_result,binding_metadata);
}
