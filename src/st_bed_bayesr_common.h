#ifndef ST_BED_BAYESR_COMMON_H
#define ST_BED_BAYESR_COMMON_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

struct MarkerMapBayesR {
 double val[4];
 double xx;
};

struct FastPackedBedMatrixBR {
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

static inline unsigned int br_get_bed_code(
  const uint8_t* packed,
  int sample
) {
 return (packed[static_cast<std::size_t>(sample >> 2)] >> (2 * (sample & 3))) & 3u;
}

static FastPackedBedMatrixBR br_read_bed_blocked(
  const std::vector<std::string>& bed_files,
  int n,
  const int* rows0,
  int n_rows,
  const std::vector<std::vector<int>>& cls_by_file,
  int read_block_size,
  int nthreads
) {
 if (n <= 0)
  throw std::runtime_error("br_read_bed_blocked: n must be positive.");
 if (read_block_size <= 0)
  throw std::runtime_error("br_read_bed_blocked: read_block_size must be positive.");
 if (bed_files.size() != cls_by_file.size())
  throw std::runtime_error("br_read_bed_blocked: bed_files and cls lengths differ.");

 const int n_used = n_rows > 0 ? n_rows : n;
 const std::size_t raw_nbytes = (static_cast<std::size_t>(n) + 3u) / 4u;
 const std::size_t out_nbytes = (static_cast<std::size_t>(n_used) + 3u) / 4u;

 int m_total = 0;
 for (std::size_t f = 0; f < cls_by_file.size(); ++f)
  m_total += static_cast<int>(cls_by_file[f].size());
 if (m_total <= 0)
  throw std::runtime_error("br_read_bed_blocked: no markers selected.");

 FastPackedBedMatrixBR G;
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
   if (r < 0 || r >= n)
    throw std::runtime_error("br_read_bed_blocked: rows0 contains index outside [0, n).");
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
  if (!fs)
   throw std::runtime_error("Could not open BED file: " + bed_files[f]);

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
     throw std::runtime_error("br_read_bed_blocked: cls contains non-positive marker index.");
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

 if (global_marker != m_total)
  throw std::runtime_error("br_read_bed_blocked: internal marker count mismatch.");

 return G;
}

static std::vector<double> br_compute_af(const FastPackedBedMatrixBR& G) {
 std::vector<double> af(static_cast<std::size_t>(G.m), 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
 for (int j = 0; j < G.m; ++j) {
  const uint8_t* row = G.row(j);
  double dosage_sum = 0.0;
  double allele_count = 0.0;

  for (int i = 0; i < G.n; ++i) {
   const unsigned int code = br_get_bed_code(row, i);
   if (code == 1u) continue;
   if (code == 0u) dosage_sum += 2.0;
   else if (code == 2u) dosage_sum += 1.0;
   allele_count += 2.0;
  }

  af[static_cast<std::size_t>(j)] =
   allele_count > 0.0 ? dosage_sum / allele_count : 0.0;
 }

 return af;
}

static std::vector<std::string> br_copy_bed_files(Rcpp::CharacterVector bed_files) {
 std::vector<std::string> out(static_cast<std::size_t>(bed_files.size()));
 for (int i = 0; i < bed_files.size(); ++i)
  out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(bed_files[i]);
 return out;
}

static std::vector<std::vector<int>> br_copy_int_list(Rcpp::List xlist) {
 std::vector<std::vector<int>> out(static_cast<std::size_t>(xlist.size()));
 for (int f = 0; f < xlist.size(); ++f) {
  Rcpp::IntegerVector x = xlist[f];
  out[static_cast<std::size_t>(f)].resize(static_cast<std::size_t>(x.size()));
  for (int i = 0; i < x.size(); ++i)
   out[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)] = x[i];
 }
 return out;
}

static std::vector<int> br_copy_rows0(
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
    if (out[static_cast<std::size_t>(i)] != i) { identity = false; break; }
   }
   if (identity) out.clear();
  }
 }
 return out;
}

static std::vector<double> br_flatten_af_list(Rcpp::Nullable<Rcpp::List> af) {
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

static inline void br_marker_map_values(
  double p,
  bool scale,
  MarkerMapBayesR& map
) {
 if (scale) {
  const double denom = std::sqrt(2.0 * p * (1.0 - p));
  if (denom <= 0.0 || !std::isfinite(denom)) {
   map.val[0] = map.val[1] = map.val[2] = map.val[3] = 0.0;
  } else {
   // PLINK 2-bit BED coding: 0->genotype 2, 1->missing, 2->genotype 1, 3->genotype 0
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

static inline double br_marker_xx(
  const FastPackedBedMatrixBR& G,
  int marker,
  const MarkerMapBayesR& map
) {
 const uint8_t* packed = G.row(marker);
 const int n = G.n;
 const std::size_t nbytes = G.nbytes;
 double xx = 0.0;

 for (std::size_t kb = 0; kb < nbytes; ++kb) {
  const unsigned char byte = packed[kb];
  const int jbase = static_cast<int>(kb << 2);
  if (jbase + 0 < n) { const double x = map.val[(byte >> 0) & 3u]; xx += x * x; }
  if (jbase + 1 < n) { const double x = map.val[(byte >> 2) & 3u]; xx += x * x; }
  if (jbase + 2 < n) { const double x = map.val[(byte >> 4) & 3u]; xx += x * x; }
  if (jbase + 3 < n) { const double x = map.val[(byte >> 6) & 3u]; xx += x * x; }
 }

 return xx;
}

static std::vector<MarkerMapBayesR> br_build_marker_maps(
  const FastPackedBedMatrixBR& G,
  const std::vector<double>& af,
  bool scale,
  int nthreads
) {
 const int m = G.m;
 if (static_cast<int>(af.size()) != m)
  throw std::runtime_error("br_build_marker_maps: af length must equal number of markers.");

 std::vector<MarkerMapBayesR> maps(static_cast<std::size_t>(m));

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
 for (int j = 0; j < m; ++j) {
  MarkerMapBayesR map;
  br_marker_map_values(af[static_cast<std::size_t>(j)], scale, map);
  map.xx = br_marker_xx(G, j, map);
  maps[static_cast<std::size_t>(j)] = map;
 }

 for (int j = 0; j < m; ++j) {
  if (!std::isfinite(maps[static_cast<std::size_t>(j)].xx) ||
      maps[static_cast<std::size_t>(j)].xx <= 0.0) {
   throw std::runtime_error(
     "br_build_marker_maps: invalid x'x for marker " + std::to_string(j) +
      ". Check allele frequencies and monomorphic markers."
   );
  }
 }

 return maps;
}

static inline double br_dot_residual(
  const FastPackedBedMatrixBR& G,
  int marker,
  const MarkerMapBayesR& map,
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

static inline void br_update_residual(
  const FastPackedBedMatrixBR& G,
  int marker,
  const MarkerMapBayesR& map,
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

static std::vector<int> br_make_marker_order(
  const std::vector<int>& sets,
  int m
) {
 if (static_cast<int>(sets.size()) != m)
  throw std::runtime_error("br_make_marker_order: sets length must equal m.");

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
 for (std::size_t b = 0; b < block_markers.size(); ++b)
  for (int j : block_markers[b]) order.push_back(j);
 return order;
}

static inline void sampleB_bayesr(
  int m,
  const std::vector<double>& c,
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
  const int di = d(iu);
  if (di > 0) {
   const double ck = c[static_cast<std::size_t>(di)];
   ssb += b(iu) * b(iu) / ck;
   dfb += 1.0;
  }
 }

 const double scale = ssb + nub * ssb_prior;
 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 vb = std::max(scale / chi2, 1e-12);
}

static inline void sampleE_bayesr(
  double nue,
  double& ve,
  const arma::vec& e,
  double sse_prior,
  std::mt19937& gen
) {
 const double sse = arma::dot(e, e);
 const double scale = sse + nue * sse_prior;

 if (!std::isfinite(scale) || scale <= 0.0)
  throw std::runtime_error("sampleE_bayesr: invalid residual scale.");

 std::chi_squared_distribution<double> rchisq(static_cast<double>(e.n_elem) + nue);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 const double ve_new = scale / chi2;

 if (!std::isfinite(ve_new) || ve_new <= 0.0)
  throw std::runtime_error("sampleE_bayesr: sampled ve is invalid.");

 ve = std::max(ve_new, 1e-12);
}

static inline void br_update_log_inv_cpo(
  const arma::vec& e,
  double ve,
  arma::vec& log_inv_cpo
) {
 const double ve_safe = std::max(ve, 1e-300);
 const double log_norm = -0.5 * std::log(2.0 * std::acos(-1.0) * ve_safe);
 const arma::uword n = e.n_elem;

 if (log_inv_cpo.n_elem != n)
  throw std::runtime_error("br_update_log_inv_cpo: dimension mismatch.");

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

static inline double br_compute_total_log_cpo(
  const arma::vec& log_inv_cpo,
  double nsamples
) {
 if (nsamples <= 0.0) return NA_REAL;
 const double log_nsamples = std::log(nsamples);
 double total = 0.0;
 for (arma::uword i = 0; i < log_inv_cpo.n_elem; ++i) {
  if (std::isfinite(log_inv_cpo(i)))
   total += log_nsamples - log_inv_cpo(i);
 }
 return total;
}

static inline double br_computeG(const arma::vec& y, const arma::vec& e) {
 double ss = 0.0;
 const int n = static_cast<int>(y.n_elem);
 for (int i = 0; i < n; ++i) {
  const double g = y(static_cast<arma::uword>(i)) - e(static_cast<arma::uword>(i));
  ss += g * g;
 }
 return ss / static_cast<double>(n);
}

static inline double br_computeLE(
  const arma::rowvec& b,
  const std::vector<MarkerMapBayesR>& maps,
  int n
) {
 double vle = 0.0;
 const int m = static_cast<int>(b.n_elem);
 const double inv_n = 1.0 / static_cast<double>(n);
 for (int j = 0; j < m; ++j) {
  const double bj = b(static_cast<arma::uword>(j));
  if (bj != 0.0)
   vle += maps[static_cast<std::size_t>(j)].xx * inv_n * bj * bj;
 }
 return vle;
}

static arma::vec br_xb(
  const FastPackedBedMatrixBR& G,
  const std::vector<MarkerMapBayesR>& maps,
  const std::vector<int>& marker_order,
  const arma::rowvec& b
) {
 arma::vec xb(G.n, arma::fill::zeros);
 for (int marker : marker_order) {
  const double bj = b(static_cast<arma::uword>(marker));
  if (bj != 0.0)
   br_update_residual(G, marker, maps[static_cast<std::size_t>(marker)], xb.memptr(), -bj);
 }
 return xb;
}

#endif

