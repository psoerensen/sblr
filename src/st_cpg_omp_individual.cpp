// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "packed_bed.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
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
// Packed BED marker-wise ST BayesC sampler with optional null-marker skipping
// =============================================================================
//
// Key idea:
//   score_j = x_j' e + x_j'x_j b_j
//   if effect changes by diff:
//      e <- e - x_j diff
//
// The functions below compute x_j'e and update e directly from the packed BED
// bytes, without decoding full marker/block matrices.
//
// Null skipping is optional. For exact full-sweep Gibbs sampling set either:
//   full_sweep_every = 1
// or
//   null_update_prob = 1.0
// or
//   skip_nulls_burnin_only = true and collect only post-burn-in samples.
// =============================================================================

struct MarkerMapST {
 double m0;
 double m1;
 double m2;
 double m3;
 double xx;
};

static std::vector<std::string> stblr_copy_bed_files_sparse(
  Rcpp::CharacterVector bed_files
) {
 std::vector<std::string> out(static_cast<std::size_t>(bed_files.size()));

 for (int i = 0; i < bed_files.size(); ++i) {
  out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(bed_files[i]);
 }

 return out;
}

static std::vector<std::vector<int>> stblr_copy_int_list_sparse(
  Rcpp::List xlist
) {
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

static std::vector<int> stblr_copy_rows0_or_empty_sparse(
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  int n
) {
 std::vector<int> out;

 if (rows.isNotNull()) {
  Rcpp::IntegerVector r(rows);
  out.resize(static_cast<std::size_t>(r.size()));

  for (int i = 0; i < r.size(); ++i) {
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

static std::vector<double> stblr_flatten_af_list_or_empty_sparse(
  Rcpp::Nullable<Rcpp::List> af
) {
 std::vector<double> out;

 if (af.isNotNull()) {
  Rcpp::List af_list = Rcpp::as<Rcpp::List>(af.get());

  for (int f = 0; f < af_list.size(); ++f) {
   Rcpp::NumericVector x = af_list[f];

   for (int i = 0; i < x.size(); ++i) {
    out.push_back(x[i]);
   }
  }
 }

 return out;
}

static inline void stblr_bed_maps_sparse(
  double p,
  bool scale,
  double& m0,
  double& m1,
  double& m2,
  double& m3
) {
 if (scale) {
  const double denom = std::sqrt(2.0 * p * (1.0 - p));

  if (denom <= 0.0 || !std::isfinite(denom)) {
   m0 = 0.0;
   m1 = 0.0;
   m2 = 0.0;
   m3 = 0.0;
  } else {
   // PLINK BED codes used here:
   //   0 -> genotype 2
   //   1 -> missing, imputed to mean
   //   2 -> genotype 1
   //   3 -> genotype 0
   m0 = (2.0 - 2.0 * p) / denom;
   m1 = 0.0;
   m2 = (1.0 - 2.0 * p) / denom;
   m3 = (0.0 - 2.0 * p) / denom;
  }
 } else {
  m0 = 2.0;
  m1 = 2.0 * p;
  m2 = 1.0;
  m3 = 0.0;
 }
}

static inline double stblr_decode_2bit_value_sparse(
  unsigned int code,
  const MarkerMapST& map
) {
 switch (code) {
 case 0u: return map.m0;
 case 1u: return map.m1;
 case 2u: return map.m2;
 default: return map.m3;
 }
}

static double stblr_marker_xx_from_packed(
  const PackedBedMatrix& G,
  int marker,
  const MarkerMapST& map
) {
 const uint8_t* packed = G.row(marker);
 const int n = G.n;
 const std::size_t nbytes = G.nbytes;

 double xx = 0.0;

 for (std::size_t kb = 0; kb < nbytes; ++kb) {
  const unsigned char byte = packed[kb];
  const int jbase = static_cast<int>(kb << 2);

  if (jbase + 0 < n) {
   const double x = stblr_decode_2bit_value_sparse((byte >> 0) & 3u, map);
   xx += x * x;
  }
  if (jbase + 1 < n) {
   const double x = stblr_decode_2bit_value_sparse((byte >> 2) & 3u, map);
   xx += x * x;
  }
  if (jbase + 2 < n) {
   const double x = stblr_decode_2bit_value_sparse((byte >> 4) & 3u, map);
   xx += x * x;
  }
  if (jbase + 3 < n) {
   const double x = stblr_decode_2bit_value_sparse((byte >> 6) & 3u, map);
   xx += x * x;
  }
 }

 return xx;
}

static std::vector<MarkerMapST> stblr_build_marker_maps_sparse(
  const PackedBedMatrix& G,
  const std::vector<double>& af,
  bool scale
) {
 const int m = G.m;

 if (static_cast<int>(af.size()) != m) {
  throw std::runtime_error("stblr_build_marker_maps_sparse: af length must equal number of markers.");
 }

 std::vector<MarkerMapST> maps(static_cast<std::size_t>(m));

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
 for (int j = 0; j < m; ++j) {
  MarkerMapST map;

  stblr_bed_maps_sparse(
   af[static_cast<std::size_t>(j)],
     scale,
     map.m0,
     map.m1,
     map.m2,
     map.m3
  );

  map.xx = stblr_marker_xx_from_packed(G, j, map);

  if (!std::isfinite(map.xx) || map.xx <= 0.0) {
   // Mark invalid markers by zero xx. The main sampler will throw with index.
   map.xx = 0.0;
  }

  maps[static_cast<std::size_t>(j)] = map;
 }

 for (int j = 0; j < m; ++j) {
  if (!std::isfinite(maps[static_cast<std::size_t>(j)].xx) ||
      maps[static_cast<std::size_t>(j)].xx <= 0.0) {
   throw std::runtime_error(
     "stblr_build_marker_maps_sparse: invalid x'x for marker " + std::to_string(j) +
      ". Check allele frequencies and monomorphic markers."
   );
  }
 }

 return maps;
}

static inline double stblr_bed_marker_dot_residual(
  const PackedBedMatrix& G,
  int marker,
  const MarkerMapST& map,
  const double* e
) {
 const uint8_t* packed = G.row(marker);
 const int n = G.n;
 const std::size_t nbytes = G.nbytes;

 double out = 0.0;

 for (std::size_t kb = 0; kb < nbytes; ++kb) {
  const unsigned char byte = packed[kb];
  const int jbase = static_cast<int>(kb << 2);

  if (jbase + 0 < n) {
   out += stblr_decode_2bit_value_sparse((byte >> 0) & 3u, map) * e[jbase + 0];
  }
  if (jbase + 1 < n) {
   out += stblr_decode_2bit_value_sparse((byte >> 2) & 3u, map) * e[jbase + 1];
  }
  if (jbase + 2 < n) {
   out += stblr_decode_2bit_value_sparse((byte >> 4) & 3u, map) * e[jbase + 2];
  }
  if (jbase + 3 < n) {
   out += stblr_decode_2bit_value_sparse((byte >> 6) & 3u, map) * e[jbase + 3];
  }
 }

 return out;
}

static inline void stblr_bed_marker_update_residual(
  const PackedBedMatrix& G,
  int marker,
  const MarkerMapST& map,
  double* e,
  double diff
) {
 const uint8_t* packed = G.row(marker);
 const int n = G.n;
 const std::size_t nbytes = G.nbytes;

 for (std::size_t kb = 0; kb < nbytes; ++kb) {
  const unsigned char byte = packed[kb];
  const int jbase = static_cast<int>(kb << 2);

  if (jbase + 0 < n) {
   e[jbase + 0] -= stblr_decode_2bit_value_sparse((byte >> 0) & 3u, map) * diff;
  }
  if (jbase + 1 < n) {
   e[jbase + 1] -= stblr_decode_2bit_value_sparse((byte >> 2) & 3u, map) * diff;
  }
  if (jbase + 2 < n) {
   e[jbase + 2] -= stblr_decode_2bit_value_sparse((byte >> 4) & 3u, map) * diff;
  }
  if (jbase + 3 < n) {
   e[jbase + 3] -= stblr_decode_2bit_value_sparse((byte >> 6) & 3u, map) * diff;
  }
 }
}

static std::vector<int> stblr_make_marker_order_from_sets_sparse(
  const std::vector<int>& sets,
  int m
) {
 if (static_cast<int>(sets.size()) != m) {
  throw std::runtime_error("stblr_make_marker_order_from_sets_sparse: sets length must equal m.");
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
  for (int j : block_markers[b]) {
   order.push_back(j);
  }
 }

 return order;
}

// -----------------------------------------------------------------------------
// Single-marker packed BED BayesC update
// -----------------------------------------------------------------------------

static inline double sampleBetaC_ST_bed_marker(
  const PackedBedMatrix& G,
  int marker,
  const MarkerMapST& map,
  const std::vector<double>& pi,
  double vb,
  double vei,
  arma::vec& e,
  double& b_j,
  int& d_j,
  std::mt19937& gen
) {
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm01(0.0, 1.0);

 const double xj2 = map.xx;
 const double vei_safe = std::max(vei, 1e-300);
 const double pi0 = std::max(pi[0], 1e-300);
 const double pi1 = std::max(pi[1], 1e-300);

 const double xte = stblr_bed_marker_dot_residual(
  G,
  marker,
  map,
  e.memptr()
 );

 const double score = xte + xj2 * b_j;
 const double denom = std::max(vei_safe + xj2 * vb, 1e-300);

 const double logBF =
  0.5 * std::log(vei_safe / denom) +
  0.5 * score * score * vb / (vei_safe * denom);

 const double logp1 = std::log(pi1) + logBF;
 const double logp0 = std::log(pi0);
 const double delta_log = logp0 - logp1;

 double p1 = 0.0;

 if (delta_log > 35.0) {
  p1 = 0.0;
 } else if (delta_log < -35.0) {
  p1 = 1.0;
 } else {
  p1 = 1.0 / (1.0 + std::exp(delta_log));
 }

 const int d_new = (runif(gen) < p1) ? 1 : 0;

 double b_new = 0.0;

 if (d_new == 1) {
  const double lhs = xj2 + vei_safe / vb;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);

  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b_j;

 if (diff != 0.0) {
  stblr_bed_marker_update_residual(
   G,
   marker,
   map,
   e.memptr(),
   diff
  );
 }

 b_j = b_new;
 d_j = d_new;

 return p1;
}

// -----------------------------------------------------------------------------
// Hyperparameter updates
// -----------------------------------------------------------------------------

static inline void sampleB_ST_bed_sparse(
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

static inline void sampleE_ST_bed_sparse(
  double nue,
  double& ve,
  const arma::vec& e,
  double sse_prior,
  std::mt19937& gen
) {
 const double sse = arma::dot(e, e);
 const double scale = sse + nue * sse_prior;

 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleE_ST_bed_sparse: invalid residual scale.");
 }

 std::chi_squared_distribution<double> rchisq(static_cast<double>(e.n_elem) + nue);
 const double chi2 = std::max(rchisq(gen), 1e-300);

 const double ve_new = scale / chi2;

 if (!std::isfinite(ve_new) || ve_new <= 0.0) {
  throw std::runtime_error("sampleE_ST_bed_sparse: sampled ve is invalid.");
 }

 ve = std::max(ve_new, 1e-12);
}

static inline double computeG_ST_bed_sparse(
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

static inline void samplePi_ST_bed_sparse(
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

static arma::vec stblr_bed_xb_from_b_sparse(
  const PackedBedMatrix& G,
  const std::vector<MarkerMapST>& maps,
  const std::vector<int>& marker_order,
  const arma::rowvec& b
) {
 arma::vec xb(G.n, arma::fill::zeros);

 for (int marker : marker_order) {
  const arma::uword ju = static_cast<arma::uword>(marker);
  const double bj = b(ju);

  if (bj != 0.0) {
   // xb += x_j * bj. This is the opposite sign of the residual update.
   stblr_bed_marker_update_residual(
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

// =============================================================================
// Exported packed BED marker-wise sampler
// =============================================================================

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_bed_marker_sparse(
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
  double null_update_prob,
  double candidate_threshold,
  int candidate_lifetime,
  bool skip_nulls_burnin_only,
  int ncores,
  int seed
) {
 if (nit <= 0 || nburn < 0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: nit must be positive and nburn non-negative.");
 }

 if (nthin <= 0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: nthin must be positive.");
 }

 if (rebuild_every < 0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: rebuild_every must be >= 0. Use 0 to disable residual rebuilding.");
 }

 if (full_sweep_every < 0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: full_sweep_every must be >= 0.");
 }

 if (!std::isfinite(null_update_prob) || null_update_prob < 0.0 || null_update_prob > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: null_update_prob must be in [0,1].");
 }

 if (!std::isfinite(candidate_threshold) || candidate_threshold < 0.0 || candidate_threshold > 1.0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: candidate_threshold must be in [0,1].");
 }

 if (candidate_lifetime < 0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: candidate_lifetime must be >= 0.");
 }

 std::vector<std::string> bed_files_cpp = stblr_copy_bed_files_sparse(bed_files);
 std::vector<std::vector<int>> cls_by_file = stblr_copy_int_list_sparse(cls);
 std::vector<int> rows0 = stblr_copy_rows0_or_empty_sparse(rows, n);

 const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
 const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());

 PackedBedMatrix G = read_bedfiles_to_packed_matrix(
  bed_files_cpp,
  n,
  rows0_ptr,
  n_rows,
  cls_by_file
 );

 const int n_used = G.n;
 const int m = G.m;
 const int nt = y.ncol();

 if (nt <= 0) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: y must have at least one trait column.");
 }

 if (y.nrow() != n_used) {
  throw std::runtime_error(
    "stblr_cpg_omp_bed_marker_sparse: y rows must equal the number of samples used after rows filtering. "
    "If rows is supplied, pass y already subset to those rows."
  );
 }

 if (static_cast<int>(sets.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: sets must have length equal to the total number of BED markers used.");
 }

 if (static_cast<int>(b_init.size()) != nt) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: b_init must have length nt.");
 }

 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(b_init[static_cast<std::size_t>(t)].size()) != m) {
   throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: each b_init[t] must have length m.");
  }
 }

 if (pi.size() != 2) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: pi must have length 2.");
 }

 if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: B must be nt x nt.");
 }

 if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: E must be nt x nt.");
 }

 if (static_cast<int>(ssb_prior.size()) != nt || static_cast<int>(sse_prior.size()) != nt) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: prior lists must have length nt.");
 }

 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(ssb_prior[t].size()) != nt ||
      static_cast<int>(sse_prior[t].size()) != nt) {
   throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: priors must be nt x nt.");
  }

  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(t, t2) = ssb_prior[t][t2];
   sse_prior_mat(t, t2) = sse_prior[t][t2];
  }
 }

 std::vector<double> af_cpp = stblr_flatten_af_list_or_empty_sparse(af);
 const bool af_computed = af_cpp.empty();

 if (af_computed) {
  af_cpp = compute_af_from_packed(G);
 }

 if (static_cast<int>(af_cpp.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_bed_marker_sparse: af must have one value per marker after flattening.");
 }

#ifdef _OPENMP
 if (ncores > 0) {
  omp_set_dynamic(0);
  omp_set_num_threads(ncores);
 }
#endif

 Rcpp::Rcout
 << "Building packed BED marker-wise STBLR: n=" << n_used
 << ", m=" << m
 << ", nt=" << nt
 << ", scale=" << scale
 << ", af_computed=" << af_computed
 << "\n";

 std::vector<MarkerMapST> marker_maps = stblr_build_marker_maps_sparse(
  G,
  af_cpp,
  scale
 );

 std::vector<int> marker_order = stblr_make_marker_order_from_sets_sparse(
  sets,
  m
 );

 Rcpp::Rcout
 << "Packed BED marker-wise sampler: full_sweep_every=" << full_sweep_every
 << ", null_update_prob=" << null_update_prob
 << ", candidate_threshold=" << candidate_threshold
 << ", candidate_lifetime=" << candidate_lifetime
 << ", skip_nulls_burnin_only=" << skip_nulls_burnin_only
 << "\n";

 arma::mat y_mat(n_used, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < n_used; ++i) {
   y_mat(static_cast<arma::uword>(i), static_cast<arma::uword>(t)) = y(i, t);
  }
 }

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::Mat<int> d_mat(nt, m, arma::fill::zeros);

 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);

 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);

 std::vector<int> failed(static_cast<std::size_t>(nt), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(nt));
 std::vector<int> thread_used(static_cast<std::size_t>(nt), 0);
 std::vector<double> trait_seconds(static_cast<std::size_t>(nt), 0.0);

 int nthreads = 1;
#ifdef _OPENMP
 nthreads = ncores > 0 ? std::max(1, std::min(ncores, nt)) : std::min(omp_get_max_threads(), nt);
 omp_set_num_threads(nthreads);
 Rcpp::Rcout
 << "STBLR packed BED OpenMP threads = " << nthreads
 << ", max threads = " << omp_get_max_threads()
 << ", num procs = " << omp_get_num_procs()
 << "\n";
#endif

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
   std::uniform_real_distribution<double> runif_skip(0.0, 1.0);

   arma::vec y_t = y_mat.col(static_cast<arma::uword>(t));
   arma::rowvec b_t(m, arma::fill::zeros);
   arma::Row<int> d_t(m, arma::fill::zeros);

   for (int j = 0; j < m; ++j) {
    b_t(static_cast<arma::uword>(j)) = b_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(j)];
    d_t(static_cast<arma::uword>(j)) = (b_t(static_cast<arma::uword>(j)) != 0.0) ? 1 : 0;
   }

   arma::vec xb_t = stblr_bed_xb_from_b_sparse(
    G,
    marker_maps,
    marker_order,
    b_t
   );

   arma::vec e_t = y_t - xb_t;

   double vb_t = B(t, t);
   double ve_t = E(t, t);
   double vg_t = computeG_ST_bed_sparse(y_t, e_t);
   double vei_t = ve_t + adjE * vg_t;

   std::vector<double> pi_t = pi;

   const double psum = pi_t[0] + pi_t[1];
   if (!std::isfinite(psum) || psum <= 0.0 || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
    throw std::runtime_error("invalid initial pi.");
   }
   pi_t[0] /= psum;
   pi_t[1] /= psum;

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);

   std::vector<unsigned char> candidate(static_cast<std::size_t>(m), 0u);
   std::vector<int> last_interesting(static_cast<std::size_t>(m), -1000000000);

   double nsamples_t = 0.0;

   for (int it = 0; it < nit + nburn; ++it) {
    const bool skipping_allowed =
     (full_sweep_every > 1 || null_update_prob < 1.0) &&
     (!skip_nulls_burnin_only || it < nburn);

    const bool full_sweep =
     !skipping_allowed ||
     full_sweep_every <= 0 ||
     ((it % full_sweep_every) == 0);

    for (int jj = 0; jj < m; ++jj) {
     const int marker = marker_order[static_cast<std::size_t>(jj)];
     const arma::uword ju = static_cast<arma::uword>(marker);

     bool update_marker = full_sweep;

     if (!update_marker && d_t(ju) > 0) {
      update_marker = true;
     }

     if (!update_marker && candidate[static_cast<std::size_t>(marker)] != 0u) {
      update_marker = true;
     }

     if (!update_marker && null_update_prob > 0.0 && runif_skip(gen_t) < null_update_prob) {
      update_marker = true;
     }

     if (!update_marker) {
      continue;
     }

     double bj = b_t(ju);
     int dj = d_t(ju);

     const double p1 = sampleBetaC_ST_bed_marker(
      G,
      marker,
      marker_maps[static_cast<std::size_t>(marker)],
                 pi_t,
                 vb_t,
                 vei_t,
                 e_t,
                 bj,
                 dj,
                 gen_t
     );

     b_t(ju) = bj;
     d_t(ju) = dj;

     if (dj > 0 || p1 >= candidate_threshold) {
      candidate[static_cast<std::size_t>(marker)] = 1u;
      last_interesting[static_cast<std::size_t>(marker)] = it;
     } else {
      if (candidate[static_cast<std::size_t>(marker)] != 0u &&
          it - last_interesting[static_cast<std::size_t>(marker)] > candidate_lifetime) {
       candidate[static_cast<std::size_t>(marker)] = 0u;
      }
     }
    }

    if (updateE && rebuild_every > 0 && ((it + 1) % rebuild_every == 0)) {
     xb_t = stblr_bed_xb_from_b_sparse(
      G,
      marker_maps,
      marker_order,
      b_t
     );
     e_t = y_t - xb_t;
    }

    if (updateB) {
     sampleB_ST_bed_sparse(
      m,
      nub,
      vb_t,
      b_t,
      d_t,
      ssb_prior_mat(t, t),
      gen_t
     );

     if (!std::isfinite(vb_t) || vb_t <= 0.0) {
      throw std::runtime_error("vb became invalid after sampleB.");
     }
    }

    if (updateE) {
     sampleE_ST_bed_sparse(
      nue,
      ve_t,
      e_t,
      sse_prior_mat(t, t),
      gen_t
     );

     if (!std::isfinite(ve_t) || ve_t <= 0.0) {
      throw std::runtime_error("ve became invalid after sampleE.");
     }
    }

    if (updatePi) {
     samplePi_ST_bed_sparse(d_t, pi_t, gen_t);

     if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
         pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
      throw std::runtime_error("pi became invalid after samplePi.");
     }
    }

    vg_t = computeG_ST_bed_sparse(y_t, e_t);
    vei_t = ve_t + adjE * vg_t;

    if (!std::isfinite(vg_t)) {
     throw std::runtime_error("vg became NaN/Inf after computeG.");
    }

    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error("adjusted residual variance became invalid.");
    }

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    pis_t(static_cast<arma::uword>(it)) = pi_t[1];

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;

     for (int j = 0; j < m; ++j) {
      const arma::uword ju = static_cast<arma::uword>(j);
      bm_t(ju) += b_t(ju);
      dm_t(ju) += static_cast<double>(d_t(ju));
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
   b_mat.row(static_cast<arma::uword>(t)) = b_t;
   d_mat.row(static_cast<arma::uword>(t)) = d_t;

   vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
   vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
   ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
   pis_mat.row(static_cast<arma::uword>(t)) = pis_t;

   final_vb(static_cast<arma::uword>(t)) = vb_t;
   final_ve(static_cast<arma::uword>(t)) = ve_t;
   final_vg(static_cast<arma::uword>(t)) = vg_t;
   final_pi(static_cast<arma::uword>(t)) = pi_t[1];
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
     "stblr_cpg_omp_bed_marker_sparse failed for trait " +
      std::to_string(t) +
      ": " +
      errors[static_cast<std::size_t>(t)]
   );
  }
 }

 // --------------------------------------------------------------------------
 // Final wy and r from packed BED.
 // --------------------------------------------------------------------------

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  arma::vec y_t = y_mat.col(static_cast<arma::uword>(t));
  arma::rowvec b_t = b_mat.row(static_cast<arma::uword>(t));

  arma::vec xb_t = stblr_bed_xb_from_b_sparse(
   G,
   marker_maps,
   marker_order,
   b_t
  );

  arma::vec e_t = y_t - xb_t;

  for (int j = 0; j < m; ++j) {
   wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(j)) =
    stblr_bed_marker_dot_residual(
     G,
     j,
     marker_maps[static_cast<std::size_t>(j)],
                y_t.memptr()
    );

   r_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(j)) =
    stblr_bed_marker_dot_residual(
     G,
     j,
     marker_maps[static_cast<std::size_t>(j)],
                e_t.memptr()
    );
  }
 }

 // --------------------------------------------------------------------------
 // Build result, preserving existing slots and adding the pi trace in slot 23.
 // --------------------------------------------------------------------------

 std::vector<std::vector<std::vector<double>>> result(23);

 for (int k = 0; k < 23; ++k) {
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
  result[18][static_cast<std::size_t>(t)].resize(4);
  result[19][static_cast<std::size_t>(t)].resize(2);
  result[22][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int j = 0; j < m; ++j) {
   const std::size_t js = static_cast<std::size_t>(j);
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword ju = static_cast<arma::uword>(j);

   result[0][ts][js] = bm_mat(tu, ju);
   result[1][ts][js] = dm_mat(tu, ju);
   result[2][ts][js] = wy_mat(tu, ju);
   result[3][ts][js] = r_mat(tu, ju);
   result[4][ts][js] = b_mat(tu, ju);
   result[5][ts][js] = static_cast<double>(d_mat(tu, ju));
   result[6][ts][js] = static_cast<double>(j);
  }
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int it = 0; it < nit + nburn; ++it) {
   const std::size_t its = static_cast<std::size_t>(it);
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword itu = static_cast<arma::uword>(it);

   result[7][ts][its] = vbs_mat(tu, itu);
   result[8][ts][its] = vgs_mat(tu, itu);
   result[9][ts][its] = ves_mat(tu, itu);
   result[22][ts][its] = pis_mat(tu, itu);
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

  if (npi > 0) {
   mean_pi /= static_cast<double>(npi);
  } else {
   mean_pi = final_pi(static_cast<arma::uword>(t));
  }

  result[17][ts][0] = 1.0 - mean_pi;
  result[17][ts][1] = mean_pi;

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
// #include <limits>
// #include <numeric>
// #include <random>
// #include <stdexcept>
// #include <vector>
// #include <unordered_map>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // -----------------------------------------------------------------------------
// // Individual-level block structures
// // -----------------------------------------------------------------------------
//
// struct STBlock {
//  std::vector<int> markers;  // global marker indices, 0-based
//  arma::fmat XtX;            // float block X'X, size p x p
//  arma::vec Xty;             // block X'y, size p
//  arma::vec xx;              // double diag(X'X), size p
// };
//
// inline std::vector<STBlock> build_st_blocks_individual(
//   const arma::mat& X,
//   const arma::vec& y,
//   const std::vector<int>& sets
// ) {
//  const int n = static_cast<int>(X.n_rows);
//  const int m = static_cast<int>(X.n_cols);
//
//  if (static_cast<int>(y.n_elem) != n) {
//   throw std::runtime_error("build_st_blocks_individual: y length must equal nrow(X).");
//  }
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("build_st_blocks_individual: sets must have length ncol(X). It should assign each marker to a block.");
//  }
//
//  std::vector<int> labels;
//  labels.reserve(m);
//  std::unordered_map<int, int> label_to_block;
//
//  for (int j = 0; j < m; ++j) {
//   const int lab = sets[j];
//
//   if (label_to_block.find(lab) == label_to_block.end()) {
//    const int block_id = static_cast<int>(labels.size());
//    label_to_block[lab] = block_id;
//    labels.push_back(lab);
//   }
//  }
//
//  std::vector<STBlock> blocks(labels.size());
//
//  for (int j = 0; j < m; ++j) {
//   blocks[label_to_block[sets[j]]].markers.push_back(j);
//  }
//
//  for (STBlock& blk : blocks) {
//   const int p = static_cast<int>(blk.markers.size());
//   arma::uvec idx(p);
//
//   for (int k = 0; k < p; ++k) {
//    idx(k) = static_cast<arma::uword>(blk.markers[k]);
//   }
//
//   arma::mat Xb = X.cols(idx);
//   arma::mat XtX_d = Xb.t() * Xb;
//
//   blk.XtX = arma::conv_to<arma::fmat>::from(XtX_d);
//   blk.Xty = Xb.t() * y;
//   blk.xx  = arma::diagvec(XtX_d);
//
//   for (int k = 0; k < p; ++k) {
//    if (!std::isfinite(blk.xx(k)) || blk.xx(k) <= 0.0) {
//     throw std::runtime_error(
//       "build_st_blocks_individual: marker has invalid X'X diagonal in block."
//     );
//    }
//   }
//  }
//
//  return blocks;
// }
//
// inline arma::vec block_xb_from_global_b(
//   const arma::mat& X,
//   const std::vector<int>& markers,
//   const arma::rowvec& b
// ) {
//  const int p = static_cast<int>(markers.size());
//  arma::vec xb(X.n_rows, arma::fill::zeros);
//
//  for (int k = 0; k < p; ++k) {
//   const int j = markers[k];
//   const double bj = b(j);
//
//   if (bj != 0.0) {
//    xb += X.col(j) * bj;
//   }
//  }
//
//  return xb;
// }
//
// inline arma::vec compute_xb_global(
//   const arma::mat& X,
//   const arma::rowvec& b
// ) {
//  return X * b.t();
// }
//
// inline arma::vec xtx_fmat_times_vec(
//   const arma::fmat& XtX,
//   const arma::vec& b
// ) {
//  const int p = static_cast<int>(b.n_elem);
//  arma::vec out(p, arma::fill::zeros);
//
//  for (int j = 0; j < p; ++j) {
//   const double bj = b(j);
//
//   if (bj == 0.0) continue;
//
//   for (int i = 0; i < p; ++i) {
//    out(i) += static_cast<double>(XtX(i, j)) * bj;
//   }
//  }
//
//  return out;
// }
//
// // -----------------------------------------------------------------------------
// // Block-local single-trait BayesC marker update
// // -----------------------------------------------------------------------------
//
// inline void sampleBetaC_ST_block_dense(
//   int k,
//   const std::vector<double>& pi,
//   double vb,
//   double vei,
//   const arma::vec& xx,
//   const arma::fmat& XtX,
//   arma::vec& r,
//   arma::vec& b_block,
//   arma::Row<int>& d_block,
//   std::mt19937& gen
// ) {
//  const double xk2 = xx(k);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm01(0.0, 1.0);
//
//  const double pi0 = std::max(pi[0], 1e-300);
//  const double pi1 = std::max(pi[1], 1e-300);
//  const double vei_safe = std::max(vei, 1e-300);
//
//  const double score = r(k) + xk2 * b_block(k);
//  const double denom = std::max(vei_safe + xk2 * vb, 1e-300);
//
//  const double logBF =
//   0.5 * std::log(vei_safe / denom) +
//   0.5 * score * score * vb / (vei_safe * denom);
//
//  const double logp1 = std::log(pi1) + logBF;
//  const double logp0 = std::log(pi0);
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
//  const int dk = (runif(gen) < p1) ? 1 : 0;
//
//  double b_new = 0.0;
//
//  if (dk == 1) {
//   const double lhs = xk2 + vei_safe / vb;
//   const double mean = score / lhs;
//   const double sd = std::sqrt(vei_safe / lhs);
//
//   b_new = mean + sd * norm01(gen);
//  }
//
//  const double diff = b_new - b_block(k);
//
//  if (diff != 0.0) {
//   const int p = static_cast<int>(r.n_elem);
//
//   for (int l = 0; l < p; ++l) {
//    r(l) -= static_cast<double>(XtX(l, k)) * diff;
//   }
//  }
//
//  b_block(k) = b_new;
//  d_block(k) = dk;
// }
//
// // -----------------------------------------------------------------------------
// // Variance and pi updates for individual-level residuals
// // -----------------------------------------------------------------------------
//
// inline void sampleB_ST_individual(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   double ssb_prior,
//   std::mt19937& gen
// ) {
//  double ssb = 0.0;
//  double dfb = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   if (d(i) > 0) {
//    ssb += b(i) * b(i);
//    dfb += 1.0;
//   }
//  }
//
//  const double scale = ssb + nub * ssb_prior;
//
//  std::chi_squared_distribution<double> rchisq(dfb + nub);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  vb = std::max(scale / chi2, 1e-12);
// }
//
// inline void sampleE_ST_individual(
//   double nue,
//   double& ve,
//   const arma::vec& e,
//   double sse_prior,
//   std::mt19937& gen
// ) {
//  const double sse = arma::dot(e, e);
//  const double scale = sse + nue * sse_prior;
//
//  if (!std::isfinite(scale) || scale <= 0.0) {
//   throw std::runtime_error("sampleE_ST_individual: invalid residual scale.");
//  }
//
//  std::chi_squared_distribution<double> rchisq(static_cast<double>(e.n_elem) + nue);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  const double ve_new = scale / chi2;
//
//  if (!std::isfinite(ve_new) || ve_new <= 0.0) {
//   throw std::runtime_error("sampleE_ST_individual: sampled ve is invalid.");
//  }
//
//  ve = std::max(ve_new, 1e-12);
// }
//
// inline double computeG_ST_individual(
//   const arma::mat& X,
//   const arma::rowvec& b
// ) {
//  arma::vec g = X * b.t();
//
//  return arma::dot(g, g) / static_cast<double>(X.n_rows);
// }
//
// inline void samplePi_ST_individual(
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
// inline arma::uvec make_order_within_block(
//   const arma::vec& local_wy,
//   const arma::vec& xx
// ) {
//  const int p = static_cast<int>(local_wy.n_elem);
//  std::vector<int> ord(p);
//
//  std::iota(ord.begin(), ord.end(), 0);
//
//  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
//   const double bha = local_wy(a) / xx(a);
//   const double bhb = local_wy(b) / xx(b);
//
//   return bha * bha > bhb * bhb;
//  });
//
//  arma::uvec out(p);
//
//  for (int k = 0; k < p; ++k) {
//   out(k) = static_cast<arma::uword>(ord[k]);
//  }
//
//  return out;
// }
//
// // -----------------------------------------------------------------------------
// // Main exported function: individual-level block STBLR using dense X
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_individual_blocks(
//   Rcpp::List X_list,
//   Rcpp::List y_list,
//   std::vector<std::vector<double>> b_init,
//   std::vector<int> sets,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   int nit,
//   int nburn,
//   int nthin,
//   int block_nit,
//   int rebuild_every,
//   int ncores,
//   int seed
// ) {
//  const int nt = static_cast<int>(X_list.size());
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: X_list must contain at least one trait.");
//  }
//
//  if (static_cast<int>(y_list.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: y_list length must equal X_list length.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: nthin must be positive.");
//  }
//
//  if (block_nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: block_nit must be positive.");
//  }
//
//  if (rebuild_every < 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: rebuild_every must be >= 0. Use 0 to disable residual rebuilding.");
//  }
//
//  if (static_cast<int>(b_init.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: b_init length must equal number of traits.");
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: pi must have length 2, c(pi0, pi1).");
//  }
//
//  if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: B must be nt x nt.");
//  }
//
//  if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: E must be nt x nt.");
//  }
//
//  arma::mat X0 = Rcpp::as<arma::mat>(X_list[0]);
//  const int m = static_cast<int>(X0.n_cols);
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: ncol(X) must be positive.");
//  }
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: sets must have length ncol(X).");
//  }
//
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   arma::mat Xt = Rcpp::as<arma::mat>(X_list[t]);
//   arma::vec yt = Rcpp::as<arma::vec>(y_list[t]);
//
//   if (static_cast<int>(Xt.n_cols) != m) {
//    throw std::runtime_error("stblr_cpg_omp_individual_blocks: all X matrices must have same number of columns.");
//   }
//
//   if (static_cast<int>(yt.n_elem) != static_cast<int>(Xt.n_rows)) {
//    throw std::runtime_error("stblr_cpg_omp_individual_blocks: y length must equal nrow(X) for each trait.");
//   }
//
//   if (static_cast<int>(b_init[t].size()) != m) {
//    throw std::runtime_error("stblr_cpg_omp_individual_blocks: each b_init[t] must have length m.");
//   }
//
//   if (static_cast<int>(ssb_prior[t].size()) != nt || static_cast<int>(sse_prior[t].size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_individual_blocks: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//  }
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
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
//  arma::vec yy_vec(nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    b_mat(t, i) = b_init[t][i];
//   }
//  }
//
// #ifdef _OPENMP
//  omp_set_dynamic(0);
//  omp_set_num_threads(ncores);
//
//  Rcpp::Rcout
//  << "STBLR individual-level OpenMP max threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
//  std::vector<int> failed(nt, 0);
//  std::vector<std::string> errors(nt);
//
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
//   try {
// #ifdef _OPENMP
// #pragma omp critical
// {
//  Rcpp::Rcout
//  << "trait " << t
//  << " running on thread "
//  << omp_get_thread_num()
//  << " of "
//  << omp_get_num_threads()
//  << "\n";
// }
// #endif
//
//    std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//    arma::mat X_t = Rcpp::as<arma::mat>(X_list[t]);
//    arma::vec y_t = Rcpp::as<arma::vec>(y_list[t]);
//
//    yy_vec(t) = arma::dot(y_t, y_t);
//
//    arma::rowvec b_t = b_mat.row(t);
//    arma::Row<int> d_t(m, arma::fill::zeros);
//
//    arma::vec e_t = y_t - compute_xb_global(X_t, b_t);
//
//    double vb_t = B(t, t);
//    double ve_t = E(t, t);
//    double vg_t = computeG_ST_individual(X_t, b_t);
//    double vei_t = ve_t + adjE * vg_t;
//
//    std::vector<double> pi_t = pi;
//
//    if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//        pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//     throw std::runtime_error("invalid initial pi.");
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
//    std::vector<STBlock> blocks = build_st_blocks_individual(X_t, y_t, sets);
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//     for (const STBlock& blk : blocks) {
//      const int p = static_cast<int>(blk.markers.size());
//
//      arma::vec b_old(p, arma::fill::zeros);
//      arma::vec b_block(p, arma::fill::zeros);
//      arma::Row<int> d_block(p, arma::fill::zeros);
//
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[k];
//
//       b_old(k) = b_t(j);
//       b_block(k) = b_t(j);
//       d_block(k) = d_t(j);
//      }
//
//      arma::vec xb_old = block_xb_from_global_b(X_t, blk.markers, b_t);
//      arma::vec e_without_block = e_t + xb_old;
//
//      arma::vec local_wy(p, arma::fill::zeros);
//
//      for (int k = 0; k < p; ++k) {
//       local_wy(k) = arma::dot(X_t.col(blk.markers[k]), e_without_block);
//      }
//
//      arma::vec r_block = local_wy - xtx_fmat_times_vec(blk.XtX, b_block);
//
//      arma::uvec order = make_order_within_block(local_wy, blk.xx);
//
//      for (int bit = 0; bit < block_nit; ++bit) {
//       for (arma::uword kk = 0; kk < order.n_elem; ++kk) {
//        const int k = static_cast<int>(order(kk));
//
//        sampleBetaC_ST_block_dense(
//         k,
//         pi_t,
//         vb_t,
//         vei_t,
//         blk.xx,
//         blk.XtX,
//         r_block,
//         b_block,
//         d_block,
//         gen_t
//        );
//       }
//      }
//
//      arma::vec db = b_block - b_old;
//
//      bool changed = false;
//
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[k];
//
//       if (db(k) != 0.0) changed = true;
//
//       b_t(j) = b_block(k);
//       d_t(j) = d_block(k);
//      }
//
//      if (changed) {
//       for (int k = 0; k < p; ++k) {
//        if (db(k) != 0.0) {
//         e_t -= X_t.col(blk.markers[k]) * db(k);
//        }
//       }
//      }
//     }
//
//     if (updateB) {
//      sampleB_ST_individual(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       ssb_prior_mat(t, t),
//       gen_t
//      );
//
//      if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//       throw std::runtime_error("vb became invalid after sampleB.");
//      }
//     }
//
//     if (updateE) {
//      if (rebuild_every > 0 && ((it + 1) % rebuild_every == 0)) {
//       e_t = y_t - compute_xb_global(X_t, b_t);
//      }
//
//      sampleE_ST_individual(
//       nue,
//       ve_t,
//       e_t,
//       sse_prior_mat(t, t),
//       gen_t
//      );
//
//      if (!std::isfinite(ve_t) || ve_t <= 0.0) {
//       throw std::runtime_error("ve became invalid after sampleE.");
//      }
//     }
//
//     if (updatePi) {
//      samplePi_ST_individual(d_t, pi_t, gen_t);
//
//      if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//          pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//       throw std::runtime_error("pi became invalid after samplePi.");
//      }
//     }
//
//     vg_t = computeG_ST_individual(X_t, b_t);
//     vei_t = ve_t + adjE * vg_t;
//
//     if (!std::isfinite(vg_t)) {
//      throw std::runtime_error("vg became NaN/Inf after computeG.");
//     }
//
//     vbs_mat(t, it) = vb_t;
//     ves_mat(t, it) = ve_t;
//     vgs_mat(t, it) = vg_t;
//     pis_mat(t, it) = pi_t[1];
//
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       bm_mat(t, i) += b_t(i);
//       dm_mat(t, i) += static_cast<double>(d_t(i));
//      }
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    for (int i = 0; i < m; ++i) {
//     bm_mat(t, i) /= nsamples_t;
//     dm_mat(t, i) /= nsamples_t;
//    }
//
//    if (!bm_mat.row(t).is_finite()) {
//     throw std::runtime_error("posterior mean bm contains NaN/Inf.");
//    }
//
//    if (!dm_mat.row(t).is_finite()) {
//     throw std::runtime_error("posterior mean dm contains NaN/Inf.");
//    }
//
//    b_mat.row(t) = b_t;
//    d_mat.row(t) = d_t;
//
//    final_vb(t) = vb_t;
//    final_ve(t) = ve_t;
//    final_vg(t) = vg_t;
//    final_pi(t) = pi_t[1];
//    nsamples_vec(t) = nsamples_t;
//
//   } catch (const std::exception& e) {
//    failed[t] = 1;
//    errors[t] = e.what();
//   } catch (...) {
//    failed[t] = 1;
//    errors[t] = "unknown error";
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (failed[t]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_individual_blocks failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[t]
//    );
//   }
//  }
//
//  std::vector<std::vector<std::vector<double>>> result(20);
//
//  for (int k = 0; k < 20; ++k) result[k].resize(nt);
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
//   result[16][t].resize(2);
//   result[17][t].resize(2);
//   result[18][t].resize(4);
//   result[19][t].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   arma::mat X_t = Rcpp::as<arma::mat>(X_list[t]);
//   arma::vec y_t = Rcpp::as<arma::vec>(y_list[t]);
//   arma::rowvec b_t = b_mat.row(t);
//   arma::vec e_t = y_t - compute_xb_global(X_t, b_t);
//   arma::vec wy_t = X_t.t() * y_t;
//   arma::vec r_t  = X_t.t() * e_t;
//
//   for (int i = 0; i < m; ++i) {
//    result[0][t][i] = bm_mat(t, i);
//    result[1][t][i] = dm_mat(t, i);
//    result[2][t][i] = wy_t(i);
//    result[3][t][i] = r_t(i);
//    result[4][t][i] = b_mat(t, i);
//    result[5][t][i] = static_cast<double>(d_mat(t, i));
//    result[6][t][i] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int it = 0; it < nit + nburn; ++it) {
//    result[7][t][it] = vbs_mat(t, it);
//    result[8][t][it] = vgs_mat(t, it);
//    result[9][t][it] = ves_mat(t, it);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    result[10][t1][t2] = 0.0;
//    result[11][t1][t2] = 0.0;
//    result[12][t1][t2] = 0.0;
//    result[13][t1][t2] = 0.0;
//    result[14][t1][t2] = 0.0;
//    result[15][t1][t2] = 0.0;
//   }
//
//   result[10][t1][t1] = final_vb(t1);
//   result[11][t1][t1] = final_vg(t1);
//   result[12][t1][t1] = final_ve(t1);
//
//   result[13][t1][t1] = final_vb(t1);
//   result[14][t1][t1] = final_vg(t1);
//   result[15][t1][t1] = final_ve(t1);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[16][t][0] = 1.0 - final_pi(t);
//   result[16][t][1] = final_pi(t);
//
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(t, it);
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi(t);
//
//   result[17][t][0] = 1.0 - mean_pi;
//   result[17][t][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) result[18][t][i] = 0.0;
//   for (int i = 0; i < 2; ++i) result[19][t][i] = 0.0;
//  }
//
//  return result;
// }
//
// // =============================================================================
// // BED-backed individual-level block STBLR extension
// // =============================================================================
//
// #include "packed_bed.h"
//
// struct STBedBlock {
//  std::vector<int> markers;  // global marker indices, 0-based
//  arma::fmat XtX;            // float cached block cross-product
//  arma::vec xx;              // double diag(XtX)
// };
//
// static std::vector<std::string> stblr_copy_bed_files(Rcpp::CharacterVector bed_files) {
//  std::vector<std::string> out(bed_files.size());
//
//  for (int i = 0; i < bed_files.size(); ++i) {
//   out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(bed_files[i]);
//  }
//
//  return out;
// }
//
// static std::vector<std::vector<int>> stblr_copy_int_list(Rcpp::List xlist) {
//  std::vector<std::vector<int>> out(xlist.size());
//
//  for (int f = 0; f < xlist.size(); ++f) {
//   Rcpp::IntegerVector x = xlist[f];
//
//   out[static_cast<std::size_t>(f)].resize(x.size());
//
//   for (int i = 0; i < x.size(); ++i) {
//    out[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)] = x[i];
//   }
//  }
//
//  return out;
// }
//
// static std::vector<int> stblr_copy_rows0_or_empty(
//   Rcpp::Nullable<Rcpp::IntegerVector> rows,
//   int n
// ) {
//  std::vector<int> out;
//
//  if (rows.isNotNull()) {
//   Rcpp::IntegerVector r(rows);
//
//   out.resize(r.size());
//
//   for (int i = 0; i < r.size(); ++i) {
//    out[static_cast<std::size_t>(i)] = r[i] - 1;
//   }
//
//   if (static_cast<int>(out.size()) == n) {
//    bool identity = true;
//
//    for (int i = 0; i < n; ++i) {
//     if (out[static_cast<std::size_t>(i)] != i) {
//      identity = false;
//      break;
//     }
//    }
//
//    if (identity) out.clear();
//   }
//  }
//
//  return out;
// }
//
// static std::vector<double> stblr_flatten_af_list_or_empty(
//   Rcpp::Nullable<Rcpp::List> af
// ) {
//  std::vector<double> out;
//
//  if (af.isNotNull()) {
//   Rcpp::List af_list = Rcpp::as<Rcpp::List>(af.get());
//
//   for (int f = 0; f < af_list.size(); ++f) {
//    Rcpp::NumericVector x = af_list[f];
//
//    for (int i = 0; i < x.size(); ++i) {
//     out.push_back(x[i]);
//    }
//   }
//  }
//
//  return out;
// }
//
// static inline void stblr_bed_maps(
//   double p,
//   bool scale,
//   double& m0,
//   double& m1,
//   double& m2,
//   double& m3
// ) {
//  if (scale) {
//   const double denom = std::sqrt(2.0 * p * (1.0 - p));
//
//   if (denom <= 0.0 || !std::isfinite(denom)) {
//    m0 = m1 = m2 = m3 = 0.0;
//   } else {
//    m0 = (2.0 - 2.0 * p) / denom;
//    m1 = 0.0;
//    m2 = (1.0 - 2.0 * p) / denom;
//    m3 = (0.0 - 2.0 * p) / denom;
//   }
//  } else {
//   m0 = 2.0;
//   m1 = 2.0 * p;
//   m2 = 1.0;
//   m3 = 0.0;
//  }
// }
//
// static inline double stblr_decode_2bit_value(
//   unsigned int code,
//   double m0,
//   double m1,
//   double m2,
//   double m3
// ) {
//  switch (code) {
//  case 0u: return m0;
//  case 1u: return m1;
//  case 2u: return m2;
//  default: return m3;
//  }
// }
//
// static arma::mat stblr_decode_bed_block_double(
//   const PackedBedMatrix& G,
//   const std::vector<int>& markers,
//   const std::vector<double>& af,
//   bool scale
// ) {
//  const int n = G.n;
//  const int p = static_cast<int>(markers.size());
//  const std::size_t nbytes = G.nbytes;
//
//  arma::mat X(n, p, arma::fill::zeros);
//
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
//  for (int kk = 0; kk < p; ++kk) {
//   const int global_i = markers[static_cast<std::size_t>(kk)];
//   const uint8_t* packed = G.row(global_i);
//
//   double m0, m1, m2, m3;
//
//   stblr_bed_maps(
//    af[static_cast<std::size_t>(global_i)],
//      scale,
//      m0,
//      m1,
//      m2,
//      m3
//   );
//
//   for (std::size_t kb = 0; kb < nbytes; ++kb) {
//    unsigned char x = packed[kb];
//    const int jbase = static_cast<int>(kb << 2);
//
//    if (jbase + 0 < n) X(jbase + 0, kk) = stblr_decode_2bit_value((x >> 0) & 3u, m0, m1, m2, m3);
//    if (jbase + 1 < n) X(jbase + 1, kk) = stblr_decode_2bit_value((x >> 2) & 3u, m0, m1, m2, m3);
//    if (jbase + 2 < n) X(jbase + 2, kk) = stblr_decode_2bit_value((x >> 4) & 3u, m0, m1, m2, m3);
//    if (jbase + 3 < n) X(jbase + 3, kk) = stblr_decode_2bit_value((x >> 6) & 3u, m0, m1, m2, m3);
//   }
//  }
//
//  return X;
// }
//
// static std::vector<STBedBlock> stblr_build_bed_blocks(
//   const PackedBedMatrix& G,
//   const std::vector<int>& sets,
//   const std::vector<double>& af,
//   bool scale
// ) {
//  const int m = G.m;
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_build_bed_blocks: sets must have length equal to number of markers.");
//  }
//
//  std::vector<int> labels;
//  labels.reserve(m);
//  std::unordered_map<int, int> label_to_block;
//
//  for (int j = 0; j < m; ++j) {
//   const int lab = sets[static_cast<std::size_t>(j)];
//
//   if (label_to_block.find(lab) == label_to_block.end()) {
//    const int block_id = static_cast<int>(labels.size());
//    label_to_block[lab] = block_id;
//    labels.push_back(lab);
//   }
//  }
//
//  std::vector<STBedBlock> blocks(labels.size());
//
//  for (int j = 0; j < m; ++j) {
//   blocks[static_cast<std::size_t>(label_to_block[sets[static_cast<std::size_t>(j)]])].markers.push_back(j);
//  }
//
//  for (STBedBlock& blk : blocks) {
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af, scale);
//   arma::mat XtX_d = Xb.t() * Xb;
//
//   blk.XtX = arma::conv_to<arma::fmat>::from(XtX_d);
//   blk.xx  = arma::diagvec(XtX_d);
//
//   for (arma::uword k = 0; k < blk.xx.n_elem; ++k) {
//    if (!std::isfinite(blk.xx(k)) || blk.xx(k) <= 0.0) {
//     throw std::runtime_error("stblr_build_bed_blocks: invalid X'X diagonal. Check allele frequencies and monomorphic markers.");
//    }
//   }
//  }
//
//  return blocks;
// }
//
// static arma::vec stblr_bed_xb_global(
//   const PackedBedMatrix& G,
//   const std::vector<STBedBlock>& blocks,
//   const std::vector<double>& af,
//   bool scale,
//   const arma::rowvec& b
// ) {
//  arma::vec xb(G.n, arma::fill::zeros);
//
//  for (const STBedBlock& blk : blocks) {
//   arma::vec b_block(blk.markers.size(), arma::fill::zeros);
//   bool any = false;
//
//   for (std::size_t k = 0; k < blk.markers.size(); ++k) {
//    const int j = blk.markers[k];
//
//    b_block(k) = b(j);
//
//    if (b_block(k) != 0.0) any = true;
//   }
//
//   if (!any) continue;
//
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af, scale);
//
//   xb += Xb * b_block;
//  }
//
//  return xb;
// }
//
// static arma::vec stblr_bed_xty(
//   const PackedBedMatrix& G,
//   const std::vector<STBedBlock>& blocks,
//   const std::vector<double>& af,
//   bool scale,
//   const arma::vec& y
// ) {
//  arma::vec xty(G.m, arma::fill::zeros);
//
//  for (const STBedBlock& blk : blocks) {
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af, scale);
//   arma::vec xty_block = Xb.t() * y;
//
//   for (std::size_t k = 0; k < blk.markers.size(); ++k) {
//    xty(blk.markers[k]) = xty_block(k);
//   }
//  }
//
//  return xty;
// }
//
// static arma::vec stblr_bed_xte(
//   const PackedBedMatrix& G,
//   const std::vector<STBedBlock>& blocks,
//   const std::vector<double>& af,
//   bool scale,
//   const arma::vec& e
// ) {
//  arma::vec xte(G.m, arma::fill::zeros);
//
//  for (const STBedBlock& blk : blocks) {
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af, scale);
//   arma::vec xte_block = Xb.t() * e;
//
//   for (std::size_t k = 0; k < blk.markers.size(); ++k) {
//    xte(blk.markers[k]) = xte_block(k);
//   }
//  }
//
//  return xte;
// }
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_bed_blocks(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   Rcpp::NumericMatrix y,
//   std::vector<std::vector<double>> b_init,
//   std::vector<int> sets,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows,
//   Rcpp::Nullable<Rcpp::List> af,
//   bool scale,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   int nit,
//   int nburn,
//   int nthin,
//   int block_nit,
//   int rebuild_every,
//   int ncores,
//   int seed
// ) {
//  if (nit <= 0 || nburn < 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: nit must be positive and nburn non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: nthin must be positive.");
//  }
//
//  if (block_nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: block_nit must be positive.");
//  }
//
//  if (rebuild_every < 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: rebuild_every must be >= 0. Use 0 to disable residual rebuilding.");
//  }
//
//  std::vector<std::string> bed_files_cpp = stblr_copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = stblr_copy_int_list(cls);
//  std::vector<int> rows0 = stblr_copy_rows0_or_empty(rows, n);
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
//  const int n_used = G.n;
//  const int m = G.m;
//  const int nt = y.ncol();
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: y must have at least one trait column.");
//  }
//
//  if (y.nrow() != n_used) {
//   throw std::runtime_error(
//     "stblr_cpg_omp_bed_blocks: y rows must equal the number of samples used after rows filtering. "
//     "If rows is supplied, pass y already subset to those rows."
//   );
//  }
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: sets must have length equal to the total number of BED markers used.");
//  }
//
//  if (static_cast<int>(b_init.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: b_init must have length nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(b_init[static_cast<std::size_t>(t)].size()) != m) {
//    throw std::runtime_error("stblr_cpg_omp_bed_blocks: each b_init[t] must have length m.");
//   }
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: pi must have length 2.");
//  }
//
//  if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: B must be nt x nt.");
//  }
//
//  if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: E must be nt x nt.");
//  }
//
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  if (static_cast<int>(ssb_prior.size()) != nt || static_cast<int>(sse_prior.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: prior lists must have length nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(ssb_prior[t].size()) != nt || static_cast<int>(sse_prior[t].size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_bed_blocks: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//  }
//
//  std::vector<double> af_cpp = stblr_flatten_af_list_or_empty(af);
//  const bool af_computed = af_cpp.empty();
//
//  if (af_computed) {
//   af_cpp = compute_af_from_packed(G);
//  }
//
//  if (static_cast<int>(af_cpp.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: af must have one value per marker after flattening.");
//  }
//
// #ifdef _OPENMP
//  if (ncores > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(ncores);
//  }
// #endif
//
//  Rcpp::Rcout << "Building BED-backed STBLR blocks: n=" << n_used
//              << ", m=" << m
//              << ", nt=" << nt
//              << ", scale=" << scale
//              << ", af_computed=" << af_computed
//              << "\n";
//
//  std::vector<STBedBlock> blocks = stblr_build_bed_blocks(G, sets, af_cpp, scale);
//
//  Rcpp::Rcout << "Number of STBLR blocks = " << blocks.size() << "\n";
//  Rcpp::Rcout << "BED sampler: safe block-first loop, one decode per block visit.\n";
//
//  arma::mat y_mat(n_used, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int j = 0; j < n_used; ++j) {
//    y_mat(j, t) = y(j, t);
//   }
//  }
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
//
//  arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
//
//  arma::vec vb_vec(nt, arma::fill::zeros);
//  arma::vec ve_vec(nt, arma::fill::zeros);
//  arma::vec vg_vec(nt, arma::fill::zeros);
//  arma::vec vei_vec(nt, arma::fill::zeros);
//  arma::vec final_pi(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  std::vector<std::vector<double>> pi_traits(nt);
//  std::vector<std::mt19937> gens;
//  gens.reserve(nt);
//
//  bool any_initial_b = false;
//
//  for (int t = 0; t < nt; ++t) {
//   vb_vec(t) = B(t, t);
//   ve_vec(t) = E(t, t);
//
//   pi_traits[t] = pi;
//
//   const double psum = pi_traits[t][0] + pi_traits[t][1];
//
//   if (!std::isfinite(psum) || psum <= 0.0 ||
//       pi_traits[t][0] <= 0.0 || pi_traits[t][1] <= 0.0) {
//    throw std::runtime_error("stblr_cpg_omp_bed_blocks: invalid initial pi.");
//   }
//
//   pi_traits[t][0] /= psum;
//   pi_traits[t][1] /= psum;
//
//   gens.emplace_back(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//   for (int i = 0; i < m; ++i) {
//    b_mat(t, i) = b_init[t][i];
//
//    if (b_mat(t, i) != 0.0) {
//     any_initial_b = true;
//    }
//   }
//  }
//
//  arma::mat e_mat(n_used, nt, arma::fill::zeros);
//
//  if (!any_initial_b) {
//   e_mat = y_mat;
//  } else {
//   arma::mat g_mat(n_used, nt, arma::fill::zeros);
//
//   for (const STBedBlock& blk : blocks) {
//    const int p = static_cast<int>(blk.markers.size());
//    arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af_cpp, scale);
//
//    for (int t = 0; t < nt; ++t) {
//     arma::vec b_block(p, arma::fill::zeros);
//     bool any = false;
//
//     for (int k = 0; k < p; ++k) {
//      const int j = blk.markers[static_cast<std::size_t>(k)];
//      b_block(k) = b_mat(t, j);
//
//      if (b_block(k) != 0.0) {
//       any = true;
//      }
//     }
//
//     if (any) {
//      arma::vec add = Xb * b_block;
//      g_mat.col(t) = g_mat.col(t) + add;
//     }
//    }
//   }
//
//   e_mat = y_mat - g_mat;
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   arma::vec g_t = y_mat.col(t) - e_mat.col(t);
//   vg_vec(t) = arma::dot(g_t, g_t) / static_cast<double>(n_used);
//   vei_vec(t) = ve_vec(t) + adjE * vg_vec(t);
//  }
//
//  for (int it = 0; it < nit + nburn; ++it) {
//
//   for (const STBedBlock& blk : blocks) {
//    const int p = static_cast<int>(blk.markers.size());
//
//    arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af_cpp, scale);
//
//    for (int t = 0; t < nt; ++t) {
//     arma::vec b_old(p, arma::fill::zeros);
//     arma::vec b_block(p, arma::fill::zeros);
//     arma::Row<int> d_block(p, arma::fill::zeros);
//
//     for (int k = 0; k < p; ++k) {
//      const int j = blk.markers[static_cast<std::size_t>(k)];
//
//      b_old(k) = b_mat(t, j);
//      b_block(k) = b_mat(t, j);
//      d_block(k) = d_mat(t, j);
//     }
//
//     arma::vec xb_old = Xb * b_old;
//     arma::vec e_without_block = e_mat.col(t) + xb_old;
//
//     arma::vec local_wy = Xb.t() * e_without_block;
//     arma::vec r_block = local_wy - xtx_fmat_times_vec(blk.XtX, b_block);
//
//     arma::uvec order = make_order_within_block(local_wy, blk.xx);
//
//     for (int bit = 0; bit < block_nit; ++bit) {
//      for (arma::uword kk = 0; kk < order.n_elem; ++kk) {
//       const int k = static_cast<int>(order(kk));
//
//       sampleBetaC_ST_block_dense(
//        k,
//        pi_traits[t],
//                 vb_vec(t),
//                 vei_vec(t),
//                 blk.xx,
//                 blk.XtX,
//                 r_block,
//                 b_block,
//                 d_block,
//                 gens[t]
//       );
//      }
//     }
//
//     arma::vec db = b_block - b_old;
//
//     bool changed = false;
//
//     for (int k = 0; k < p; ++k) {
//      const int j = blk.markers[static_cast<std::size_t>(k)];
//
//      if (db(k) != 0.0) {
//       changed = true;
//      }
//
//      b_mat(t, j) = b_block(k);
//      d_mat(t, j) = d_block(k);
//     }
//
//     //if (changed) {
//     // arma::vec delta = Xb * db;
//     // e_mat.col(t) = e_mat.col(t) - delta;
//     //}
//     if (changed) {
//      for (int k = 0; k < p; ++k) {
//       const double diff = db(k);
//
//       if (diff != 0.0) {
//        arma::vec delta = Xb.col(k) * diff;
//        e_mat.col(t) = e_mat.col(t) - delta;
//       }
//      }
//     }
//    }
//   }
//
//   if (updateE && rebuild_every > 0 && ((it + 1) % rebuild_every == 0)) {
//    arma::mat g_mat(n_used, nt, arma::fill::zeros);
//
//    for (const STBedBlock& blk : blocks) {
//     const int p = static_cast<int>(blk.markers.size());
//     arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af_cpp, scale);
//
//     for (int t = 0; t < nt; ++t) {
//      arma::vec b_block(p, arma::fill::zeros);
//      bool any = false;
//
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[static_cast<std::size_t>(k)];
//       b_block(k) = b_mat(t, j);
//
//       if (b_block(k) != 0.0) {
//        any = true;
//       }
//      }
//
//      if (any) {
//       arma::vec add = Xb * b_block;
//       g_mat.col(t) = g_mat.col(t) + add;
//      }
//     }
//    }
//
//    e_mat = y_mat - g_mat;
//   }
//
//   for (int t = 0; t < nt; ++t) {
//    if (updateB) {
//     arma::rowvec b_t = b_mat.row(t);
//     arma::Row<int> d_t = d_mat.row(t);
//
//     sampleB_ST_individual(
//      m,
//      nub,
//      vb_vec(t),
//      b_t,
//      d_t,
//      ssb_prior_mat(t, t),
//      gens[t]
//     );
//    }
//
//    if (updateE) {
//     arma::vec e_t = e_mat.col(t);
//
//     sampleE_ST_individual(
//      nue,
//      ve_vec(t),
//      e_t,
//      sse_prior_mat(t, t),
//      gens[t]
//     );
//    }
//
//    if (updatePi) {
//     arma::Row<int> d_t = d_mat.row(t);
//
//     samplePi_ST_individual(
//      d_t,
//      pi_traits[t],
//               gens[t]
//     );
//    }
//
//    arma::vec g_t = y_mat.col(t) - e_mat.col(t);
//
//    vg_vec(t) = arma::dot(g_t, g_t) / static_cast<double>(n_used);
//    vei_vec(t) = ve_vec(t) + adjE * vg_vec(t);
//
//    if (!std::isfinite(vb_vec(t)) || vb_vec(t) <= 0.0) {
//     throw std::runtime_error("stblr_cpg_omp_bed_blocks: invalid vb.");
//    }
//
//    if (!std::isfinite(ve_vec(t)) || ve_vec(t) <= 0.0) {
//     throw std::runtime_error("stblr_cpg_omp_bed_blocks: invalid ve.");
//    }
//
//    if (!std::isfinite(vg_vec(t))) {
//     throw std::runtime_error("stblr_cpg_omp_bed_blocks: invalid vg.");
//    }
//
//    if (!std::isfinite(pi_traits[t][0]) || !std::isfinite(pi_traits[t][1]) ||
//        pi_traits[t][0] <= 0.0 || pi_traits[t][1] <= 0.0) {
//     throw std::runtime_error("stblr_cpg_omp_bed_blocks: invalid pi.");
//    }
//
//    vbs_mat(t, it) = vb_vec(t);
//    ves_mat(t, it) = ve_vec(t);
//    vgs_mat(t, it) = vg_vec(t);
//    pis_mat(t, it) = pi_traits[t][1];
//
//    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//     nsamples_vec(t) += 1.0;
//
//     for (int i = 0; i < m; ++i) {
//      bm_mat(t, i) += b_mat(t, i);
//      dm_mat(t, i) += static_cast<double>(d_mat(t, i));
//     }
//    }
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (nsamples_vec(t) <= 0.0) {
//    nsamples_vec(t) = 1.0;
//   }
//
//   for (int i = 0; i < m; ++i) {
//    bm_mat(t, i) /= nsamples_vec(t);
//    dm_mat(t, i) /= nsamples_vec(t);
//   }
//
//   final_pi(t) = pi_traits[t][1];
//  }
//
//  arma::mat wy_mat(m, nt, arma::fill::zeros);
//  arma::mat r_mat(m, nt, arma::fill::zeros);
//
//  for (const STBedBlock& blk : blocks) {
//   const int p = static_cast<int>(blk.markers.size());
//
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af_cpp, scale);
//
//   arma::mat wy_block = Xb.t() * y_mat;
//   arma::mat r_block  = Xb.t() * e_mat;
//
//   for (int k = 0; k < p; ++k) {
//    const int j = blk.markers[static_cast<std::size_t>(k)];
//
//    for (int t = 0; t < nt; ++t) {
//     wy_mat(j, t) = wy_block(k, t);
//     r_mat(j, t)  = r_block(k, t);
//    }
//   }
//  }
//
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
//   result[16][t].resize(2);
//   result[17][t].resize(2);
//   result[18][t].resize(4);
//   result[19][t].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    result[0][t][i] = bm_mat(t, i);
//    result[1][t][i] = dm_mat(t, i);
//    result[2][t][i] = wy_mat(i, t);
//    result[3][t][i] = r_mat(i, t);
//    result[4][t][i] = b_mat(t, i);
//    result[5][t][i] = static_cast<double>(d_mat(t, i));
//    result[6][t][i] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int it = 0; it < nit + nburn; ++it) {
//    result[7][t][it] = vbs_mat(t, it);
//    result[8][t][it] = vgs_mat(t, it);
//    result[9][t][it] = ves_mat(t, it);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    result[10][t1][t2] = 0.0;
//    result[11][t1][t2] = 0.0;
//    result[12][t1][t2] = 0.0;
//    result[13][t1][t2] = 0.0;
//    result[14][t1][t2] = 0.0;
//    result[15][t1][t2] = 0.0;
//   }
//
//   result[10][t1][t1] = vb_vec(t1);
//   result[11][t1][t1] = vg_vec(t1);
//   result[12][t1][t1] = ve_vec(t1);
//
//   result[13][t1][t1] = vb_vec(t1);
//   result[14][t1][t1] = vg_vec(t1);
//   result[15][t1][t1] = ve_vec(t1);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[16][t][0] = 1.0 - final_pi(t);
//   result[16][t][1] = final_pi(t);
//
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(t, it);
//    ++npi;
//   }
//
//   if (npi > 0) {
//    mean_pi /= static_cast<double>(npi);
//   } else {
//    mean_pi = final_pi(t);
//   }
//
//   result[17][t][0] = 1.0 - mean_pi;
//   result[17][t][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) {
//    result[18][t][i] = 0.0;
//   }
//
//   for (int i = 0; i < 2; ++i) {
//    result[19][t][i] = 0.0;
//   }
//  }
//
//  return result;
// }
//
// // =============================================================================
// // Pure residual BED-backed block STBLR
// // =============================================================================
// // This version does NOT cache block X'X.
// // It only stores block marker indices and xx = diag(X'X).
// //
// // For each marker:
// //   score = x_i' e + xx_i * b_i
// //   sample b_i, d_i
// //   if diff != 0:
// //      e <- e - x_i * diff
// //
// // This is closest to your old individual-level BayesC update.
// // =============================================================================
//
// struct STBedResidualBlock {
//  std::vector<int> markers;  // global marker indices, 0-based
//  arma::vec xx;              // x_i'x_i for markers in block
// };
//
// static std::vector<STBedResidualBlock> stblr_build_bed_residual_blocks(
//   const PackedBedMatrix& G,
//   const std::vector<int>& sets,
//   const std::vector<double>& af,
//   bool scale
// ) {
//  const int m = G.m;
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_build_bed_residual_blocks: sets must have length equal to number of markers.");
//  }
//
//  std::vector<int> labels;
//  labels.reserve(m);
//
//  std::unordered_map<int, int> label_to_block;
//
//  for (int j = 0; j < m; ++j) {
//   const int lab = sets[static_cast<std::size_t>(j)];
//
//   if (label_to_block.find(lab) == label_to_block.end()) {
//    const int block_id = static_cast<int>(labels.size());
//    label_to_block[lab] = block_id;
//    labels.push_back(lab);
//   }
//  }
//
//  std::vector<STBedResidualBlock> blocks(labels.size());
//
//  for (int j = 0; j < m; ++j) {
//   const int block_id = label_to_block[sets[static_cast<std::size_t>(j)]];
//   blocks[static_cast<std::size_t>(block_id)].markers.push_back(j);
//  }
//
//  for (STBedResidualBlock& blk : blocks) {
//   const int p = static_cast<int>(blk.markers.size());
//
//   arma::mat Xb = stblr_decode_bed_block_double(
//    G,
//    blk.markers,
//    af,
//    scale
//   );
//
//   blk.xx.set_size(p);
//
//   for (int k = 0; k < p; ++k) {
//    blk.xx(k) = arma::dot(Xb.col(k), Xb.col(k));
//
//    if (!std::isfinite(blk.xx(k)) || blk.xx(k) <= 0.0) {
//     throw std::runtime_error(
//       "stblr_build_bed_residual_blocks: invalid X'X diagonal. "
//       "Check allele frequencies and monomorphic markers."
//     );
//    }
//   }
//  }
//
//  return blocks;
// }
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_bed_residual_blocks(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   Rcpp::NumericMatrix y,
//   std::vector<std::vector<double>> b_init,
//   std::vector<int> sets,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows,
//   Rcpp::Nullable<Rcpp::List> af,
//   bool scale,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   int nit,
//   int nburn,
//   int nthin,
//   int block_nit,
//   int rebuild_every,
//   int ncores,
//   int seed
// ) {
//  if (nit <= 0 || nburn < 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: nit must be positive and nburn non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: nthin must be positive.");
//  }
//
//  if (block_nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: block_nit must be positive.");
//  }
//
//  if (rebuild_every < 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: rebuild_every must be >= 0. Use 0 to disable residual rebuilding.");
//  }
//
//  std::vector<std::string> bed_files_cpp = stblr_copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = stblr_copy_int_list(cls);
//  std::vector<int> rows0 = stblr_copy_rows0_or_empty(rows, n);
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
//  const int n_used = G.n;
//  const int m = G.m;
//  const int nt = y.ncol();
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: y must have at least one trait column.");
//  }
//
//  if (y.nrow() != n_used) {
//   throw std::runtime_error(
//     "stblr_cpg_omp_bed_residual_blocks: y rows must equal the number of samples used after rows filtering. "
//     "If rows is supplied, pass y already subset to those rows."
//   );
//  }
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: sets must have length equal to the total number of BED markers used.");
//  }
//
//  if (static_cast<int>(b_init.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: b_init must have length nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(b_init[static_cast<std::size_t>(t)].size()) != m) {
//    throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: each b_init[t] must have length m.");
//   }
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: pi must have length 2.");
//  }
//
//  if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: B must be nt x nt.");
//  }
//
//  if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: E must be nt x nt.");
//  }
//
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  if (static_cast<int>(ssb_prior.size()) != nt ||
//      static_cast<int>(sse_prior.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: prior lists must have length nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(ssb_prior[t].size()) != nt ||
//       static_cast<int>(sse_prior[t].size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//  }
//
//  std::vector<double> af_cpp = stblr_flatten_af_list_or_empty(af);
//  const bool af_computed = af_cpp.empty();
//
//  if (af_computed) {
//   af_cpp = compute_af_from_packed(G);
//  }
//
//  if (static_cast<int>(af_cpp.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: af must have one value per marker after flattening.");
//  }
//
// #ifdef _OPENMP
//  if (ncores > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(ncores);
//  }
// #endif
//
//  Rcpp::Rcout << "Building pure residual BED STBLR blocks: n=" << n_used
//              << ", m=" << m
//              << ", nt=" << nt
//              << ", scale=" << scale
//              << ", af_computed=" << af_computed
//              << "\n";
//
//  std::vector<STBedResidualBlock> blocks = stblr_build_bed_residual_blocks(
//   G,
//   sets,
//   af_cpp,
//   scale
//  );
//
//  Rcpp::Rcout << "Number of pure residual STBLR blocks = " << blocks.size() << "\n";
//  Rcpp::Rcout << "BED sampler: pure residual marker-wise updates, no block X'X cache.\n";
//
//  arma::mat y_mat(n_used, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int j = 0; j < n_used; ++j) {
//    y_mat(j, t) = y(j, t);
//   }
//  }
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
//
//  arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
//
//  arma::vec vb_vec(nt, arma::fill::zeros);
//  arma::vec ve_vec(nt, arma::fill::zeros);
//  arma::vec vg_vec(nt, arma::fill::zeros);
//  arma::vec vei_vec(nt, arma::fill::zeros);
//  arma::vec final_pi(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  std::vector<std::vector<double>> pi_traits(nt);
//  std::vector<std::mt19937> gens;
//  gens.reserve(nt);
//
//  bool any_initial_b = false;
//
//  for (int t = 0; t < nt; ++t) {
//   vb_vec(t) = B(t, t);
//   ve_vec(t) = E(t, t);
//
//   pi_traits[t] = pi;
//
//   const double psum = pi_traits[t][0] + pi_traits[t][1];
//
//   if (!std::isfinite(psum) || psum <= 0.0 ||
//       pi_traits[t][0] <= 0.0 || pi_traits[t][1] <= 0.0) {
//    throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid initial pi.");
//   }
//
//   pi_traits[t][0] /= psum;
//   pi_traits[t][1] /= psum;
//
//   gens.emplace_back(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//   for (int i = 0; i < m; ++i) {
//    b_mat(t, i) = b_init[t][i];
//
//    if (b_mat(t, i) != 0.0) {
//     any_initial_b = true;
//    }
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Initial residuals e = y - Xb.
//  // --------------------------------------------------------------------------
//
//  arma::mat e_mat(n_used, nt, arma::fill::zeros);
//
//  if (!any_initial_b) {
//   e_mat = y_mat;
//  } else {
//   arma::mat g_mat(n_used, nt, arma::fill::zeros);
//
//   for (const STBedResidualBlock& blk : blocks) {
//    const int p = static_cast<int>(blk.markers.size());
//
//    arma::mat Xb = stblr_decode_bed_block_double(
//     G,
//     blk.markers,
//     af_cpp,
//     scale
//    );
//
//    for (int t = 0; t < nt; ++t) {
//     arma::vec b_block(p, arma::fill::zeros);
//     bool any = false;
//
//     for (int k = 0; k < p; ++k) {
//      const int j = blk.markers[static_cast<std::size_t>(k)];
//      b_block(k) = b_mat(t, j);
//
//      if (b_block(k) != 0.0) {
//       any = true;
//      }
//     }
//
//     if (any) {
//      g_mat.col(t) = g_mat.col(t) + Xb * b_block;
//     }
//    }
//   }
//
//   e_mat = y_mat - g_mat;
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   arma::vec g_t = y_mat.col(t) - e_mat.col(t);
//   vg_vec(t) = arma::dot(g_t, g_t) / static_cast<double>(n_used);
//   vei_vec(t) = ve_vec(t) + adjE * vg_vec(t);
//  }
//
//  // // --------------------------------------------------------------------------
//  // // MCMC.
//  // // --------------------------------------------------------------------------
//  //
//  // for (int it = 0; it < nit + nburn; ++it) {
//  //
//  //  for (const STBedResidualBlock& blk : blocks) {
//  //   const int p = static_cast<int>(blk.markers.size());
//  //
//  //   arma::mat Xb = stblr_decode_bed_block_double(
//  //    G,
//  //    blk.markers,
//  //    af_cpp,
//  //    scale
//  //   );
//  //
//  //   for (int t = 0; t < nt; ++t) {
//  //    std::uniform_real_distribution<double> runif(0.0, 1.0);
//  //    std::normal_distribution<double> norm01(0.0, 1.0);
//  //
//  //    arma::vec e_t = e_mat.col(t);
//  //
//  //    for (int bit = 0; bit < block_nit; ++bit) {
//  //     for (int k = 0; k < p; ++k) {
//  //      const int j = blk.markers[static_cast<std::size_t>(k)];
//  //
//  //      const double xk2 = blk.xx(k);
//  //      const double b_old = b_mat(t, j);
//  //
//  //      const double pi0 = std::max(pi_traits[t][0], 1e-300);
//  //      const double pi1 = std::max(pi_traits[t][1], 1e-300);
//  //      const double vei_safe = std::max(vei_vec(t), 1e-300);
//  //
//  //      // Current residual includes the current marker effect.
//  //      // Therefore add back xx_i * b_i to get x_i' residual_without_i.
//  //      const double score = arma::dot(Xb.col(k), e_t) + xk2 * b_old;
//  //
//  //      const double denom = std::max(vei_safe + xk2 * vb_vec(t), 1e-300);
//  //
//  //      const double logBF =
//  //       0.5 * std::log(vei_safe / denom) +
//  //       0.5 * score * score * vb_vec(t) / (vei_safe * denom);
//  //
//  //      const double logp1 = std::log(pi1) + logBF;
//  //      const double logp0 = std::log(pi0);
//  //      const double delta_log = logp0 - logp1;
//  //
//  //      double p1_incl = 0.0;
//  //
//  //      if (delta_log > 35.0) {
//  //       p1_incl = 0.0;
//  //      } else if (delta_log < -35.0) {
//  //       p1_incl = 1.0;
//  //      } else {
//  //       p1_incl = 1.0 / (1.0 + std::exp(delta_log));
//  //      }
//  //
//  //      const int d_new = (runif(gens[t]) < p1_incl) ? 1 : 0;
//  //
//  //      double b_new = 0.0;
//  //
//  //      if (d_new == 1) {
//  //       const double lhs = xk2 + vei_safe / vb_vec(t);
//  //       const double mean = score / lhs;
//  //       const double sd = std::sqrt(vei_safe / lhs);
//  //
//  //       b_new = mean + sd * norm01(gens[t]);
//  //      }
//  //
//  //      const double diff = b_new - b_old;
//  //
//  //      if (diff != 0.0) {
//  //       // Pure individual-level residual update:
//  //       // e <- e - x_i * diff
//  //       e_t -= Xb.col(k) * diff;
//  //      }
//  //
//  //      b_mat(t, j) = b_new;
//  //      d_mat(t, j) = d_new;
//  //     }
//  //    }
//  //
//  //    e_mat.col(t) = e_t;
//  //   }
//  //  }
//  //
//  //  // --------------------------------------------------------------------------
//  //  // Optional residual rebuild from BED.
//  //  // --------------------------------------------------------------------------
//  //
//  //  if (updateE && rebuild_every > 0 && ((it + 1) % rebuild_every == 0)) {
//  //   arma::mat g_mat(n_used, nt, arma::fill::zeros);
//  //
//  //   for (const STBedResidualBlock& blk : blocks) {
//  //    const int p = static_cast<int>(blk.markers.size());
//  //
//  //    arma::mat Xb = stblr_decode_bed_block_double(
//  //     G,
//  //     blk.markers,
//  //     af_cpp,
//  //     scale
//  //    );
//  //
//  //    for (int t = 0; t < nt; ++t) {
//  //     arma::vec b_block(p, arma::fill::zeros);
//  //     bool any = false;
//  //
//  //     for (int k = 0; k < p; ++k) {
//  //      const int j = blk.markers[static_cast<std::size_t>(k)];
//  //      b_block(k) = b_mat(t, j);
//  //
//  //      if (b_block(k) != 0.0) {
//  //       any = true;
//  //      }
//  //     }
//  //
//  //     if (any) {
//  //      g_mat.col(t) = g_mat.col(t) + Xb * b_block;
//  //     }
//  //    }
//  //   }
//  //
//  //   e_mat = y_mat - g_mat;
//  //  }
//  //
//  //  // --------------------------------------------------------------------------
//  //  // Hyperparameter updates.
//  //  // --------------------------------------------------------------------------
//  //
//  //  for (int t = 0; t < nt; ++t) {
//  //   if (updateB) {
//  //    arma::rowvec b_t = b_mat.row(t);
//  //    arma::Row<int> d_t = d_mat.row(t);
//  //
//  //    sampleB_ST_individual(
//  //     m,
//  //     nub,
//  //     vb_vec(t),
//  //     b_t,
//  //     d_t,
//  //     ssb_prior_mat(t, t),
//  //     gens[t]
//  //    );
//  //   }
//  //
//  //   if (updateE) {
//  //    arma::vec e_t = e_mat.col(t);
//  //
//  //    sampleE_ST_individual(
//  //     nue,
//  //     ve_vec(t),
//  //     e_t,
//  //     sse_prior_mat(t, t),
//  //     gens[t]
//  //    );
//  //   }
//  //
//  //   if (updatePi) {
//  //    arma::Row<int> d_t = d_mat.row(t);
//  //
//  //    samplePi_ST_individual(
//  //     d_t,
//  //     pi_traits[t],
//  //              gens[t]
//  //    );
//  //   }
//  //
//  //   arma::vec g_t = y_mat.col(t) - e_mat.col(t);
//  //
//  //   vg_vec(t) = arma::dot(g_t, g_t) / static_cast<double>(n_used);
//  //   vei_vec(t) = ve_vec(t) + adjE * vg_vec(t);
//  //
//  //   if (!std::isfinite(vb_vec(t)) || vb_vec(t) <= 0.0) {
//  //    throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid vb.");
//  //   }
//  //
//  //   if (!std::isfinite(ve_vec(t)) || ve_vec(t) <= 0.0) {
//  //    throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid ve.");
//  //   }
//  //
//  //   if (!std::isfinite(vg_vec(t))) {
//  //    throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid vg.");
//  //   }
//  //
//  //   if (!std::isfinite(pi_traits[t][0]) || !std::isfinite(pi_traits[t][1]) ||
//  //       pi_traits[t][0] <= 0.0 || pi_traits[t][1] <= 0.0) {
//  //    throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid pi.");
//  //   }
//  //
//  //   vbs_mat(t, it) = vb_vec(t);
//  //   ves_mat(t, it) = ve_vec(t);
//  //   vgs_mat(t, it) = vg_vec(t);
//  //   pis_mat(t, it) = pi_traits[t][1];
//  //
//  //   if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//  //    nsamples_vec(t) += 1.0;
//  //
//  //    for (int i = 0; i < m; ++i) {
//  //     bm_mat(t, i) += b_mat(t, i);
//  //     dm_mat(t, i) += static_cast<double>(d_mat(t, i));
//  //    }
//  //   }
//  //  }
//  // }
//
//  // --------------------------------------------------------------------------
//  // MCMC.
//  // --------------------------------------------------------------------------
//  //
//  // Hybrid residual-score update:
//  //
//  // For each decoded block:
//  //   S = Xb' e_mat              // p x nt, all trait scores at once
//  //
//  // For marker k, trait t:
//  //   score = S(k,t) + xx_k b_k
//  //
//  // If b_k changes by diff:
//  //   e_t      <- e_t - x_k diff
//  //   S(:, t) <- S(:, t) - Xb' x_k diff
//  //
//  // This keeps scores exact within the block without storing block X'X.
//  // Only changed markers require the local Xb' x_k update.
//  // --------------------------------------------------------------------------
//
//  for (int it = 0; it < nit + nburn; ++it) {
//
//   for (const STBedResidualBlock& blk : blocks) {
//    const int p = static_cast<int>(blk.markers.size());
//
//    arma::mat Xb = stblr_decode_bed_block_double(
//     G,
//     blk.markers,
//     af_cpp,
//     scale
//    );
//
//    // Current block scores for all traits:
//    // S(k,t) = x_k' e_t.
//    arma::mat S = Xb.t() * e_mat;
//
//    for (int t = 0; t < nt; ++t) {
//     std::uniform_real_distribution<double> runif(0.0, 1.0);
//     std::normal_distribution<double> norm01(0.0, 1.0);
//
//     arma::vec e_t = e_mat.col(t);
//
//     const double pi0_base = std::max(pi_traits[t][0], 1e-300);
//     const double pi1_base = std::max(pi_traits[t][1], 1e-300);
//
//     for (int bit = 0; bit < block_nit; ++bit) {
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[static_cast<std::size_t>(k)];
//
//       const double xk2 = blk.xx(k);
//       const double b_old = b_mat(t, j);
//
//       const double vei_safe = std::max(vei_vec(t), 1e-300);
//
//       // Current residual e_t includes the current marker effect.
//       // Therefore:
//       //   x_k' residual_without_k = x_k' e_t + x_k'x_k b_k
//       //
//       // S(k,t) is maintained exactly as x_k' e_t.
//       const double score = S(k, t) + xk2 * b_old;
//
//       const double denom = std::max(
//        vei_safe + xk2 * vb_vec(t),
//        1e-300
//       );
//
//       const double logBF =
//        0.5 * std::log(vei_safe / denom) +
//        0.5 * score * score * vb_vec(t) / (vei_safe * denom);
//
//       const double logp1 = std::log(pi1_base) + logBF;
//       const double logp0 = std::log(pi0_base);
//       const double delta_log = logp0 - logp1;
//
//       double p1_incl = 0.0;
//
//       if (delta_log > 35.0) {
//        p1_incl = 0.0;
//       } else if (delta_log < -35.0) {
//        p1_incl = 1.0;
//       } else {
//        p1_incl = 1.0 / (1.0 + std::exp(delta_log));
//       }
//
//       const int d_new = (runif(gens[t]) < p1_incl) ? 1 : 0;
//
//       double b_new = 0.0;
//
//       if (d_new == 1) {
//        const double lhs = xk2 + vei_safe / vb_vec(t);
//        const double mean = score / lhs;
//        const double sd = std::sqrt(vei_safe / lhs);
//
//        b_new = mean + sd * norm01(gens[t]);
//       }
//
//       const double diff = b_new - b_old;
//
//       if (diff != 0.0) {
//        // Update individual residual:
//        //   e_t <- e_t - x_k diff
//        arma::vec xk = Xb.col(k);
//        e_t -= xk * diff;
//
//        // Keep block score vector current:
//        //   S(:,t) <- S(:,t) - Xb' x_k diff
//        //
//        // This computes only the needed local LD column on demand.
//        arma::vec xkx = Xb.t() * xk;
//        S.col(t) = S.col(t) - xkx * diff;
//       }
//
//       b_mat(t, j) = b_new;
//       d_mat(t, j) = d_new;
//      }
//     }
//
//     e_mat.col(t) = e_t;
//    }
//   }
//
//   // --------------------------------------------------------------------------
//   // Optional residual rebuild from BED.
//   // --------------------------------------------------------------------------
//
//   if (updateE && rebuild_every > 0 && ((it + 1) % rebuild_every == 0)) {
//    arma::mat g_mat(n_used, nt, arma::fill::zeros);
//
//    for (const STBedResidualBlock& blk : blocks) {
//     const int p = static_cast<int>(blk.markers.size());
//
//     arma::mat Xb = stblr_decode_bed_block_double(
//      G,
//      blk.markers,
//      af_cpp,
//      scale
//     );
//
//     for (int t = 0; t < nt; ++t) {
//      arma::vec b_block(p, arma::fill::zeros);
//      bool any = false;
//
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[static_cast<std::size_t>(k)];
//       b_block(k) = b_mat(t, j);
//
//       if (b_block(k) != 0.0) {
//        any = true;
//       }
//      }
//
//      if (any) {
//       g_mat.col(t) = g_mat.col(t) + Xb * b_block;
//      }
//     }
//    }
//
//    e_mat = y_mat - g_mat;
//   }
//
//   // --------------------------------------------------------------------------
//   // Hyperparameter updates.
//   // --------------------------------------------------------------------------
//
//   for (int t = 0; t < nt; ++t) {
//    if (updateB) {
//     arma::rowvec b_t = b_mat.row(t);
//     arma::Row<int> d_t = d_mat.row(t);
//
//     sampleB_ST_individual(
//      m,
//      nub,
//      vb_vec(t),
//      b_t,
//      d_t,
//      ssb_prior_mat(t, t),
//      gens[t]
//     );
//    }
//
//    if (updateE) {
//     arma::vec e_t = e_mat.col(t);
//
//     sampleE_ST_individual(
//      nue,
//      ve_vec(t),
//      e_t,
//      sse_prior_mat(t, t),
//      gens[t]
//     );
//    }
//
//    if (updatePi) {
//     arma::Row<int> d_t = d_mat.row(t);
//
//     samplePi_ST_individual(
//      d_t,
//      pi_traits[t],
//               gens[t]
//     );
//    }
//
//    arma::vec g_t = y_mat.col(t) - e_mat.col(t);
//
//    vg_vec(t) = arma::dot(g_t, g_t) / static_cast<double>(n_used);
//    vei_vec(t) = ve_vec(t) + adjE * vg_vec(t);
//
//    if (!std::isfinite(vb_vec(t)) || vb_vec(t) <= 0.0) {
//     throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid vb.");
//    }
//
//    if (!std::isfinite(ve_vec(t)) || ve_vec(t) <= 0.0) {
//     throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid ve.");
//    }
//
//    if (!std::isfinite(vg_vec(t))) {
//     throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid vg.");
//    }
//
//    if (!std::isfinite(pi_traits[t][0]) || !std::isfinite(pi_traits[t][1]) ||
//        pi_traits[t][0] <= 0.0 || pi_traits[t][1] <= 0.0) {
//     throw std::runtime_error("stblr_cpg_omp_bed_residual_blocks: invalid pi.");
//    }
//
//    vbs_mat(t, it) = vb_vec(t);
//    ves_mat(t, it) = ve_vec(t);
//    vgs_mat(t, it) = vg_vec(t);
//    pis_mat(t, it) = pi_traits[t][1];
//
//    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//     nsamples_vec(t) += 1.0;
//
//     for (int i = 0; i < m; ++i) {
//      bm_mat(t, i) += b_mat(t, i);
//      dm_mat(t, i) += static_cast<double>(d_mat(t, i));
//     }
//    }
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (nsamples_vec(t) <= 0.0) {
//    nsamples_vec(t) = 1.0;
//   }
//
//   for (int i = 0; i < m; ++i) {
//    bm_mat(t, i) /= nsamples_vec(t);
//    dm_mat(t, i) /= nsamples_vec(t);
//   }
//
//   final_pi(t) = pi_traits[t][1];
//  }
//
//  // --------------------------------------------------------------------------
//  // Final wy and r. Same output convention as your existing ST sampler.
//  // --------------------------------------------------------------------------
//
//  arma::mat wy_mat(m, nt, arma::fill::zeros);
//  arma::mat r_mat(m, nt, arma::fill::zeros);
//
//  for (const STBedResidualBlock& blk : blocks) {
//   const int p = static_cast<int>(blk.markers.size());
//
//   arma::mat Xb = stblr_decode_bed_block_double(
//    G,
//    blk.markers,
//    af_cpp,
//    scale
//   );
//
//   arma::mat wy_block = Xb.t() * y_mat;
//   arma::mat r_block  = Xb.t() * e_mat;
//
//   for (int k = 0; k < p; ++k) {
//    const int j = blk.markers[static_cast<std::size_t>(k)];
//
//    for (int t = 0; t < nt; ++t) {
//     wy_mat(j, t) = wy_block(k, t);
//     r_mat(j, t)  = r_block(k, t);
//    }
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Build 20-slot result.
//  // --------------------------------------------------------------------------
//
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
//   result[16][t].resize(2);
//   result[17][t].resize(2);
//   result[18][t].resize(4);
//   result[19][t].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    result[0][t][i] = bm_mat(t, i);
//    result[1][t][i] = dm_mat(t, i);
//    result[2][t][i] = wy_mat(i, t);
//    result[3][t][i] = r_mat(i, t);
//    result[4][t][i] = b_mat(t, i);
//    result[5][t][i] = static_cast<double>(d_mat(t, i));
//    result[6][t][i] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int it = 0; it < nit + nburn; ++it) {
//    result[7][t][it] = vbs_mat(t, it);
//    result[8][t][it] = vgs_mat(t, it);
//    result[9][t][it] = ves_mat(t, it);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    result[10][t1][t2] = 0.0;
//    result[11][t1][t2] = 0.0;
//    result[12][t1][t2] = 0.0;
//    result[13][t1][t2] = 0.0;
//    result[14][t1][t2] = 0.0;
//    result[15][t1][t2] = 0.0;
//   }
//
//   result[10][t1][t1] = vb_vec(t1);
//   result[11][t1][t1] = vg_vec(t1);
//   result[12][t1][t1] = ve_vec(t1);
//
//   result[13][t1][t1] = vb_vec(t1);
//   result[14][t1][t1] = vg_vec(t1);
//   result[15][t1][t1] = ve_vec(t1);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[16][t][0] = 1.0 - final_pi(t);
//   result[16][t][1] = final_pi(t);
//
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(t, it);
//    ++npi;
//   }
//
//   if (npi > 0) {
//    mean_pi /= static_cast<double>(npi);
//   } else {
//    mean_pi = final_pi(t);
//   }
//
//   result[17][t][0] = 1.0 - mean_pi;
//   result[17][t][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) {
//    result[18][t][i] = 0.0;
//   }
//
//   for (int i = 0; i < 2; ++i) {
//    result[19][t][i] = 0.0;
//   }
//  }
//
//  return result;
// }
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_bed_blocks(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   Rcpp::NumericMatrix y,
//   std::vector<std::vector<double>> b_init,
//   std::vector<int> sets,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows,
//   Rcpp::Nullable<Rcpp::List> af,
//   bool scale,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   int nit,
//   int nburn,
//   int nthin,
//   int block_nit,
//   int rebuild_every,
//   int ncores,
//   int seed
// ) {
//  if (nit <= 0 || nburn < 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: nit must be positive and nburn non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: nthin must be positive.");
//  }
//
//  if (block_nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: block_nit must be positive.");
//  }
//
//  if (rebuild_every < 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: rebuild_every must be >= 0. Use 0 to disable residual rebuilding.");
//  }
//
//  std::vector<std::string> bed_files_cpp = stblr_copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = stblr_copy_int_list(cls);
//  std::vector<int> rows0 = stblr_copy_rows0_or_empty(rows, n);
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
//  const int n_used = G.n;
//  const int m = G.m;
//  const int nt = y.ncol();
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: y must have at least one trait column.");
//  }
//
//  if (y.nrow() != n_used) {
//   throw std::runtime_error(
//     "stblr_cpg_omp_bed_blocks: y rows must equal the number of samples used after rows filtering. "
//     "If rows is supplied, pass y already subset to those rows."
//   );
//  }
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: sets must have length equal to the total number of BED markers used.");
//  }
//
//  if (static_cast<int>(b_init.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: b_init must have length nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(b_init[static_cast<std::size_t>(t)].size()) != m) {
//    throw std::runtime_error("stblr_cpg_omp_bed_blocks: each b_init[t] must have length m.");
//   }
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: pi must have length 2.");
//  }
//
//  if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: B must be nt x nt.");
//  }
//
//  if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: E must be nt x nt.");
//  }
//
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  if (static_cast<int>(ssb_prior.size()) != nt || static_cast<int>(sse_prior.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: prior lists must have length nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(ssb_prior[t].size()) != nt || static_cast<int>(sse_prior[t].size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_bed_blocks: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//  }
//
//  std::vector<double> af_cpp = stblr_flatten_af_list_or_empty(af);
//  const bool af_computed = af_cpp.empty();
//
//  if (af_computed) af_cpp = compute_af_from_packed(G);
//
//  if (static_cast<int>(af_cpp.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: af must have one value per marker after flattening.");
//  }
//
// #ifdef _OPENMP
//  if (ncores > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(ncores);
//  }
// #endif
//
//  Rcpp::Rcout << "Building BED-backed STBLR blocks: n=" << n_used
//              << ", m=" << m
//              << ", nt=" << nt
//              << ", scale=" << scale
//              << ", af_computed=" << af_computed
//              << "\n";
//
//  std::vector<STBedBlock> blocks = stblr_build_bed_blocks(G, sets, af_cpp, scale);
//
//  Rcpp::Rcout << "Number of STBLR blocks = " << blocks.size() << "\n";
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
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
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    b_mat(t, i) = b_init[t][i];
//   }
//  }
//
//  std::vector<int> failed(nt, 0);
//  std::vector<std::string> errors(nt);
//
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
//   try {
//    std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//    arma::vec y_t(n_used, arma::fill::zeros);
//
//    for (int j = 0; j < n_used; ++j) {
//     y_t(j) = y(j, t);
//    }
//
//    arma::rowvec b_t = b_mat.row(t);
//    arma::Row<int> d_t(m, arma::fill::zeros);
//
//    arma::vec e_t = y_t - stblr_bed_xb_global(G, blocks, af_cpp, scale, b_t);
//
//    double vb_t = B(t, t);
//    double ve_t = E(t, t);
//    double vg_t = arma::dot(y_t - e_t, y_t - e_t) / static_cast<double>(n_used);
//    double vei_t = ve_t + adjE * vg_t;
//
//    std::vector<double> pi_t = pi;
//
//    {
//     const double psum = pi_t[0] + pi_t[1];
//
//     if (!std::isfinite(psum) || psum <= 0.0 || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//      throw std::runtime_error("invalid initial pi.");
//     }
//
//     pi_t[0] /= psum;
//     pi_t[1] /= psum;
//    }
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//     for (const STBedBlock& blk : blocks) {
//      const int p = static_cast<int>(blk.markers.size());
//
//      arma::vec b_old(p, arma::fill::zeros);
//      arma::vec b_block(p, arma::fill::zeros);
//      arma::Row<int> d_block(p, arma::fill::zeros);
//
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[static_cast<std::size_t>(k)];
//
//       b_old(k) = b_t(j);
//       b_block(k) = b_t(j);
//       d_block(k) = d_t(j);
//      }
//
//      arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af_cpp, scale);
//
//      arma::vec xb_old = Xb * b_old;
//      arma::vec e_without_block = e_t + xb_old;
//
//      arma::vec local_wy = Xb.t() * e_without_block;
//      arma::vec r_block = local_wy - xtx_fmat_times_vec(blk.XtX, b_block);
//
//      arma::uvec order = make_order_within_block(local_wy, blk.xx);
//
//      for (int bit = 0; bit < block_nit; ++bit) {
//       for (arma::uword kk = 0; kk < order.n_elem; ++kk) {
//        const int k = static_cast<int>(order(kk));
//
//        sampleBetaC_ST_block_dense(
//         k,
//         pi_t,
//         vb_t,
//         vei_t,
//         blk.xx,
//         blk.XtX,
//         r_block,
//         b_block,
//         d_block,
//         gen_t
//        );
//       }
//      }
//
//      arma::vec db = b_block - b_old;
//
//      bool changed = false;
//
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[static_cast<std::size_t>(k)];
//
//       if (db(k) != 0.0) changed = true;
//
//       b_t(j) = b_block(k);
//       d_t(j) = d_block(k);
//      }
//
//      if (changed) {
//       e_t -= Xb * db;
//      }
//     }
//
//     if (updateB) {
//      sampleB_ST_individual(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       ssb_prior_mat(t, t),
//       gen_t
//      );
//     }
//
//     if (updateE) {
//      if (rebuild_every > 0 && ((it + 1) % rebuild_every == 0)) {
//       e_t = y_t - stblr_bed_xb_global(G, blocks, af_cpp, scale, b_t);
//      }
//
//      sampleE_ST_individual(
//       nue,
//       ve_t,
//       e_t,
//       sse_prior_mat(t, t),
//       gen_t
//      );
//     }
//
//     if (updatePi) {
//      samplePi_ST_individual(d_t, pi_t, gen_t);
//     }
//
//     arma::vec g_t = y_t - e_t;
//     vg_t = arma::dot(g_t, g_t) / static_cast<double>(n_used);
//     vei_t = ve_t + adjE * vg_t;
//
//     if (!std::isfinite(vb_t) || vb_t <= 0.0) throw std::runtime_error("invalid vb.");
//     if (!std::isfinite(ve_t) || ve_t <= 0.0) throw std::runtime_error("invalid ve.");
//     if (!std::isfinite(vg_t)) throw std::runtime_error("invalid vg.");
//
//     if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//      throw std::runtime_error("invalid pi.");
//     }
//
//     vbs_mat(t, it) = vb_t;
//     ves_mat(t, it) = ve_t;
//     vgs_mat(t, it) = vg_t;
//     pis_mat(t, it) = pi_t[1];
//
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       bm_mat(t, i) += b_t(i);
//       dm_mat(t, i) += static_cast<double>(d_t(i));
//      }
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    for (int i = 0; i < m; ++i) {
//     bm_mat(t, i) /= nsamples_t;
//     dm_mat(t, i) /= nsamples_t;
//    }
//
//    b_mat.row(t) = b_t;
//    d_mat.row(t) = d_t;
//
//    final_vb(t) = vb_t;
//    final_ve(t) = ve_t;
//    final_vg(t) = vg_t;
//    final_pi(t) = pi_t[1];
//    nsamples_vec(t) = nsamples_t;
//
//   } catch (const std::exception& e) {
//    failed[t] = 1;
//    errors[t] = e.what();
//   } catch (...) {
//    failed[t] = 1;
//    errors[t] = "unknown error";
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (failed[t]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_bed_blocks failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[t]
//    );
//   }
//  }
//
//  std::vector<std::vector<std::vector<double>>> result(20);
//
//  for (int k = 0; k < 20; ++k) result[k].resize(nt);
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
//   result[16][t].resize(2);
//   result[17][t].resize(2);
//   result[18][t].resize(4);
//   result[19][t].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   arma::vec y_t(n_used, arma::fill::zeros);
//
//   for (int j = 0; j < n_used; ++j) {
//    y_t(j) = y(j, t);
//   }
//
//   arma::rowvec b_t = b_mat.row(t);
//   arma::vec e_t = y_t - stblr_bed_xb_global(G, blocks, af_cpp, scale, b_t);
//   arma::vec wy_t = stblr_bed_xty(G, blocks, af_cpp, scale, y_t);
//   arma::vec r_t  = stblr_bed_xte(G, blocks, af_cpp, scale, e_t);
//
//   for (int i = 0; i < m; ++i) {
//    result[0][t][i] = bm_mat(t, i);
//    result[1][t][i] = dm_mat(t, i);
//    result[2][t][i] = wy_t(i);
//    result[3][t][i] = r_t(i);
//    result[4][t][i] = b_mat(t, i);
//    result[5][t][i] = static_cast<double>(d_mat(t, i));
//    result[6][t][i] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int it = 0; it < nit + nburn; ++it) {
//    result[7][t][it] = vbs_mat(t, it);
//    result[8][t][it] = vgs_mat(t, it);
//    result[9][t][it] = ves_mat(t, it);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    result[10][t1][t2] = 0.0;
//    result[11][t1][t2] = 0.0;
//    result[12][t1][t2] = 0.0;
//    result[13][t1][t2] = 0.0;
//    result[14][t1][t2] = 0.0;
//    result[15][t1][t2] = 0.0;
//   }
//
//   result[10][t1][t1] = final_vb(t1);
//   result[11][t1][t1] = final_vg(t1);
//   result[12][t1][t1] = final_ve(t1);
//
//   result[13][t1][t1] = final_vb(t1);
//   result[14][t1][t1] = final_vg(t1);
//   result[15][t1][t1] = final_ve(t1);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[16][t][0] = 1.0 - final_pi(t);
//   result[16][t][1] = final_pi(t);
//
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(t, it);
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi(t);
//
//   result[17][t][0] = 1.0 - mean_pi;
//   result[17][t][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) result[18][t][i] = 0.0;
//   for (int i = 0; i < 2; ++i) result[19][t][i] = 0.0;
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
// #include <unordered_map>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // -----------------------------------------------------------------------------
// // Individual-level block structures
// // -----------------------------------------------------------------------------
//
// struct STBlock {
//  std::vector<int> markers;  // global marker indices, 0-based
//  arma::mat XtX;             // block X'X, size p x p
//  arma::vec Xty;             // block X'y, size p
//  arma::vec xx;              // diag(X'X), size p
// };
//
// inline std::vector<STBlock> build_st_blocks_individual(
//   const arma::mat& X,
//   const arma::vec& y,
//   const std::vector<int>& sets
// ) {
//  const int n = static_cast<int>(X.n_rows);
//  const int m = static_cast<int>(X.n_cols);
//
//  if (static_cast<int>(y.n_elem) != n) {
//   throw std::runtime_error("build_st_blocks_individual: y length must equal nrow(X).");
//  }
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("build_st_blocks_individual: sets must have length ncol(X). It should assign each marker to a block.");
//  }
//
//  // Preserve the first-occurrence order of set labels. This keeps adjacent blocks
//  // in the order supplied by the user.
//  std::vector<int> labels;
//  labels.reserve(m);
//  std::unordered_map<int, int> label_to_block;
//
//  for (int j = 0; j < m; ++j) {
//   const int lab = sets[j];
//   if (label_to_block.find(lab) == label_to_block.end()) {
//    const int block_id = static_cast<int>(labels.size());
//    label_to_block[lab] = block_id;
//    labels.push_back(lab);
//   }
//  }
//
//  std::vector<STBlock> blocks(labels.size());
//
//  for (int j = 0; j < m; ++j) {
//   blocks[label_to_block[sets[j]]].markers.push_back(j);
//  }
//
//  for (STBlock& blk : blocks) {
//   const int p = static_cast<int>(blk.markers.size());
//   arma::uvec idx(p);
//
//   for (int k = 0; k < p; ++k) {
//    idx(k) = static_cast<arma::uword>(blk.markers[k]);
//   }
//
//   arma::mat Xb = X.cols(idx);
//   blk.XtX = Xb.t() * Xb;
//   blk.Xty = Xb.t() * y;
//   blk.xx  = arma::diagvec(blk.XtX);
//
//   for (int k = 0; k < p; ++k) {
//    if (!std::isfinite(blk.xx(k)) || blk.xx(k) <= 0.0) {
//     throw std::runtime_error(
//       "build_st_blocks_individual: marker has invalid X'X diagonal in block."
//     );
//    }
//   }
//  }
//
//  return blocks;
// }
//
// inline arma::vec block_xb_from_global_b(
//   const arma::mat& X,
//   const std::vector<int>& markers,
//   const arma::rowvec& b
// ) {
//  const int p = static_cast<int>(markers.size());
//  arma::vec xb(X.n_rows, arma::fill::zeros);
//
//  for (int k = 0; k < p; ++k) {
//   const int j = markers[k];
//   const double bj = b(j);
//   if (bj != 0.0) {
//    xb += X.col(j) * bj;
//   }
//  }
//
//  return xb;
// }
//
// inline arma::vec compute_xb_global(
//   const arma::mat& X,
//   const arma::rowvec& b
// ) {
//  return X * b.t();
// }
//
// // -----------------------------------------------------------------------------
// // Block-local single-trait BayesC marker update
// //
// // Here r = local_wy - local_XtX * b_block.
// // For marker k inside the block:
// //   score = x_k' residual_without_marker_k = r_k + xx_k * b_k
// // -----------------------------------------------------------------------------
//
// inline void sampleBetaC_ST_block_dense(
//   int k,
//   const std::vector<double>& pi,
//   double vb,
//   double vei,
//   const arma::vec& xx,
//   const arma::mat& XtX,
//   arma::vec& r,
//   arma::vec& b_block,
//   arma::Row<int>& d_block,
//   std::mt19937& gen
// ) {
//  const double xk2 = xx(k);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm01(0.0, 1.0);
//
//  const double pi0 = std::max(pi[0], 1e-300);
//  const double pi1 = std::max(pi[1], 1e-300);
//  const double vei_safe = std::max(vei, 1e-300);
//
//  const double score = r(k) + xk2 * b_block(k);
//  const double denom = std::max(vei_safe + xk2 * vb, 1e-300);
//
//  const double logBF =
//   0.5 * std::log(vei_safe / denom) +
//   0.5 * score * score * vb / (vei_safe * denom);
//
//  const double logp1 = std::log(pi1) + logBF;
//  const double logp0 = std::log(pi0);
//  const double delta_log = logp0 - logp1;
//
//  double p1 = 0.0;
//  if (delta_log > 35.0) {
//   p1 = 0.0;
//  } else if (delta_log < -35.0) {
//   p1 = 1.0;
//  } else {
//   p1 = 1.0 / (1.0 + std::exp(delta_log));
//  }
//
//  const int dk = (runif(gen) < p1) ? 1 : 0;
//
//  double b_new = 0.0;
//  if (dk == 1) {
//   const double lhs = xk2 + vei_safe / vb;
//   const double mean = score / lhs;
//   const double sd = std::sqrt(vei_safe / lhs);
//   b_new = mean + sd * norm01(gen);
//  }
//
//  const double diff = b_new - b_block(k);
//
//  if (diff != 0.0) {
//   // r := local_wy - XtX * b_block
//   // Changing b_k by diff subtracts XtX[, k] * diff from r.
//   r -= XtX.col(k) * diff;
//  }
//
//  b_block(k) = b_new;
//  d_block(k) = dk;
// }
//
// // -----------------------------------------------------------------------------
// // Variance and pi updates for individual-level residuals
// // -----------------------------------------------------------------------------
//
// inline void sampleB_ST_individual(
//   int m,
//   double nub,
//   double& vb,
//   const arma::rowvec& b,
//   const arma::Row<int>& d,
//   double ssb_prior,
//   std::mt19937& gen
// ) {
//  double ssb = 0.0;
//  double dfb = 0.0;
//
//  for (int i = 0; i < m; ++i) {
//   if (d(i) > 0) {
//    ssb += b(i) * b(i);
//    dfb += 1.0;
//   }
//  }
//
//  const double scale = ssb + nub * ssb_prior;
//  std::chi_squared_distribution<double> rchisq(dfb + nub);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  vb = std::max(scale / chi2, 1e-12);
// }
//
// inline void sampleE_ST_individual(
//   double nue,
//   double& ve,
//   const arma::vec& e,
//   double sse_prior,
//   std::mt19937& gen
// ) {
//  const double sse = arma::dot(e, e);
//  const double scale = sse + nue * sse_prior;
//
//  if (!std::isfinite(scale) || scale <= 0.0) {
//   throw std::runtime_error("sampleE_ST_individual: invalid residual scale.");
//  }
//
//  std::chi_squared_distribution<double> rchisq(static_cast<double>(e.n_elem) + nue);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//
//  const double ve_new = scale / chi2;
//  if (!std::isfinite(ve_new) || ve_new <= 0.0) {
//   throw std::runtime_error("sampleE_ST_individual: sampled ve is invalid.");
//  }
//
//  ve = std::max(ve_new, 1e-12);
// }
//
// inline double computeG_ST_individual(
//   const arma::mat& X,
//   const arma::rowvec& b
// ) {
//  arma::vec g = X * b.t();
//  return arma::dot(g, g) / static_cast<double>(X.n_rows);
// }
//
// inline void samplePi_ST_individual(
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
// inline arma::uvec make_order_within_block(
//   const arma::vec& local_wy,
//   const arma::vec& xx
// ) {
//  const int p = static_cast<int>(local_wy.n_elem);
//  std::vector<int> ord(p);
//  std::iota(ord.begin(), ord.end(), 0);
//
//  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
//   const double bha = local_wy(a) / xx(a);
//   const double bhb = local_wy(b) / xx(b);
//   return bha * bha > bhb * bhb;
//  });
//
//  arma::uvec out(p);
//  for (int k = 0; k < p; ++k) out(k) = static_cast<arma::uword>(ord[k]);
//  return out;
// }
//
// // -----------------------------------------------------------------------------
// // Main exported function: individual-level block STBLR
// //
// // X_list: list of n_t x m genotype/covariate matrices, one per trait.
// // y_list: list of length n_t response vectors, one per trait.
// // sets:   length m vector assigning each adjacent marker to a block.
// //         Example in R for block_size=1000:
// //         sets <- rep(seq_len(ceiling(m / block_size)), each = block_size)[seq_len(m)]
// //
// // Notes:
// // - X and y should be centered/standardized consistently before calling this.
// // - The sampler keeps an individual-level residual e = y - Xb.
// // - For each block, it computes local_wy = X_block' * residual_without_block,
// //   then runs block_nit local BayesC sweeps using the cached block XtX.
// // - The individual-level residual is updated once per block using X_block * db.
// // -----------------------------------------------------------------------------
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_individual_blocks(
//   Rcpp::List X_list,
//   Rcpp::List y_list,
//   std::vector<std::vector<double>> b_init,
//   std::vector<int> sets,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   int nit,
//   int nburn,
//   int nthin,
//   int block_nit,
//   int ncores,
//   int seed
// ) {
//  const int nt = static_cast<int>(X_list.size());
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: X_list must contain at least one trait.");
//  }
//
//  if (static_cast<int>(y_list.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: y_list length must equal X_list length.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: nthin must be positive.");
//  }
//
//  if (block_nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: block_nit must be positive.");
//  }
//
//  if (static_cast<int>(b_init.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: b_init length must equal number of traits.");
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: pi must have length 2, c(pi0, pi1).");
//  }
//
//  if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: B must be nt x nt.");
//  }
//
//  if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: E must be nt x nt.");
//  }
//
//  arma::mat X0 = Rcpp::as<arma::mat>(X_list[0]);
//  const int m = static_cast<int>(X0.n_cols);
//
//  if (m <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: ncol(X) must be positive.");
//  }
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_individual_blocks: sets must have length ncol(X).");
//  }
//
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   arma::mat Xt = Rcpp::as<arma::mat>(X_list[t]);
//   arma::vec yt = Rcpp::as<arma::vec>(y_list[t]);
//
//   if (static_cast<int>(Xt.n_cols) != m) {
//    throw std::runtime_error("stblr_cpg_omp_individual_blocks: all X matrices must have same number of columns.");
//   }
//
//   if (static_cast<int>(yt.n_elem) != static_cast<int>(Xt.n_rows)) {
//    throw std::runtime_error("stblr_cpg_omp_individual_blocks: y length must equal nrow(X) for each trait.");
//   }
//
//   if (static_cast<int>(b_init[t].size()) != m) {
//    throw std::runtime_error("stblr_cpg_omp_individual_blocks: each b_init[t] must have length m.");
//   }
//
//   if (static_cast<int>(ssb_prior[t].size()) != nt || static_cast<int>(sse_prior[t].size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_individual_blocks: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//  }
//
//  // Output storage.
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
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
//  arma::vec yy_vec(nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) b_mat(t, i) = b_init[t][i];
//  }
//
// #ifdef _OPENMP
//  omp_set_dynamic(0);
//  omp_set_num_threads(ncores);
//  Rcpp::Rcout
//  << "STBLR individual-level OpenMP max threads = "
//  << omp_get_max_threads()
//  << ", num procs = "
//  << omp_get_num_procs()
//  << "\n";
// #endif
//
//  std::vector<int> failed(nt, 0);
//  std::vector<std::string> errors(nt);
//
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
//   try {
// #ifdef _OPENMP
// #pragma omp critical
// {
//  Rcpp::Rcout
//  << "trait " << t
//  << " running on thread "
//  << omp_get_thread_num()
//  << " of "
//  << omp_get_num_threads()
//  << "\n";
// }
// #endif
//
//    std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//    arma::mat X_t = Rcpp::as<arma::mat>(X_list[t]);
//    arma::vec y_t = Rcpp::as<arma::vec>(y_list[t]);
//    const int n_t = static_cast<int>(X_t.n_rows);
//
//    yy_vec(t) = arma::dot(y_t, y_t);
//
//    arma::rowvec b_t = b_mat.row(t);
//    arma::Row<int> d_t(m, arma::fill::zeros);
//
//    // Initial individual residual e = y - Xb.
//    arma::vec e_t = y_t - compute_xb_global(X_t, b_t);
//
//    double vb_t = B(t, t);
//    double ve_t = E(t, t);
//    double vg_t = computeG_ST_individual(X_t, b_t);
//    double vei_t = ve_t + adjE * vg_t;
//
//    std::vector<double> pi_t = pi;
//    if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//        pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//     throw std::runtime_error("invalid initial pi.");
//    }
//
//    {
//     const double psum = pi_t[0] + pi_t[1];
//     if (!std::isfinite(psum) || psum <= 0.0) {
//      throw std::runtime_error("invalid initial pi sum.");
//     }
//     pi_t[0] /= psum;
//     pi_t[1] /= psum;
//    }
//
//    // Build and store block-level X'X and X'y in memory for this trait.
//    std::vector<STBlock> blocks = build_st_blocks_individual(X_t, y_t, sets);
//
//    // Fixed block order: first occurrence order in sets. Within each block we
//    // recompute the local marker order from residual_without_block at each visit.
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//     for (const STBlock& blk : blocks) {
//      const int p = static_cast<int>(blk.markers.size());
//
//      arma::vec b_old(p, arma::fill::zeros);
//      arma::vec b_block(p, arma::fill::zeros);
//      arma::Row<int> d_block(p, arma::fill::zeros);
//
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[k];
//       b_old(k) = b_t(j);
//       b_block(k) = b_t(j);
//       d_block(k) = d_t(j);
//      }
//
//      // residual_without_block = e + X_block * b_old.
//      arma::vec xb_old = block_xb_from_global_b(X_t, blk.markers, b_t);
//      arma::vec e_without_block = e_t + xb_old;
//
//      // local_wy = X_block' * residual_without_block.
//      // Equivalent to blk.Xty - X_block'X_outside b_outside, but avoids
//      // storing cross-block products.
//      arma::vec local_wy(p, arma::fill::zeros);
//      for (int k = 0; k < p; ++k) {
//       local_wy(k) = arma::dot(X_t.col(blk.markers[k]), e_without_block);
//      }
//
//      // r_block = local_wy - XtX_block * b_block.
//      arma::vec r_block = local_wy - blk.XtX * b_block;
//
//      for (int bit = 0; bit < block_nit; ++bit) {
//       arma::uvec order = make_order_within_block(local_wy, blk.xx);
//
//       for (arma::uword kk = 0; kk < order.n_elem; ++kk) {
//        const int k = static_cast<int>(order(kk));
//
//        sampleBetaC_ST_block_dense(
//         k,
//         pi_t,
//         vb_t,
//         vei_t,
//         blk.xx,
//         blk.XtX,
//         r_block,
//         b_block,
//         d_block,
//         gen_t
//        );
//       }
//      }
//
//      // Write back block b and d.
//      arma::vec db = b_block - b_old;
//
//      bool changed = false;
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[k];
//       if (db(k) != 0.0) changed = true;
//       b_t(j) = b_block(k);
//       d_t(j) = d_block(k);
//      }
//
//      // Update individual residual once per block: e := e - X_block * db.
//      if (changed) {
//       for (int k = 0; k < p; ++k) {
//        if (db(k) != 0.0) {
//         e_t -= X_t.col(blk.markers[k]) * db(k);
//        }
//       }
//      }
//     }
//
//     if (updateB) {
//      sampleB_ST_individual(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       ssb_prior_mat(t, t),
//       gen_t
//      );
//
//      if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//       throw std::runtime_error("vb became invalid after sampleB.");
//      }
//     }
//
//     if (updateE) {
//      // Rebuild residual before sampleE to avoid accumulated numerical drift.
//      e_t = y_t - compute_xb_global(X_t, b_t);
//
//      sampleE_ST_individual(
//       nue,
//       ve_t,
//       e_t,
//       sse_prior_mat(t, t),
//       gen_t
//      );
//
//      if (!std::isfinite(ve_t) || ve_t <= 0.0) {
//       throw std::runtime_error("ve became invalid after sampleE.");
//      }
//     }
//
//     if (updatePi) {
//      samplePi_ST_individual(d_t, pi_t, gen_t);
//
//      if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
//          pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//       throw std::runtime_error("pi became invalid after samplePi.");
//      }
//     }
//
//     vg_t = computeG_ST_individual(X_t, b_t);
//     vei_t = ve_t + adjE * vg_t;
//
//     if (!std::isfinite(vg_t)) {
//      throw std::runtime_error("vg became NaN/Inf after computeG.");
//     }
//
//     vbs_mat(t, it) = vb_t;
//     ves_mat(t, it) = ve_t;
//     vgs_mat(t, it) = vg_t;
//     pis_mat(t, it) = pi_t[1];
//
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//
//      for (int i = 0; i < m; ++i) {
//       bm_mat(t, i) += b_t(i);
//       dm_mat(t, i) += static_cast<double>(d_t(i));
//      }
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//
//    for (int i = 0; i < m; ++i) {
//     bm_mat(t, i) /= nsamples_t;
//     dm_mat(t, i) /= nsamples_t;
//    }
//
//    if (!bm_mat.row(t).is_finite()) {
//     throw std::runtime_error("posterior mean bm contains NaN/Inf.");
//    }
//
//    if (!dm_mat.row(t).is_finite()) {
//     throw std::runtime_error("posterior mean dm contains NaN/Inf.");
//    }
//
//    b_mat.row(t) = b_t;
//    d_mat.row(t) = d_t;
//
//    final_vb(t) = vb_t;
//    final_ve(t) = ve_t;
//    final_vg(t) = vg_t;
//    final_pi(t) = pi_t[1];
//    nsamples_vec(t) = nsamples_t;
//
//   } catch (const std::exception& e) {
//    failed[t] = 1;
//    errors[t] = e.what();
//   } catch (...) {
//    failed[t] = 1;
//    errors[t] = "unknown error";
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (failed[t]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_individual_blocks failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[t]
//    );
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Build result with same 20-slot style as your existing ST output.
//  // Slot 2 now stores X'y computed from individual data, not input wy.
//  // Slot 3 stores X'e final residual, computed block-free at the end.
//  // --------------------------------------------------------------------------
//
//  std::vector<std::vector<std::vector<double>>> result(20);
//  for (int k = 0; k < 20; ++k) result[k].resize(nt);
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][t].resize(m);             // bm
//   result[1][t].resize(m);             // dm
//   result[2][t].resize(m);             // wy = X'y
//   result[3][t].resize(m);             // r = X'e final residual
//   result[4][t].resize(m);             // final b
//   result[5][t].resize(m);             // final d
//   result[6][t].resize(m);             // marker index
//
//   result[7][t].resize(nit + nburn);   // vbs
//   result[8][t].resize(nit + nburn);   // vgs
//   result[9][t].resize(nit + nburn);   // ves
//
//   result[10][t].resize(nt);           // covb
//   result[11][t].resize(nt);           // covg
//   result[12][t].resize(nt);           // cove
//   result[13][t].resize(nt);           // final B
//   result[14][t].resize(nt);           // final G
//   result[15][t].resize(nt);           // final E
//
//   result[16][t].resize(2);            // final pi
//   result[17][t].resize(2);            // posterior mean pi approx
//   result[18][t].resize(4);
//   result[19][t].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   arma::mat X_t = Rcpp::as<arma::mat>(X_list[t]);
//   arma::vec y_t = Rcpp::as<arma::vec>(y_list[t]);
//   arma::rowvec b_t = b_mat.row(t);
//   arma::vec e_t = y_t - compute_xb_global(X_t, b_t);
//   arma::vec wy_t = X_t.t() * y_t;
//   arma::vec r_t  = X_t.t() * e_t;
//
//   for (int i = 0; i < m; ++i) {
//    result[0][t][i] = bm_mat(t, i);
//    result[1][t][i] = dm_mat(t, i);
//    result[2][t][i] = wy_t(i);
//    result[3][t][i] = r_t(i);
//    result[4][t][i] = b_mat(t, i);
//    result[5][t][i] = static_cast<double>(d_mat(t, i));
//    result[6][t][i] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int it = 0; it < nit + nburn; ++it) {
//    result[7][t][it] = vbs_mat(t, it);
//    result[8][t][it] = vgs_mat(t, it);
//    result[9][t][it] = ves_mat(t, it);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    result[10][t1][t2] = 0.0;
//    result[11][t1][t2] = 0.0;
//    result[12][t1][t2] = 0.0;
//    result[13][t1][t2] = 0.0;
//    result[14][t1][t2] = 0.0;
//    result[15][t1][t2] = 0.0;
//   }
//
//   result[10][t1][t1] = final_vb(t1);
//   result[11][t1][t1] = final_vg(t1);
//   result[12][t1][t1] = final_ve(t1);
//   result[13][t1][t1] = final_vb(t1);
//   result[14][t1][t1] = final_vg(t1);
//   result[15][t1][t1] = final_ve(t1);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[16][t][0] = 1.0 - final_pi(t);
//   result[16][t][1] = final_pi(t);
//
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(t, it);
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi(t);
//
//   result[17][t][0] = 1.0 - mean_pi;
//   result[17][t][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) result[18][t][i] = 0.0;
//   for (int i = 0; i < 2; ++i) result[19][t][i] = 0.0;
//  }
//
//  return result;
// }
//
// // =============================================================================
// // BED-backed individual-level block STBLR extension
// // =============================================================================
// // This version avoids storing the full dense X matrix. Genotypes are kept in the
// // compact PackedBedMatrix representation from packed_bed.h and decoded block-wise
// // when computing
// //   - block X'X cached in memory,
// //   - block X'e for the current individual-level residual,
// //   - block residual updates e <- e - X_block * db.
// //
// // It uses the same BED coding and missing handling as bed_xtx_diag_xty_core():
// //   00 -> dosage 2
// //   01 -> missing
// //   10 -> dosage 1
// //   11 -> dosage 0
// // If scale = true:
// //   x = (dosage - 2p) / sqrt(2p(1-p)); missing -> 0
// // If scale = false:
// //   raw dosage; missing -> 2p
// //
// // Required external header from your sparseLD code:
// //   packed_bed.h
// // providing PackedBedMatrix, read_bedfiles_to_packed_matrix(),
// // and compute_af_from_packed().
// // =============================================================================
//
// #include "packed_bed.h"
//
// struct STBedBlock {
//  std::vector<int> markers;  // global marker indices, 0-based
//  arma::mat XtX;             // p x p cached block cross-product
//  arma::vec xx;              // diag(XtX)
// };
//
// static std::vector<std::string> stblr_copy_bed_files(Rcpp::CharacterVector bed_files) {
//  std::vector<std::string> out(bed_files.size());
//  for (int i = 0; i < bed_files.size(); ++i) {
//   out[static_cast<std::size_t>(i)] = Rcpp::as<std::string>(bed_files[i]);
//  }
//  return out;
// }
//
// static std::vector<std::vector<int>> stblr_copy_int_list(Rcpp::List xlist) {
//  std::vector<std::vector<int>> out(xlist.size());
//  for (int f = 0; f < xlist.size(); ++f) {
//   Rcpp::IntegerVector x = xlist[f];
//   out[static_cast<std::size_t>(f)].resize(x.size());
//   for (int i = 0; i < x.size(); ++i) {
//    out[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)] = x[i];
//   }
//  }
//  return out;
// }
//
// static std::vector<int> stblr_copy_rows0_or_empty(
//   Rcpp::Nullable<Rcpp::IntegerVector> rows,
//   int n
// ) {
//  std::vector<int> out;
//  if (rows.isNotNull()) {
//   Rcpp::IntegerVector r(rows);
//   out.resize(r.size());
//   for (int i = 0; i < r.size(); ++i) {
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
// static std::vector<double> stblr_flatten_af_list_or_empty(
//   Rcpp::Nullable<Rcpp::List> af
// ) {
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
// static inline void stblr_bed_maps(
//   double p,
//   bool scale,
//   double& m0,
//   double& m1,
//   double& m2,
//   double& m3
// ) {
//  if (scale) {
//   const double denom = std::sqrt(2.0 * p * (1.0 - p));
//   if (denom <= 0.0 || !std::isfinite(denom)) {
//    m0 = m1 = m2 = m3 = 0.0;
//   } else {
//    m0 = (2.0 - 2.0 * p) / denom;
//    m1 = 0.0;
//    m2 = (1.0 - 2.0 * p) / denom;
//    m3 = (0.0 - 2.0 * p) / denom;
//   }
//  } else {
//   m0 = 2.0;
//   m1 = 2.0 * p;
//   m2 = 1.0;
//   m3 = 0.0;
//  }
// }
//
// static inline double stblr_decode_2bit_value(
//   unsigned int code,
//   double m0,
//   double m1,
//   double m2,
//   double m3
// ) {
//  switch (code) {
//  case 0u: return m0;
//  case 1u: return m1;
//  case 2u: return m2;
//  default: return m3;
//  }
// }
//
// static arma::mat stblr_decode_bed_block_double(
//   const PackedBedMatrix& G,
//   const std::vector<int>& markers,
//   const std::vector<double>& af,
//   bool scale
// ) {
//  const int n = G.n;
//  const int p = static_cast<int>(markers.size());
//  const std::size_t nbytes = G.nbytes;
//
//  arma::mat X(n, p, arma::fill::zeros);
//
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
//  for (int kk = 0; kk < p; ++kk) {
//   const int global_i = markers[static_cast<std::size_t>(kk)];
//   const uint8_t* packed = G.row(global_i);
//
//   double m0, m1, m2, m3;
//   stblr_bed_maps(af[static_cast<std::size_t>(global_i)], scale, m0, m1, m2, m3);
//
//   for (std::size_t kb = 0; kb < nbytes; ++kb) {
//    unsigned char x = packed[kb];
//    const int jbase = static_cast<int>(kb << 2);
//
//    if (jbase + 0 < n) X(jbase + 0, kk) = stblr_decode_2bit_value((x >> 0) & 3u, m0, m1, m2, m3);
//    if (jbase + 1 < n) X(jbase + 1, kk) = stblr_decode_2bit_value((x >> 2) & 3u, m0, m1, m2, m3);
//    if (jbase + 2 < n) X(jbase + 2, kk) = stblr_decode_2bit_value((x >> 4) & 3u, m0, m1, m2, m3);
//    if (jbase + 3 < n) X(jbase + 3, kk) = stblr_decode_2bit_value((x >> 6) & 3u, m0, m1, m2, m3);
//   }
//  }
//
//  return X;
// }
//
// static std::vector<STBedBlock> stblr_build_bed_blocks(
//   const PackedBedMatrix& G,
//   const std::vector<int>& sets,
//   const std::vector<double>& af,
//   bool scale
// ) {
//  const int m = G.m;
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_build_bed_blocks: sets must have length equal to number of markers.");
//  }
//
//  std::vector<int> labels;
//  labels.reserve(m);
//  std::unordered_map<int, int> label_to_block;
//
//  for (int j = 0; j < m; ++j) {
//   const int lab = sets[static_cast<std::size_t>(j)];
//   if (label_to_block.find(lab) == label_to_block.end()) {
//    const int block_id = static_cast<int>(labels.size());
//    label_to_block[lab] = block_id;
//    labels.push_back(lab);
//   }
//  }
//
//  std::vector<STBedBlock> blocks(labels.size());
//
//  for (int j = 0; j < m; ++j) {
//   blocks[static_cast<std::size_t>(label_to_block[sets[static_cast<std::size_t>(j)]])].markers.push_back(j);
//  }
//
//  for (STBedBlock& blk : blocks) {
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af, scale);
//   blk.XtX = Xb.t() * Xb;
//   blk.xx = arma::diagvec(blk.XtX);
//
//   for (arma::uword k = 0; k < blk.xx.n_elem; ++k) {
//    if (!std::isfinite(blk.xx(k)) || blk.xx(k) <= 0.0) {
//     throw std::runtime_error("stblr_build_bed_blocks: invalid X'X diagonal. Check allele frequencies and monomorphic markers.");
//    }
//   }
//  }
//
//  return blocks;
// }
//
// static arma::vec stblr_bed_xb_global(
//   const PackedBedMatrix& G,
//   const std::vector<STBedBlock>& blocks,
//   const std::vector<double>& af,
//   bool scale,
//   const arma::rowvec& b
// ) {
//  arma::vec xb(G.n, arma::fill::zeros);
//
//  for (const STBedBlock& blk : blocks) {
//   arma::vec b_block(blk.markers.size(), arma::fill::zeros);
//   bool any = false;
//
//   for (std::size_t k = 0; k < blk.markers.size(); ++k) {
//    const int j = blk.markers[k];
//    b_block(k) = b(j);
//    if (b_block(k) != 0.0) any = true;
//   }
//
//   if (!any) continue;
//
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af, scale);
//   xb += Xb * b_block;
//  }
//
//  return xb;
// }
//
// static arma::vec stblr_bed_xty(
//   const PackedBedMatrix& G,
//   const std::vector<STBedBlock>& blocks,
//   const std::vector<double>& af,
//   bool scale,
//   const arma::vec& y
// ) {
//  arma::vec xty(G.m, arma::fill::zeros);
//
//  for (const STBedBlock& blk : blocks) {
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af, scale);
//   arma::vec xty_block = Xb.t() * y;
//
//   for (std::size_t k = 0; k < blk.markers.size(); ++k) {
//    xty(blk.markers[k]) = xty_block(k);
//   }
//  }
//
//  return xty;
// }
//
// static arma::vec stblr_bed_xte(
//   const PackedBedMatrix& G,
//   const std::vector<STBedBlock>& blocks,
//   const std::vector<double>& af,
//   bool scale,
//   const arma::vec& e
// ) {
//  arma::vec xte(G.m, arma::fill::zeros);
//
//  for (const STBedBlock& blk : blocks) {
//   arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af, scale);
//   arma::vec xte_block = Xb.t() * e;
//
//   for (std::size_t k = 0; k < blk.markers.size(); ++k) {
//    xte(blk.markers[k]) = xte_block(k);
//   }
//  }
//
//  return xte;
// }
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_bed_blocks(
//   Rcpp::CharacterVector bed_files,
//   int n,
//   Rcpp::List cls,
//   Rcpp::NumericMatrix y,
//   std::vector<std::vector<double>> b_init,
//   std::vector<int> sets,
//   Rcpp::Nullable<Rcpp::IntegerVector> rows,
//   Rcpp::Nullable<Rcpp::List> af,
//   bool scale,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   double adjE,
//   int nit,
//   int nburn,
//   int nthin,
//   int block_nit,
//   int ncores,
//   int seed
// ) {
//  if (nit <= 0 || nburn < 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: nit must be positive and nburn non-negative.");
//  }
//
//  if (nthin <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: nthin must be positive.");
//  }
//
//  if (block_nit <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: block_nit must be positive.");
//  }
//
//  std::vector<std::string> bed_files_cpp = stblr_copy_bed_files(bed_files);
//  std::vector<std::vector<int>> cls_by_file = stblr_copy_int_list(cls);
//  std::vector<int> rows0 = stblr_copy_rows0_or_empty(rows, n);
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
//  const int n_used = G.n;
//  const int m = G.m;
//  const int nt = y.ncol();
//
//  if (nt <= 0) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: y must have at least one trait column.");
//  }
//
//  if (y.nrow() != n_used) {
//   throw std::runtime_error(
//     "stblr_cpg_omp_bed_blocks: y rows must equal the number of samples used after rows filtering. "
//     "If rows is supplied, pass y already subset to those rows."
//   );
//  }
//
//  if (static_cast<int>(sets.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: sets must have length equal to the total number of BED markers used.");
//  }
//
//  if (static_cast<int>(b_init.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: b_init must have length nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(b_init[static_cast<std::size_t>(t)].size()) != m) {
//    throw std::runtime_error("stblr_cpg_omp_bed_blocks: each b_init[t] must have length m.");
//   }
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: pi must have length 2.");
//  }
//
//  // Keep B and E explicit in the exported signature. Rcpp::compileAttributes()
//  // cannot parse defaults such as arma::mat(), and if it partially parses the
//  // function declaration it can generate a wrapper with the wrong C++ signature.
//
//  if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: B must be nt x nt.");
//  }
//
//  if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: E must be nt x nt.");
//  }
//
//  // Keep priors explicit for the same reason as B/E. If you want R-level
//  // defaults, define them in an R wrapper instead of in the exported C++
//  // function signature.
//
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  if (static_cast<int>(ssb_prior.size()) != nt || static_cast<int>(sse_prior.size()) != nt) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: prior lists must have length nt.");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (static_cast<int>(ssb_prior[t].size()) != nt || static_cast<int>(sse_prior[t].size()) != nt) {
//    throw std::runtime_error("stblr_cpg_omp_bed_blocks: priors must be nt x nt.");
//   }
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//  }
//
//  std::vector<double> af_cpp = stblr_flatten_af_list_or_empty(af);
//  const bool af_computed = af_cpp.empty();
//  if (af_computed) af_cpp = compute_af_from_packed(G);
//
//  if (static_cast<int>(af_cpp.size()) != m) {
//   throw std::runtime_error("stblr_cpg_omp_bed_blocks: af must have one value per marker after flattening.");
//  }
//
// #ifdef _OPENMP
//  if (ncores > 0) {
//   omp_set_dynamic(0);
//   omp_set_num_threads(ncores);
//  }
// #endif
//
//  Rcpp::Rcout << "Building BED-backed STBLR blocks: n=" << n_used
//              << ", m=" << m
//              << ", nt=" << nt
//              << ", scale=" << scale
//              << ", af_computed=" << af_computed
//              << "\n";
//
//  // Cached genotype block cross-products are shared across traits.
//  std::vector<STBedBlock> blocks = stblr_build_bed_blocks(G, sets, af_cpp, scale);
//
//  Rcpp::Rcout << "Number of STBLR blocks = " << blocks.size() << "\n";
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
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
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) b_mat(t, i) = b_init[t][i];
//  }
//
//  std::vector<int> failed(nt, 0);
//  std::vector<std::string> errors(nt);
//
// #ifdef _OPENMP
// #pragma omp parallel for schedule(static)
// #endif
//  for (int t = 0; t < nt; ++t) {
//   try {
//    std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));
//
//    arma::vec y_t(n_used, arma::fill::zeros);
//    for (int j = 0; j < n_used; ++j) y_t(j) = y(j, t);
//
//    arma::rowvec b_t = b_mat.row(t);
//    arma::Row<int> d_t(m, arma::fill::zeros);
//
//    arma::vec e_t = y_t - stblr_bed_xb_global(G, blocks, af_cpp, scale, b_t);
//
//    double vb_t = B(t, t);
//    double ve_t = E(t, t);
//    double vg_t = arma::dot(y_t - e_t, y_t - e_t) / static_cast<double>(n_used);
//    double vei_t = ve_t + adjE * vg_t;
//
//    std::vector<double> pi_t = pi;
//    {
//     const double psum = pi_t[0] + pi_t[1];
//     if (!std::isfinite(psum) || psum <= 0.0 || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//      throw std::runtime_error("invalid initial pi.");
//     }
//     pi_t[0] /= psum;
//     pi_t[1] /= psum;
//    }
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < nit + nburn; ++it) {
//     for (const STBedBlock& blk : blocks) {
//      const int p = static_cast<int>(blk.markers.size());
//
//      arma::vec b_old(p, arma::fill::zeros);
//      arma::vec b_block(p, arma::fill::zeros);
//      arma::Row<int> d_block(p, arma::fill::zeros);
//
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[static_cast<std::size_t>(k)];
//       b_old(k) = b_t(j);
//       b_block(k) = b_t(j);
//       d_block(k) = d_t(j);
//      }
//
//      // Decode only this block from the packed BED representation.
//      arma::mat Xb = stblr_decode_bed_block_double(G, blk.markers, af_cpp, scale);
//
//      // e_without_block = e + X_block * b_old.
//      arma::vec xb_old = Xb * b_old;
//      arma::vec e_without_block = e_t + xb_old;
//
//      // local_wy is X_block' residual_without_block.
//      arma::vec local_wy = Xb.t() * e_without_block;
//      arma::vec r_block = local_wy - blk.XtX * b_block;
//
//      for (int bit = 0; bit < block_nit; ++bit) {
//       arma::uvec order = make_order_within_block(local_wy, blk.xx);
//
//       for (arma::uword kk = 0; kk < order.n_elem; ++kk) {
//        const int k = static_cast<int>(order(kk));
//
//        sampleBetaC_ST_block_dense(
//         k,
//         pi_t,
//         vb_t,
//         vei_t,
//         blk.xx,
//         blk.XtX,
//         r_block,
//         b_block,
//         d_block,
//         gen_t
//        );
//       }
//      }
//
//      arma::vec db = b_block - b_old;
//
//      bool changed = false;
//      for (int k = 0; k < p; ++k) {
//       const int j = blk.markers[static_cast<std::size_t>(k)];
//       if (db(k) != 0.0) changed = true;
//       b_t(j) = b_block(k);
//       d_t(j) = d_block(k);
//      }
//
//      // Direct BED-backed residual update using the decoded current block.
//      if (changed) {
//       e_t -= Xb * db;
//      }
//     }
//
//     if (updateB) {
//      sampleB_ST_individual(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       ssb_prior_mat(t, t),
//       gen_t
//      );
//     }
//
//     if (updateE) {
//      // Rebuild from BED-backed Xb occasionally/each iteration to avoid drift.
//      e_t = y_t - stblr_bed_xb_global(G, blocks, af_cpp, scale, b_t);
//
//      sampleE_ST_individual(
//       nue,
//       ve_t,
//       e_t,
//       sse_prior_mat(t, t),
//       gen_t
//      );
//     }
//
//     if (updatePi) {
//      samplePi_ST_individual(d_t, pi_t, gen_t);
//     }
//
//     arma::vec g_t = y_t - e_t;
//     vg_t = arma::dot(g_t, g_t) / static_cast<double>(n_used);
//     vei_t = ve_t + adjE * vg_t;
//
//     if (!std::isfinite(vb_t) || vb_t <= 0.0) throw std::runtime_error("invalid vb.");
//     if (!std::isfinite(ve_t) || ve_t <= 0.0) throw std::runtime_error("invalid ve.");
//     if (!std::isfinite(vg_t)) throw std::runtime_error("invalid vg.");
//     if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//      throw std::runtime_error("invalid pi.");
//     }
//
//     vbs_mat(t, it) = vb_t;
//     ves_mat(t, it) = ve_t;
//     vgs_mat(t, it) = vg_t;
//     pis_mat(t, it) = pi_t[1];
//
//     if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
//      nsamples_t += 1.0;
//      for (int i = 0; i < m; ++i) {
//       bm_mat(t, i) += b_t(i);
//       dm_mat(t, i) += static_cast<double>(d_t(i));
//      }
//     }
//    }
//
//    if (nsamples_t <= 0.0) nsamples_t = 1.0;
//    for (int i = 0; i < m; ++i) {
//     bm_mat(t, i) /= nsamples_t;
//     dm_mat(t, i) /= nsamples_t;
//    }
//
//    b_mat.row(t) = b_t;
//    d_mat.row(t) = d_t;
//
//    final_vb(t) = vb_t;
//    final_ve(t) = ve_t;
//    final_vg(t) = vg_t;
//    final_pi(t) = pi_t[1];
//    nsamples_vec(t) = nsamples_t;
//
//   } catch (const std::exception& e) {
//    failed[t] = 1;
//    errors[t] = e.what();
//   } catch (...) {
//    failed[t] = 1;
//    errors[t] = "unknown error";
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if (failed[t]) {
//    throw std::runtime_error(
//      "stblr_cpg_omp_bed_blocks failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[t]
//    );
//   }
//  }
//
//  std::vector<std::vector<std::vector<double>>> result(20);
//  for (int k = 0; k < 20; ++k) result[k].resize(nt);
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][t].resize(m);
//   result[1][t].resize(m);
//   result[2][t].resize(m);
//   result[3][t].resize(m);
//   result[4][t].resize(m);
//   result[5][t].resize(m);
//   result[6][t].resize(m);
//   result[7][t].resize(nit + nburn);
//   result[8][t].resize(nit + nburn);
//   result[9][t].resize(nit + nburn);
//   result[10][t].resize(nt);
//   result[11][t].resize(nt);
//   result[12][t].resize(nt);
//   result[13][t].resize(nt);
//   result[14][t].resize(nt);
//   result[15][t].resize(nt);
//   result[16][t].resize(2);
//   result[17][t].resize(2);
//   result[18][t].resize(4);
//   result[19][t].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   arma::vec y_t(n_used, arma::fill::zeros);
//   for (int j = 0; j < n_used; ++j) y_t(j) = y(j, t);
//
//   arma::rowvec b_t = b_mat.row(t);
//   arma::vec e_t = y_t - stblr_bed_xb_global(G, blocks, af_cpp, scale, b_t);
//   arma::vec wy_t = stblr_bed_xty(G, blocks, af_cpp, scale, y_t);
//   arma::vec r_t  = stblr_bed_xte(G, blocks, af_cpp, scale, e_t);
//
//   for (int i = 0; i < m; ++i) {
//    result[0][t][i] = bm_mat(t, i);
//    result[1][t][i] = dm_mat(t, i);
//    result[2][t][i] = wy_t(i);
//    result[3][t][i] = r_t(i);
//    result[4][t][i] = b_mat(t, i);
//    result[5][t][i] = static_cast<double>(d_mat(t, i));
//    result[6][t][i] = static_cast<double>(i);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int it = 0; it < nit + nburn; ++it) {
//    result[7][t][it] = vbs_mat(t, it);
//    result[8][t][it] = vgs_mat(t, it);
//    result[9][t][it] = ves_mat(t, it);
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    result[10][t1][t2] = 0.0;
//    result[11][t1][t2] = 0.0;
//    result[12][t1][t2] = 0.0;
//    result[13][t1][t2] = 0.0;
//    result[14][t1][t2] = 0.0;
//    result[15][t1][t2] = 0.0;
//   }
//
//   result[10][t1][t1] = final_vb(t1);
//   result[11][t1][t1] = final_vg(t1);
//   result[12][t1][t1] = final_ve(t1);
//   result[13][t1][t1] = final_vb(t1);
//   result[14][t1][t1] = final_vg(t1);
//   result[15][t1][t1] = final_ve(t1);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   result[16][t][0] = 1.0 - final_pi(t);
//   result[16][t][1] = final_pi(t);
//
//   double mean_pi = 0.0;
//   int npi = 0;
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(t, it);
//    ++npi;
//   }
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi(t);
//
//   result[17][t][0] = 1.0 - mean_pi;
//   result[17][t][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) result[18][t][i] = 0.0;
//   for (int i = 0; i < 2; ++i) result[19][t][i] = 0.0;
//  }
//
//  return result;
// }
//
