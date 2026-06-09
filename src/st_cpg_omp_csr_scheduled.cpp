// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "st_csr_common.h"


#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

// =============================================================================
// Scheduled single-trait summary-stat BayesC with flat symmetric CSR LD
// =============================================================================
//
// This file provides a new exported function:
//   stblr_cpg_omp_csr_scheduled(...)
//
// It uses the updated 22-slot return layout, adding VLE/VLD traces, and changes the marker update
// loop to support approximate sparse scheduling:
//   - full sweeps every `full_sweep_every` iterations
//   - active markers updated every iteration
//   - candidate markers updated every iteration
//   - null markers updated only when scheduled
//   - optional LD-neighbor wakeup when marker effects change
//   - explicit Beta(pi_prior_a, pi_prior_b) prior for pi
//   - pi is updated only on full-sweep iterations
//   - VLE = sum_j ww_j b_j^2 / n and VLD = VG - VLE are returned
//
// Exact mode:
//   full_sweep_every = 1
//   null_skip_base   = 1
//
// Fast approximate mode similar to BED scheduled sampler:
//   full_sweep_every = 10
//   null_skip_base   = 50
//   null_skip_max    = 200
//   candidate_threshold = 1e-3
//   candidate_lifetime  = 20
//   wakeup_ld_neighbors = true
//   wakeup_diff_threshold = 0.0
// =============================================================================


// -----------------------------------------------------------------------------
// Marker update, returning p1 and diff for scheduling logic
// -----------------------------------------------------------------------------

struct STMarkerUpdateResult {
 double p1;
 double diff;
 int d_new;
};

inline STMarkerUpdateResult sampleBetaC_ST_csr_scheduled_one(
  int i,
  const std::vector<double>& pi,
  double vb,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const STLDCSR& ld,
  std::mt19937& gen
) {
 static thread_local std::uniform_real_distribution<double> runif(0.0, 1.0);
 static thread_local std::normal_distribution<double> norm01(0.0, 1.0);

 const arma::uword iu = static_cast<arma::uword>(i);
 const double wi = ww(iu);

 const double pi0 = std::max(pi[0], 1e-300);
 const double pi1 = std::max(pi[1], 1e-300);
 const double vei_safe = std::max(vei_i, 1e-300);

 const double score = r(iu) + wi * b(iu);
 const double denom = std::max(vei_safe + wi * vb, 1e-300);

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

 const int di = (runif(gen) < p1) ? 1 : 0;

 double b_new = 0.0;

 if (di == 1) {
  const double lhs = wi + vei_safe / vb;
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

 STMarkerUpdateResult out;
 out.p1 = p1;
 out.diff = diff;
 out.d_new = di;
 return out;
}

// -----------------------------------------------------------------------------
// Variance and pi updates
// -----------------------------------------------------------------------------

inline void sampleB_ST_csr_scheduled(
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

inline void sampleE_ST_csr_scheduled(
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
  throw std::runtime_error("sampleE_ST_csr_scheduled: invalid residual scale.");
 }

 std::chi_squared_distribution<double> rchisq(n + nue);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 const double ve_new = scale / chi2;

 if (!std::isfinite(ve_new) || ve_new <= 0.0) {
  throw std::runtime_error("sampleE_ST_csr_scheduled: sampled ve is invalid.");
 }

 ve = std::max(ve_new, 1e-12);
}

inline double computeG_ST_csr_scheduled(
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

inline double computeLE_ST_csr_scheduled(
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

inline void samplePi_ST_scheduled(
  const arma::Row<int>& d,
  std::vector<double>& pi,
  double pi_prior_a,
  double pi_prior_b,
  std::mt19937& gen
) {
 // pi[1] is inclusion probability; pi[0] is exclusion probability.
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

// -----------------------------------------------------------------------------
// Scheduling helpers
// -----------------------------------------------------------------------------

inline int adaptive_skip_length_csr_scheduled(
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

// =============================================================================
// Main exported scheduled CSR function
// =============================================================================

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_scheduled(
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
  int full_sweep_every,
  int null_skip_base,
  int null_skip_max,
  double candidate_threshold,
  int candidate_lifetime,
  bool skip_nulls_burnin_only,
  bool wakeup_ld_neighbors,
  double wakeup_diff_threshold,
  int wakeup_max_neighbors,
  double pi_prior_a,
  double pi_prior_b,
  int ncores,
  int seed
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nt must be positive.");
 const int m = static_cast<int>(wy[0].size());
 if (m <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: m must be positive.");
 if (nit <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nit must be positive.");
 if (nburn < 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nburn must be non-negative.");
 if (nthin <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nthin must be positive.");
 if (full_sweep_every < 0) throw std::runtime_error("full_sweep_every must be >= 0.");
 if (null_skip_base <= 0) throw std::runtime_error("null_skip_base must be positive.");
 if (null_skip_max < 0) throw std::runtime_error("null_skip_max must be >= 0.");
 if (!std::isfinite(candidate_threshold) || candidate_threshold < 0.0 || candidate_threshold > 1.0) {
  throw std::runtime_error("candidate_threshold must be in [0,1].");
 }
 if (candidate_lifetime < 0) throw std::runtime_error("candidate_lifetime must be >= 0.");
 if (!std::isfinite(wakeup_diff_threshold) || wakeup_diff_threshold < 0.0) {
  throw std::runtime_error("wakeup_diff_threshold must be >= 0.");
 }
 if (wakeup_max_neighbors < 0) throw std::runtime_error("wakeup_max_neighbors must be >= 0.");
 if (!std::isfinite(pi_prior_a) || pi_prior_a <= 0.0) {
  throw std::runtime_error("pi_prior_a must be finite and positive.");
 }
 if (!std::isfinite(pi_prior_b) || pi_prior_b <= 0.0) {
  throw std::runtime_error("pi_prior_b must be finite and positive.");
 }

 if ((int)ww.size() != nt || (int)b_init.size() != nt ||
     (int)yy.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_scheduled: inconsistent trait dimensions.");
 }

 if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_scheduled: priors must be nt x nt.");
 }

 if (pi.size() != 2) {
  throw std::runtime_error("stblr_cpg_omp_csr_scheduled: pi must have length 2, c(pi0, pi1).");
 }

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m || (int)ww[t].size() != m || (int)b_init[t].size() != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_scheduled: inconsistent marker dimensions.");
  }
 }

 if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_scheduled: B must be nt x nt.");
 }

 if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_scheduled: E must be nt x nt.");
 }

 if (use_r_init) {
  if (static_cast<int>(r_init.size()) != nt) {
   throw std::runtime_error("r_init must have length nt when use_r_init = true.");
  }
  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(r_init[t].size()) != m) {
    throw std::runtime_error("each r_init[t] must have length m.");
   }
  }
 }

 if (use_d_init) {
  if (static_cast<int>(d_init.size()) != nt) {
   throw std::runtime_error("d_init must have length nt when use_d_init = true.");
  }
  for (int t = 0; t < nt; ++t) {
   if (static_cast<int>(d_init[t].size()) != m) {
    throw std::runtime_error("each d_init[t] must have length m.");
   }
  }
 }

 // --------------------------------------------------------------------------
 // Convert inputs to Armadillo.
 // --------------------------------------------------------------------------

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);
 arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
 arma::vec yy_vec(nt, arma::fill::zeros);
 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  yy_vec(static_cast<arma::uword>(t)) = yy[t];

  for (int i = 0; i < m; ++i) {
   wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
   ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
   b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i))  = b_init[t][i];
  }

  if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("stblr_cpg_omp_csr_scheduled: priors must be nt x nt.");
  }

  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
   sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
  }
 }

 // --------------------------------------------------------------------------
 // Validate shared scaling and build LD.
 // --------------------------------------------------------------------------

 for (int t = 1; t < nt; ++t) {
  if (n[t] != n[0]) {
   throw std::runtime_error("stblr_cpg_omp_csr_scheduled: shared-LD scaling assumes equal n across traits.");
  }
 }

 for (int t = 1; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   const double w0 = ww_mat(0, static_cast<arma::uword>(i));
   const double wt = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   const double tol = 1e-8 * std::max(1.0, std::abs(w0));

   if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
    throw std::runtime_error("ww contains invalid value before LD pre-scaling.");
   }

   if (std::abs(w0 - wt) > tol) {
    throw std::runtime_error("ww differs across traits; pre-scaled shared ST LD is invalid.");
   }
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
 for (int i = 0; i < m; ++i) {
  const double wi = ww_mat(0, static_cast<arma::uword>(i));
  if (!std::isfinite(wi) || wi <= 0.0) {
   throw std::runtime_error("ww contains invalid value in trait 0.");
  }
  xx[static_cast<std::size_t>(i)] = wi;
 }

 //STLDCSR ld = read_and_build_st_ld_csr_scheduled(ld_prefix, m, xx);
 STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);

 // --------------------------------------------------------------------------
 // Marker order.
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

 std::sort(order.begin(), order.end(), [&](int a, int b) {
  return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
 });

 // --------------------------------------------------------------------------
 // Output storage.
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
 // Parallel over traits.
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
 << "STBLR scheduled CSR OpenMP requested threads = " << nthreads
 << ", omp_get_max_threads = " << omp_get_max_threads()
 << ", num procs = " << omp_get_num_procs()
 << ", pi_prior_a=" << pi_prior_a
 << ", pi_prior_b=" << pi_prior_b
 << "\n";
#endif

 Rcpp::Rcout
 << "Scheduled CSR sampler: full_sweep_every=" << full_sweep_every
 << ", null_skip_base=" << null_skip_base
 << ", null_skip_max=" << null_skip_max
 << ", candidate_threshold=" << candidate_threshold
 << ", candidate_lifetime=" << candidate_lifetime
 << ", skip_nulls_burnin_only=" << skip_nulls_burnin_only
 << ", wakeup_ld_neighbors=" << wakeup_ld_neighbors
 << ", wakeup_diff_threshold=" << wakeup_diff_threshold
 << ", wakeup_max_neighbors=" << wakeup_max_neighbors
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
   std::uniform_int_distribution<int> jitter_dist(0, std::max(0, null_skip_base - 1));

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));

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
     d_t(static_cast<arma::uword>(i)) = b_t(static_cast<arma::uword>(i)) != 0.0 ? 1 : 0;
    }
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) {
     r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
    }

    if (!r_t.is_finite()) {
     throw std::runtime_error("r_init contains NaN/Inf.");
    }
   } else {
    //rebuild_residual_st_csr_scheduled(m, wy_t, ww_t, b_t, r_t, ld);
    rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = computeG_ST_csr_scheduled(b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_ST_csr_scheduled(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   std::vector<double> pi_t = pi;
   if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
    throw std::runtime_error("invalid initial pi.");
   }
   {
    const double psum = pi_t[0] + pi_t[1];
    if (!std::isfinite(psum) || psum <= 0.0) throw std::runtime_error("invalid initial pi sum.");
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

   const int total_it = nit + nburn;
   const int bucket_count = total_it + std::max(2, null_skip_max > 0 ? null_skip_max : 4 * null_skip_base) + null_skip_base + 10;

   std::vector<std::vector<int>> scheduled(static_cast<std::size_t>(bucket_count));
   std::vector<int> scheduled_at(static_cast<std::size_t>(m), -1);
   std::vector<int> last_updated(static_cast<std::size_t>(m), -1);
   std::vector<unsigned char> is_candidate(static_cast<std::size_t>(m), 0u);
   std::vector<int> candidate_list;
   std::vector<unsigned char> in_candidate_list(static_cast<std::size_t>(m), 0u);
   std::vector<int> active_list;
   std::vector<unsigned char> in_active_list(static_cast<std::size_t>(m), 0u);
   std::vector<int> last_interesting(static_cast<std::size_t>(m), -1000000000);

   candidate_list.reserve(static_cast<std::size_t>(std::min(m, 10000)));
   active_list.reserve(static_cast<std::size_t>(std::min(m, 10000)));

   auto add_candidate = [&](int marker, int it) {
    if (marker < 0 || marker >= m) return;
    is_candidate[static_cast<std::size_t>(marker)] = 1u;
    last_interesting[static_cast<std::size_t>(marker)] = it;
    if (in_candidate_list[static_cast<std::size_t>(marker)] == 0u) {
     candidate_list.push_back(marker);
     in_candidate_list[static_cast<std::size_t>(marker)] = 1u;
    }
   };

   auto add_active = [&](int marker) {
    if (marker < 0 || marker >= m) return;
    if (in_active_list[static_cast<std::size_t>(marker)] == 0u) {
     active_list.push_back(marker);
     in_active_list[static_cast<std::size_t>(marker)] = 1u;
    }
   };

   auto schedule_marker = [&](int marker, int target_it) {
    if (marker < 0 || marker >= m) return;
    if (target_it >= bucket_count) target_it = bucket_count - 1;
    if (target_it < 0) target_it = 0;
    scheduled_at[static_cast<std::size_t>(marker)] = target_it;
    scheduled[static_cast<std::size_t>(target_it)].push_back(marker);
   };

   auto wakeup_neighbors = [&](int marker, int it) {
    if (!wakeup_ld_neighbors) return;
    const uint64_t start = ld.ptr[static_cast<std::size_t>(marker)];
    const uint64_t end   = ld.ptr[static_cast<std::size_t>(marker + 1)];
    int n_wake = 0;

    for (uint64_t p = start; p < end; ++p) {
     const int j = ld.idx[static_cast<std::size_t>(p)];
     add_candidate(j, it);
     scheduled_at[static_cast<std::size_t>(j)] = -1;
     ++n_wake;

     if (wakeup_max_neighbors > 0 && n_wake >= wakeup_max_neighbors) break;
    }
   };

   for (int i = 0; i < m; ++i) {
    if (d_t(static_cast<arma::uword>(i)) > 0) {
     add_active(i);
     add_candidate(i, 0);
    } else {
     const int skip = null_skip_base + (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
     schedule_marker(i, skip);
    }
   }

   auto update_one_marker = [&](int marker, int it) {
    if (marker < 0 || marker >= m) return;
    if (last_updated[static_cast<std::size_t>(marker)] == it) return;
    last_updated[static_cast<std::size_t>(marker)] = it;

    const STMarkerUpdateResult res = sampleBetaC_ST_csr_scheduled_one(
     marker,
     pi_t,
     vb_t,
     vei_t,
     ww_t,
     r_t,
     b_t,
     d_t,
     ld,
     gen_t
    );

    if (res.d_new > 0) {
     add_active(marker);
     add_candidate(marker, it);
     scheduled_at[static_cast<std::size_t>(marker)] = -1;
    } else if (res.p1 >= candidate_threshold) {
     add_candidate(marker, it);
     scheduled_at[static_cast<std::size_t>(marker)] = -1;
    } else {
     if (is_candidate[static_cast<std::size_t>(marker)] != 0u &&
         it - last_interesting[static_cast<std::size_t>(marker)] > candidate_lifetime) {
      is_candidate[static_cast<std::size_t>(marker)] = 0u;
     }

     if (is_candidate[static_cast<std::size_t>(marker)] == 0u) {
      const int skip = adaptive_skip_length_csr_scheduled(res.p1, null_skip_base, null_skip_max) +
       (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
      schedule_marker(marker, it + skip);
     }
    }

    if (std::abs(res.diff) > wakeup_diff_threshold) {
     wakeup_neighbors(marker, it);
    }
   };

   double nsamples_t = 0.0;

   for (int it = 0; it < total_it; ++it) {
    const bool skipping_allowed =
     null_skip_base > 1 &&
     (!skip_nulls_burnin_only || it < nburn);

    const bool full_sweep =
     !skipping_allowed ||
     full_sweep_every <= 0 ||
     ((it % full_sweep_every) == 0);

    if (full_sweep) {
     for (int isort = 0; isort < m; ++isort) {
      update_one_marker(order[static_cast<std::size_t>(isort)], it);
     }
    } else {
     for (int marker : active_list) {
      if (d_t(static_cast<arma::uword>(marker)) > 0) update_one_marker(marker, it);
     }

     for (int marker : candidate_list) {
      if (is_candidate[static_cast<std::size_t>(marker)] != 0u) update_one_marker(marker, it);
     }

     if (it < bucket_count) {
      const std::vector<int>& due = scheduled[static_cast<std::size_t>(it)];
      for (int marker : due) {
       if (scheduled_at[static_cast<std::size_t>(marker)] == it &&
           d_t(static_cast<arma::uword>(marker)) == 0 &&
           is_candidate[static_cast<std::size_t>(marker)] == 0u) {
        update_one_marker(marker, it);
       }
      }
     }
    }

    // Periodic compaction of stale active/candidate lists.
    if ((it + 1) % 50 == 0) {
     std::vector<int> active_new;
     active_new.reserve(active_list.size());
     std::fill(in_active_list.begin(), in_active_list.end(), 0u);

     for (int marker : active_list) {
      if (d_t(static_cast<arma::uword>(marker)) > 0 &&
          in_active_list[static_cast<std::size_t>(marker)] == 0u) {
       active_new.push_back(marker);
       in_active_list[static_cast<std::size_t>(marker)] = 1u;
      }
     }
     active_list.swap(active_new);

     std::vector<int> cand_new;
     cand_new.reserve(candidate_list.size());
     std::fill(in_candidate_list.begin(), in_candidate_list.end(), 0u);

     for (int marker : candidate_list) {
      if (is_candidate[static_cast<std::size_t>(marker)] != 0u &&
          in_candidate_list[static_cast<std::size_t>(marker)] == 0u) {
       cand_new.push_back(marker);
       in_candidate_list[static_cast<std::size_t>(marker)] = 1u;
      }
     }
     candidate_list.swap(cand_new);
    }

    if (updateB) {
     sampleB_ST_csr_scheduled(
      m,
      nub,
      vb_t,
      b_t,
      d_t,
      ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      gen_t
     );

     if (!std::isfinite(vb_t) || vb_t <= 0.0) {
      throw std::runtime_error("vb became invalid after sampleB. iter=" + std::to_string(it));
     }
    }

    if (updateE) {
     if (rebuild_r_before_updateE) {
      //rebuild_residual_st_csr_scheduled(m, wy_t, ww_t, b_t, r_t, ld);
      rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
     }

     sampleE_ST_csr_scheduled(
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

    if (updatePi && full_sweep) {
     samplePi_ST_scheduled(d_t, pi_t, pi_prior_a, pi_prior_b, gen_t);

     if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
      throw std::runtime_error("pi became invalid after samplePi. iter=" + std::to_string(it));
     }
    }

    vg_t = computeG_ST_csr_scheduled(b_t, wy_t, r_t, n[t]);
    vle_t = computeLE_ST_csr_scheduled(m, b_t, ww_t, n[t]);
    vld_t = vg_t - vle_t;

    if (!std::isfinite(vg_t)) {
     throw std::runtime_error("vg became NaN/Inf after computeG. iter=" + std::to_string(it));
    }
    if (!std::isfinite(vle_t)) {
     throw std::runtime_error("vle became NaN/Inf after computeLE. iter=" + std::to_string(it));
    }
    if (!std::isfinite(vld_t)) {
     throw std::runtime_error("vld became NaN/Inf after computeLE. iter=" + std::to_string(it));
    }

    vei_t = ve_t + adjE * vg_t;

    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error("adjusted residual variance vei became invalid. iter=" + std::to_string(it));
    }

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    pis_t(static_cast<arma::uword>(it)) = pi_t[1];
    vles_t(static_cast<arma::uword>(it)) = vle_t;
    vlds_t(static_cast<arma::uword>(it)) = vld_t;

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

   if (!bm_t.is_finite()) throw std::runtime_error("posterior mean bm contains NaN/Inf.");
   if (!dm_t.is_finite()) throw std::runtime_error("posterior mean dm contains NaN/Inf.");

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
   final_vle(static_cast<arma::uword>(t)) = vle_t;
   final_vld(static_cast<arma::uword>(t)) = vld_t;
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
     "stblr_cpg_omp_csr_scheduled failed for trait " +
      std::to_string(t) +
      ": " +
      errors[static_cast<std::size_t>(t)]
   );
  }
 }

 // --------------------------------------------------------------------------
 // Build 22-slot result.
 // --------------------------------------------------------------------------

 std::vector<std::vector<std::vector<double>>> result(22);
 for (int k = 0; k < 22; ++k) {
  result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);

  for (int k = 0; k <= 6; ++k) {
   result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(m));
  }

  for (int k = 7; k <= 9; ++k) {
   result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(nit + nburn));
  }

  for (int k = 10; k <= 15; ++k) {
   result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(nt));
  }

  result[16][ts].resize(2);
  result[17][ts].resize(2);
  result[18][ts].resize(4);
  result[19][ts].resize(2);
  result[20][ts].resize(static_cast<std::size_t>(nit + nburn)); // VLE
  result[21][ts].resize(static_cast<std::size_t>(nit + nburn)); // VLD = VG - VLE
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  const arma::uword tu = static_cast<arma::uword>(t);

  for (int i = 0; i < m; ++i) {
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
  const arma::uword tu = static_cast<arma::uword>(t);

  for (int it = 0; it < nit + nburn; ++it) {
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
  const arma::uword tu = static_cast<arma::uword>(t);

  result[16][ts][0] = 1.0 - final_pi(tu);
  result[16][ts][1] = final_pi(tu);

  double mean_pi = 0.0;
  int npi = 0;

  for (int it = nburn; it < nit + nburn; ++it) {
   mean_pi += pis_mat(tu, static_cast<arma::uword>(it));
   ++npi;
  }

  if (npi > 0) mean_pi /= static_cast<double>(npi);
  else mean_pi = final_pi(tu);

  result[17][ts][0] = 1.0 - mean_pi;
  result[17][ts][1] = mean_pi;

  for (int i = 0; i < 4; ++i) result[18][ts][static_cast<std::size_t>(i)] = 0.0;
  for (int i = 0; i < 2; ++i) result[19][ts][static_cast<std::size_t>(i)] = 0.0;
 }

 return result;
}

// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
//
// #include "cpg_samplers.h"
// #include "distributions.h"
// #include "st_csr_common.h"
//
//
// #include <algorithm>
// #include <cmath>
// #include <cstdio>
// #include <cstdint>
// #include <fstream>
// #include <limits>
// #include <numeric>
// #include <random>
// #include <stdexcept>
// #include <string>
// #include <vector>
// #include <unordered_map>
// #include <unordered_set>
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
// using namespace arma;
//
// // =============================================================================
// // Scheduled single-trait summary-stat BayesC with flat symmetric CSR LD
// // =============================================================================
// //
// // This file provides a new exported function:
// //   stblr_cpg_omp_csr_scheduled(...)
// //
// // It keeps the original 20-slot return layout, but changes the marker update
// // loop to support approximate sparse scheduling:
// //   - full sweeps every `full_sweep_every` iterations
// //   - active markers updated every iteration
// //   - candidate markers updated every iteration
// //   - null markers updated only when scheduled
// //   - optional LD-neighbor wakeup when marker effects change
// //
// // Exact mode:
// //   full_sweep_every = 1
// //   null_skip_base   = 1
// //
// // Fast approximate mode similar to BED scheduled sampler:
// //   full_sweep_every = 10
// //   null_skip_base   = 50
// //   null_skip_max    = 200
// //   candidate_threshold = 1e-3
// //   candidate_lifetime  = 20
// //   wakeup_ld_neighbors = true
// //   wakeup_diff_threshold = 0.0
// // =============================================================================
//
//
// // -----------------------------------------------------------------------------
// // Marker update, returning p1 and diff for scheduling logic
// // -----------------------------------------------------------------------------
//
// struct STMarkerUpdateResult {
//  double p1;
//  double diff;
//  int d_new;
// };
//
// inline STMarkerUpdateResult sampleBetaC_ST_csr_scheduled_one(
//   int i,
//   const std::vector<double>& pi,
//   double vb,
//   double vei_i,
//   const arma::rowvec& ww,
//   arma::rowvec& r,
//   arma::rowvec& b,
//   arma::Row<int>& d,
//   const STLDCSR& ld,
//   std::mt19937& gen
// ) {
//  static thread_local std::uniform_real_distribution<double> runif(0.0, 1.0);
//  static thread_local std::normal_distribution<double> norm01(0.0, 1.0);
//
//  const arma::uword iu = static_cast<arma::uword>(i);
//  const double wi = ww(iu);
//
//  const double pi0 = std::max(pi[0], 1e-300);
//  const double pi1 = std::max(pi[1], 1e-300);
//  const double vei_safe = std::max(vei_i, 1e-300);
//
//  const double score = r(iu) + wi * b(iu);
//  const double denom = std::max(vei_safe + wi * vb, 1e-300);
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
//  const int di = (runif(gen) < p1) ? 1 : 0;
//
//  double b_new = 0.0;
//
//  if (di == 1) {
//   const double lhs = wi + vei_safe / vb;
//   const double mean = score / lhs;
//   const double sd = std::sqrt(vei_safe / lhs);
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
//
//  STMarkerUpdateResult out;
//  out.p1 = p1;
//  out.diff = diff;
//  out.d_new = di;
//  return out;
// }
//
// // -----------------------------------------------------------------------------
// // Variance and pi updates
// // -----------------------------------------------------------------------------
//
// inline void sampleB_ST_csr_scheduled(
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
//   const arma::uword iu = static_cast<arma::uword>(i);
//   if (d(iu) > 0) {
//    ssb += b(iu) * b(iu);
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
// inline void sampleE_ST_csr_scheduled(
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
//   throw std::runtime_error("sampleE_ST_csr_scheduled: invalid residual scale.");
//  }
//
//  std::chi_squared_distribution<double> rchisq(n + nue);
//  const double chi2 = std::max(rchisq(gen), 1e-300);
//  const double ve_new = scale / chi2;
//
//  if (!std::isfinite(ve_new) || ve_new <= 0.0) {
//   throw std::runtime_error("sampleE_ST_csr_scheduled: sampled ve is invalid.");
//  }
//
//  ve = std::max(ve_new, 1e-12);
// }
//
// inline double computeG_ST_csr_scheduled(
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
// inline void samplePi_ST_scheduled(
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
// // -----------------------------------------------------------------------------
// // Scheduling helpers
// // -----------------------------------------------------------------------------
//
// inline int adaptive_skip_length_csr_scheduled(
//   double p1,
//   int null_skip_base,
//   int null_skip_max
// ) {
//  if (null_skip_base <= 1) return 1;
//
//  int skip = null_skip_base;
//
//  if (p1 < 1e-6) skip = 4 * null_skip_base;
//  else if (p1 < 1e-5) skip = 2 * null_skip_base;
//  else if (p1 < 1e-4) skip = null_skip_base;
//  else if (p1 < 1e-3) skip = std::max(1, null_skip_base / 2);
//  else skip = 1;
//
//  if (null_skip_max > 0) skip = std::min(skip, null_skip_max);
//  return std::max(1, skip);
// }
//
// // =============================================================================
// // Main exported scheduled CSR function
// // =============================================================================
//
// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_scheduled(
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
//   int full_sweep_every,
//   int null_skip_base,
//   int null_skip_max,
//   double candidate_threshold,
//   int candidate_lifetime,
//   bool skip_nulls_burnin_only,
//   bool wakeup_ld_neighbors,
//   double wakeup_diff_threshold,
//   int wakeup_max_neighbors,
//   int ncores,
//   int seed
// ) {
//  const int nt = static_cast<int>(wy.size());
//
//  if (nt <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nt must be positive.");
//  const int m = static_cast<int>(wy[0].size());
//  if (m <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: m must be positive.");
//  if (nit <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nit must be positive.");
//  if (nburn < 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nburn must be non-negative.");
//  if (nthin <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nthin must be positive.");
//  if (full_sweep_every < 0) throw std::runtime_error("full_sweep_every must be >= 0.");
//  if (null_skip_base <= 0) throw std::runtime_error("null_skip_base must be positive.");
//  if (null_skip_max < 0) throw std::runtime_error("null_skip_max must be >= 0.");
//  if (!std::isfinite(candidate_threshold) || candidate_threshold < 0.0 || candidate_threshold > 1.0) {
//   throw std::runtime_error("candidate_threshold must be in [0,1].");
//  }
//  if (candidate_lifetime < 0) throw std::runtime_error("candidate_lifetime must be >= 0.");
//  if (!std::isfinite(wakeup_diff_threshold) || wakeup_diff_threshold < 0.0) {
//   throw std::runtime_error("wakeup_diff_threshold must be >= 0.");
//  }
//  if (wakeup_max_neighbors < 0) throw std::runtime_error("wakeup_max_neighbors must be >= 0.");
//
//  if ((int)ww.size() != nt || (int)b_init.size() != nt ||
//      (int)yy.size() != nt || (int)n.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_scheduled: inconsistent trait dimensions.");
//  }
//
//  if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_scheduled: priors must be nt x nt.");
//  }
//
//  if (pi.size() != 2) {
//   throw std::runtime_error("stblr_cpg_omp_csr_scheduled: pi must have length 2, c(pi0, pi1).");
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   if ((int)wy[t].size() != m || (int)ww[t].size() != m || (int)b_init[t].size() != m) {
//    throw std::runtime_error("stblr_cpg_omp_csr_scheduled: inconsistent marker dimensions.");
//   }
//  }
//
//  if ((int)B.n_rows != nt || (int)B.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_scheduled: B must be nt x nt.");
//  }
//
//  if ((int)E.n_rows != nt || (int)E.n_cols != nt) {
//   throw std::runtime_error("stblr_cpg_omp_csr_scheduled: E must be nt x nt.");
//  }
//
//  if (use_r_init) {
//   if (static_cast<int>(r_init.size()) != nt) {
//    throw std::runtime_error("r_init must have length nt when use_r_init = true.");
//   }
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(r_init[t].size()) != m) {
//     throw std::runtime_error("each r_init[t] must have length m.");
//    }
//   }
//  }
//
//  if (use_d_init) {
//   if (static_cast<int>(d_init.size()) != nt) {
//    throw std::runtime_error("d_init must have length nt when use_d_init = true.");
//   }
//   for (int t = 0; t < nt; ++t) {
//    if (static_cast<int>(d_init[t].size()) != m) {
//     throw std::runtime_error("each d_init[t] must have length m.");
//    }
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Convert inputs to Armadillo.
//  // --------------------------------------------------------------------------
//
//  arma::mat wy_mat(nt, m, arma::fill::zeros);
//  arma::mat ww_mat(nt, m, arma::fill::zeros);
//  arma::mat b_mat(nt, m, arma::fill::zeros);
//  arma::mat r_mat(nt, m, arma::fill::zeros);
//  arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
//  arma::vec yy_vec(nt, arma::fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   yy_vec(static_cast<arma::uword>(t)) = yy[t];
//
//   for (int i = 0; i < m; ++i) {
//    wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
//    ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
//    b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i))  = b_init[t][i];
//   }
//
//   if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
//    throw std::runtime_error("stblr_cpg_omp_csr_scheduled: priors must be nt x nt.");
//   }
//
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
//    sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Validate shared scaling and build LD.
//  // --------------------------------------------------------------------------
//
//  for (int t = 1; t < nt; ++t) {
//   if (n[t] != n[0]) {
//    throw std::runtime_error("stblr_cpg_omp_csr_scheduled: shared-LD scaling assumes equal n across traits.");
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
//     throw std::runtime_error("ww contains invalid value before LD pre-scaling.");
//    }
//
//    if (std::abs(w0 - wt) > tol) {
//     throw std::runtime_error("ww differs across traits; pre-scaled shared ST LD is invalid.");
//    }
//   }
//  }
//
//  std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
//  for (int i = 0; i < m; ++i) {
//   const double wi = ww_mat(0, static_cast<arma::uword>(i));
//   if (!std::isfinite(wi) || wi <= 0.0) {
//    throw std::runtime_error("ww contains invalid value in trait 0.");
//   }
//   xx[static_cast<std::size_t>(i)] = wi;
//  }
//
//  //STLDCSR ld = read_and_build_st_ld_csr_scheduled(ld_prefix, m, xx);
//  STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
//
//  // --------------------------------------------------------------------------
//  // Marker order.
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
//  std::sort(order.begin(), order.end(), [&](int a, int b) {
//   return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
//  });
//
//  // --------------------------------------------------------------------------
//  // Output storage.
//  // --------------------------------------------------------------------------
//
//  arma::mat bm_mat(nt, m, arma::fill::zeros);
//  arma::mat dm_mat(nt, m, arma::fill::zeros);
//  arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
//  arma::vec final_vb(nt, arma::fill::zeros);
//  arma::vec final_vg(nt, arma::fill::zeros);
//  arma::vec final_ve(nt, arma::fill::zeros);
//  arma::vec final_pi(nt, arma::fill::zeros);
//  arma::vec nsamples_vec(nt, arma::fill::zeros);
//
//  // --------------------------------------------------------------------------
//  // Parallel over traits.
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
//  << "STBLR scheduled CSR OpenMP requested threads = " << nthreads
//  << ", omp_get_max_threads = " << omp_get_max_threads()
//  << ", num procs = " << omp_get_num_procs()
//  << "\n";
// #endif
//
//  Rcpp::Rcout
//  << "Scheduled CSR sampler: full_sweep_every=" << full_sweep_every
//  << ", null_skip_base=" << null_skip_base
//  << ", null_skip_max=" << null_skip_max
//  << ", candidate_threshold=" << candidate_threshold
//  << ", candidate_lifetime=" << candidate_lifetime
//  << ", skip_nulls_burnin_only=" << skip_nulls_burnin_only
//  << ", wakeup_ld_neighbors=" << wakeup_ld_neighbors
//  << ", wakeup_diff_threshold=" << wakeup_diff_threshold
//  << ", wakeup_max_neighbors=" << wakeup_max_neighbors
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
//    std::uniform_int_distribution<int> jitter_dist(0, std::max(0, null_skip_base - 1));
//
//    arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
//    arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));
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
//      d_t(static_cast<arma::uword>(i)) = b_t(static_cast<arma::uword>(i)) != 0.0 ? 1 : 0;
//     }
//    }
//
//    if (use_r_init) {
//     for (int i = 0; i < m; ++i) {
//      r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
//     }
//
//     if (!r_t.is_finite()) {
//      throw std::runtime_error("r_init contains NaN/Inf.");
//     }
//    } else {
//     //rebuild_residual_st_csr_scheduled(m, wy_t, ww_t, b_t, r_t, ld);
//     rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
//    }
//
//    double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
//    double vg_t = computeG_ST_csr_scheduled(b_t, wy_t, r_t, n[t]);
//    double vei_t = ve_t + adjE * vg_t;
//
//    std::vector<double> pi_t = pi;
//    if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//     throw std::runtime_error("invalid initial pi.");
//    }
//    {
//     const double psum = pi_t[0] + pi_t[1];
//     if (!std::isfinite(psum) || psum <= 0.0) throw std::runtime_error("invalid initial pi sum.");
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
//    const int total_it = nit + nburn;
//    const int bucket_count = total_it + std::max(2, null_skip_max > 0 ? null_skip_max : 4 * null_skip_base) + null_skip_base + 10;
//
//    std::vector<std::vector<int>> scheduled(static_cast<std::size_t>(bucket_count));
//    std::vector<int> scheduled_at(static_cast<std::size_t>(m), -1);
//    std::vector<int> last_updated(static_cast<std::size_t>(m), -1);
//    std::vector<unsigned char> is_candidate(static_cast<std::size_t>(m), 0u);
//    std::vector<int> candidate_list;
//    std::vector<unsigned char> in_candidate_list(static_cast<std::size_t>(m), 0u);
//    std::vector<int> active_list;
//    std::vector<unsigned char> in_active_list(static_cast<std::size_t>(m), 0u);
//    std::vector<int> last_interesting(static_cast<std::size_t>(m), -1000000000);
//
//    candidate_list.reserve(static_cast<std::size_t>(std::min(m, 10000)));
//    active_list.reserve(static_cast<std::size_t>(std::min(m, 10000)));
//
//    auto add_candidate = [&](int marker, int it) {
//     if (marker < 0 || marker >= m) return;
//     is_candidate[static_cast<std::size_t>(marker)] = 1u;
//     last_interesting[static_cast<std::size_t>(marker)] = it;
//     if (in_candidate_list[static_cast<std::size_t>(marker)] == 0u) {
//      candidate_list.push_back(marker);
//      in_candidate_list[static_cast<std::size_t>(marker)] = 1u;
//     }
//    };
//
//    auto add_active = [&](int marker) {
//     if (marker < 0 || marker >= m) return;
//     if (in_active_list[static_cast<std::size_t>(marker)] == 0u) {
//      active_list.push_back(marker);
//      in_active_list[static_cast<std::size_t>(marker)] = 1u;
//     }
//    };
//
//    auto schedule_marker = [&](int marker, int target_it) {
//     if (marker < 0 || marker >= m) return;
//     if (target_it >= bucket_count) target_it = bucket_count - 1;
//     if (target_it < 0) target_it = 0;
//     scheduled_at[static_cast<std::size_t>(marker)] = target_it;
//     scheduled[static_cast<std::size_t>(target_it)].push_back(marker);
//    };
//
//    auto wakeup_neighbors = [&](int marker, int it) {
//     if (!wakeup_ld_neighbors) return;
//     const uint64_t start = ld.ptr[static_cast<std::size_t>(marker)];
//     const uint64_t end   = ld.ptr[static_cast<std::size_t>(marker + 1)];
//     int n_wake = 0;
//
//     for (uint64_t p = start; p < end; ++p) {
//      const int j = ld.idx[static_cast<std::size_t>(p)];
//      add_candidate(j, it);
//      scheduled_at[static_cast<std::size_t>(j)] = -1;
//      ++n_wake;
//
//      if (wakeup_max_neighbors > 0 && n_wake >= wakeup_max_neighbors) break;
//     }
//    };
//
//    for (int i = 0; i < m; ++i) {
//     if (d_t(static_cast<arma::uword>(i)) > 0) {
//      add_active(i);
//      add_candidate(i, 0);
//     } else {
//      const int skip = null_skip_base + (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
//      schedule_marker(i, skip);
//     }
//    }
//
//    auto update_one_marker = [&](int marker, int it) {
//     if (marker < 0 || marker >= m) return;
//     if (last_updated[static_cast<std::size_t>(marker)] == it) return;
//     last_updated[static_cast<std::size_t>(marker)] = it;
//
//     const STMarkerUpdateResult res = sampleBetaC_ST_csr_scheduled_one(
//      marker,
//      pi_t,
//      vb_t,
//      vei_t,
//      ww_t,
//      r_t,
//      b_t,
//      d_t,
//      ld,
//      gen_t
//     );
//
//     if (res.d_new > 0) {
//      add_active(marker);
//      add_candidate(marker, it);
//      scheduled_at[static_cast<std::size_t>(marker)] = -1;
//     } else if (res.p1 >= candidate_threshold) {
//      add_candidate(marker, it);
//      scheduled_at[static_cast<std::size_t>(marker)] = -1;
//     } else {
//      if (is_candidate[static_cast<std::size_t>(marker)] != 0u &&
//          it - last_interesting[static_cast<std::size_t>(marker)] > candidate_lifetime) {
//       is_candidate[static_cast<std::size_t>(marker)] = 0u;
//      }
//
//      if (is_candidate[static_cast<std::size_t>(marker)] == 0u) {
//       const int skip = adaptive_skip_length_csr_scheduled(res.p1, null_skip_base, null_skip_max) +
//        (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
//       schedule_marker(marker, it + skip);
//      }
//     }
//
//     if (std::abs(res.diff) > wakeup_diff_threshold) {
//      wakeup_neighbors(marker, it);
//     }
//    };
//
//    double nsamples_t = 0.0;
//
//    for (int it = 0; it < total_it; ++it) {
//     const bool skipping_allowed =
//      null_skip_base > 1 &&
//      (!skip_nulls_burnin_only || it < nburn);
//
//     const bool full_sweep =
//      !skipping_allowed ||
//      full_sweep_every <= 0 ||
//      ((it % full_sweep_every) == 0);
//
//     if (full_sweep) {
//      for (int isort = 0; isort < m; ++isort) {
//       update_one_marker(order[static_cast<std::size_t>(isort)], it);
//      }
//     } else {
//      for (int marker : active_list) {
//       if (d_t(static_cast<arma::uword>(marker)) > 0) update_one_marker(marker, it);
//      }
//
//      for (int marker : candidate_list) {
//       if (is_candidate[static_cast<std::size_t>(marker)] != 0u) update_one_marker(marker, it);
//      }
//
//      if (it < bucket_count) {
//       const std::vector<int>& due = scheduled[static_cast<std::size_t>(it)];
//       for (int marker : due) {
//        if (scheduled_at[static_cast<std::size_t>(marker)] == it &&
//            d_t(static_cast<arma::uword>(marker)) == 0 &&
//            is_candidate[static_cast<std::size_t>(marker)] == 0u) {
//         update_one_marker(marker, it);
//        }
//       }
//      }
//     }
//
//     // Periodic compaction of stale active/candidate lists.
//     if ((it + 1) % 50 == 0) {
//      std::vector<int> active_new;
//      active_new.reserve(active_list.size());
//      std::fill(in_active_list.begin(), in_active_list.end(), 0u);
//
//      for (int marker : active_list) {
//       if (d_t(static_cast<arma::uword>(marker)) > 0 &&
//           in_active_list[static_cast<std::size_t>(marker)] == 0u) {
//        active_new.push_back(marker);
//        in_active_list[static_cast<std::size_t>(marker)] = 1u;
//       }
//      }
//      active_list.swap(active_new);
//
//      std::vector<int> cand_new;
//      cand_new.reserve(candidate_list.size());
//      std::fill(in_candidate_list.begin(), in_candidate_list.end(), 0u);
//
//      for (int marker : candidate_list) {
//       if (is_candidate[static_cast<std::size_t>(marker)] != 0u &&
//           in_candidate_list[static_cast<std::size_t>(marker)] == 0u) {
//        cand_new.push_back(marker);
//        in_candidate_list[static_cast<std::size_t>(marker)] = 1u;
//       }
//      }
//      candidate_list.swap(cand_new);
//     }
//
//     if (updateB) {
//      sampleB_ST_csr_scheduled(
//       m,
//       nub,
//       vb_t,
//       b_t,
//       d_t,
//       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
//       gen_t
//      );
//
//      if (!std::isfinite(vb_t) || vb_t <= 0.0) {
//       throw std::runtime_error("vb became invalid after sampleB. iter=" + std::to_string(it));
//      }
//     }
//
//     if (updateE) {
//      if (rebuild_r_before_updateE) {
//       //rebuild_residual_st_csr_scheduled(m, wy_t, ww_t, b_t, r_t, ld);
//       rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
//      }
//
//      sampleE_ST_csr_scheduled(
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
//      samplePi_ST_scheduled(d_t, pi_t, gen_t);
//
//      if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
//       throw std::runtime_error("pi became invalid after samplePi. iter=" + std::to_string(it));
//      }
//     }
//
//     vg_t = computeG_ST_csr_scheduled(b_t, wy_t, r_t, n[t]);
//
//     if (!std::isfinite(vg_t)) {
//      throw std::runtime_error("vg became NaN/Inf after computeG. iter=" + std::to_string(it));
//     }
//
//     vei_t = ve_t + adjE * vg_t;
//
//     if (!std::isfinite(vei_t) || vei_t <= 0.0) {
//      throw std::runtime_error("adjusted residual variance vei became invalid. iter=" + std::to_string(it));
//     }
//
//     vbs_t(static_cast<arma::uword>(it)) = vb_t;
//     ves_t(static_cast<arma::uword>(it)) = ve_t;
//     vgs_t(static_cast<arma::uword>(it)) = vg_t;
//     pis_t(static_cast<arma::uword>(it)) = pi_t[1];
//
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
//    bm_t /= nsamples_t;
//    dm_t /= nsamples_t;
//
//    if (!bm_t.is_finite()) throw std::runtime_error("posterior mean bm contains NaN/Inf.");
//    if (!dm_t.is_finite()) throw std::runtime_error("posterior mean dm contains NaN/Inf.");
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
//      "stblr_cpg_omp_csr_scheduled failed for trait " +
//       std::to_string(t) +
//       ": " +
//       errors[static_cast<std::size_t>(t)]
//    );
//   }
//  }
//
//  // --------------------------------------------------------------------------
//  // Build 20-slot result.
//  // --------------------------------------------------------------------------
//
//  std::vector<std::vector<std::vector<double>>> result(20);
//  for (int k = 0; k < 20; ++k) {
//   result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//
//   for (int k = 0; k <= 6; ++k) {
//    result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(m));
//   }
//
//   for (int k = 7; k <= 9; ++k) {
//    result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(nit + nburn));
//   }
//
//   for (int k = 10; k <= 15; ++k) {
//    result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(nt));
//   }
//
//   result[16][ts].resize(2);
//   result[17][ts].resize(2);
//   result[18][ts].resize(4);
//   result[19][ts].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//   const arma::uword tu = static_cast<arma::uword>(t);
//
//   for (int i = 0; i < m; ++i) {
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
//   const arma::uword tu = static_cast<arma::uword>(t);
//
//   for (int it = 0; it < nit + nburn; ++it) {
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
//    result[13][t1s][t2s] = 0.0;
//    result[14][t1s][t2s] = 0.0;
//    result[15][t1s][t2s] = 0.0;
//   }
//
//   result[10][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
//   result[11][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
//   result[12][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
//   result[13][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
//   result[14][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
//   result[15][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   const std::size_t ts = static_cast<std::size_t>(t);
//   const arma::uword tu = static_cast<arma::uword>(t);
//
//   result[16][ts][0] = 1.0 - final_pi(tu);
//   result[16][ts][1] = final_pi(tu);
//
//   double mean_pi = 0.0;
//   int npi = 0;
//
//   for (int it = nburn; it < nit + nburn; ++it) {
//    mean_pi += pis_mat(tu, static_cast<arma::uword>(it));
//    ++npi;
//   }
//
//   if (npi > 0) mean_pi /= static_cast<double>(npi);
//   else mean_pi = final_pi(tu);
//
//   result[17][ts][0] = 1.0 - mean_pi;
//   result[17][ts][1] = mean_pi;
//
//   for (int i = 0; i < 4; ++i) result[18][ts][static_cast<std::size_t>(i)] = 0.0;
//   for (int i = 0; i < 2; ++i) result[19][ts][static_cast<std::size_t>(i)] = 0.0;
//  }
//
//  return result;
// }
