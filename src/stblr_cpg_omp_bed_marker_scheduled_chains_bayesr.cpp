// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "packed_bed.h"
#include "st_bed_bayesr_common.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
// Packed BED marker-wise ST BayesR sampler with scheduled null-marker updates
// and multiple independent chains.
//
// BayesR extension of stblr_cpg_omp_bed_marker_scheduled_chains():
//   - K-component mixture prior: c[0]=0 (null), c[1..K-1] are non-null scalars
//   - d[j] is a component index in {0,...,K-1} instead of 0/1
//   - pi is a K-vector sampled from Dirichlet(alpha + counts)
//   - sampleB uses b_j^2 / c[d_j] as each marker's ssb contribution
//   - Additional result slot [22]: per-component PIP matrix, flat K*m per trait
//   - Standard dm is P(component > 0); posterior mean component index is
//     preserved in result slot [29].
//
// All BED reading, packing, scheduling logic, chain parallelism, and R/C++
// signatures for BED/packed matrix routines are unchanged from the BayesC version.
// =============================================================================

struct ChainResultBayesR {
 arma::rowvec bm;
 arma::rowvec dm;
 arma::rowvec component_mean;
 arma::rowvec b;
 arma::rowvec d_as_double;
 arma::rowvec vbs;
 arma::rowvec vgs;
 arma::rowvec ves;
 arma::mat    pip_k;           // K x m: per-component PIP, averaged over post-burnin
 arma::rowvec vles;
 arma::rowvec vlds;
 double final_vb = 0.0;
 double final_vg = 0.0;
 double final_ve = 0.0;
 std::vector<double> final_pi; // length K
 double final_vle = 0.0;
 double final_vld = 0.0;
 double log_cpo = NA_REAL;
 double mean_log_cpo = NA_REAL;
 std::vector<double> mean_pi;  // length K
 double nsamples = 0.0;
 double seconds = 0.0;
 int failed = 0;
 std::string error;
};

// FastPackedBedMatrix: identical definition to st_cpg_omp_individual_scheduled_chains.cpp.
// Static functions below have internal linkage so there is no ODR conflict.
// BayesR marker sampler: K-component mixture.
// Returns 1 - prob[0] (non-null probability) for the scheduling heuristic.
static inline double sample_marker_bayesr(
  const FastPackedBedMatrixBR& G,
  int marker,
  const MarkerMapBayesR& map,
  const std::vector<double>& pi,
  const std::vector<double>& c,
  double vb,
  double vei,
  arma::vec& e,
  double& b_j,
  int& d_j,
  std::mt19937& gen
) {
 static thread_local std::uniform_real_distribution<double> runif(0.0, 1.0);
 static thread_local std::normal_distribution<double> norm01(0.0, 1.0);

 const int K = static_cast<int>(c.size());
 const double xx = map.xx;
 const double vei_safe = std::max(vei, 1e-300);

 const double xte = br_dot_residual(G, marker, map, e.memptr());
 const double score = xte + xx * b_j;

 // K log-posteriors
 std::vector<double> logp(static_cast<std::size_t>(K));
 logp[0] = std::log(std::max(pi[0], 1e-300));
 for (int k = 1; k < K; ++k) {
  const double vbk = vb * c[static_cast<std::size_t>(k)];
  const double denom = std::max(vei_safe + xx * vbk, 1e-300);
  logp[static_cast<std::size_t>(k)] =
   std::log(std::max(pi[static_cast<std::size_t>(k)], 1e-300)) +
   0.5 * std::log(vei_safe / denom) +
   0.5 * score * score * vbk / (vei_safe * denom);
 }

 // Softmax: subtract max for numerical stability, then normalize
 double logp_max = logp[0];
 for (int k = 1; k < K; ++k)
  if (logp[static_cast<std::size_t>(k)] > logp_max)
   logp_max = logp[static_cast<std::size_t>(k)];

 std::vector<double> prob(static_cast<std::size_t>(K));
 double psum = 0.0;
 for (int k = 0; k < K; ++k) {
  prob[static_cast<std::size_t>(k)] = std::exp(logp[static_cast<std::size_t>(k)] - logp_max);
  psum += prob[static_cast<std::size_t>(k)];
 }
 for (int k = 0; k < K; ++k)
  prob[static_cast<std::size_t>(k)] /= psum;

 // Sample component via inverse CDF; default to K-1 if u falls at the tail
 int d_new = K - 1;
 const double u = runif(gen);
 double cumsum = 0.0;
 for (int k = 0; k < K - 1; ++k) {
  cumsum += prob[static_cast<std::size_t>(k)];
  if (u < cumsum) { d_new = k; break; }
 }

 double b_new = 0.0;
 if (d_new > 0) {
  const double vbk = vb * c[static_cast<std::size_t>(d_new)];
  const double lhs = xx + vei_safe / vbk;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b_j;
 if (diff != 0.0)
  br_update_residual(G, marker, map, e.memptr(), diff);

 b_j = b_new;
 d_j = d_new;
 return 1.0 - prob[0]; // non-null probability, mirrors p1 used in scheduling
}

// BayesR vb sampler: ssb contribution for marker j is b_j^2 / c[d_j].
// BayesR pi sampler: Dirichlet(alpha[k] + n_k) via K independent Gamma samples.
static inline void samplePi_bayesr(
  const arma::Row<int>& d,
  std::vector<double>& pi,
  const std::vector<double>& alpha,
  std::mt19937& gen
) {
 const int K = static_cast<int>(pi.size());

 std::vector<double> counts(static_cast<std::size_t>(K), 0.0);
 for (arma::uword i = 0; i < d.n_elem; ++i) {
  const int di = d(i);
  if (di >= 0 && di < K)
   counts[static_cast<std::size_t>(di)] += 1.0;
 }

 double gsum = 0.0;
 std::vector<double> g(static_cast<std::size_t>(K));
 for (int k = 0; k < K; ++k) {
  const double shape = alpha[static_cast<std::size_t>(k)] + counts[static_cast<std::size_t>(k)];
  std::gamma_distribution<double> rgamma(shape, 1.0);
  g[static_cast<std::size_t>(k)] = std::max(rgamma(gen), 1e-300);
  gsum += g[static_cast<std::size_t>(k)];
 }

 for (int k = 0; k < K; ++k)
  pi[static_cast<std::size_t>(k)] = g[static_cast<std::size_t>(k)] / gsum;
}

static inline int br_adaptive_skip(double p_nonnull, int null_skip_base, int null_skip_max) {
 if (null_skip_base <= 1) return 1;

 int skip = null_skip_base;
 if (p_nonnull < 1e-6) skip = 4 * null_skip_base;
 else if (p_nonnull < 1e-5) skip = 2 * null_skip_base;
 else if (p_nonnull < 1e-4) skip = null_skip_base;
 else if (p_nonnull < 1e-3) skip = std::max(1, null_skip_base / 2);
 else skip = 1;

 if (null_skip_max > 0) skip = std::min(skip, null_skip_max);
 return std::max(1, skip);
}

static ChainResultBayesR run_one_bayesr_chain(
  const FastPackedBedMatrixBR& G,
  const std::vector<MarkerMapBayesR>& marker_maps,
  const std::vector<int>& marker_order,
  const arma::mat& y_mat,
  const std::vector<std::vector<double>>& b_init,
  const arma::mat& B,
  const arma::mat& E,
  const arma::mat& ssb_prior_mat,
  const arma::mat& sse_prior_mat,
  const std::vector<double>& pi_init,
  const std::vector<double>& c,
  const std::vector<double>& alpha,
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
  int t,
  int chain,
  int seed,
  int progress_every
) {
#ifdef _OPENMP
 const double wall_start = omp_get_wtime();
#else
 const double wall_start = 0.0;
#endif

 const int m = G.m;
 const int K = static_cast<int>(c.size());
 ChainResultBayesR out;
 out.bm = arma::rowvec(m, arma::fill::zeros);
 out.dm = arma::rowvec(m, arma::fill::zeros);
 out.component_mean = arma::rowvec(m, arma::fill::zeros);
 out.b = arma::rowvec(m, arma::fill::zeros);
 out.d_as_double = arma::rowvec(m, arma::fill::zeros);
 out.vbs = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.vgs = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.ves = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.pip_k = arma::mat(static_cast<arma::uword>(K), static_cast<arma::uword>(m), arma::fill::zeros);
 out.vles = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.vlds = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.final_pi.assign(static_cast<std::size_t>(K), 0.0);
 out.mean_pi.assign(static_cast<std::size_t>(K), 0.0);

 try {
  const unsigned int chain_seed = static_cast<unsigned int>(
   seed + 1000003 * (t + 1) + 9176 * (chain + 1)
  );

  std::mt19937 gen_t(chain_seed);
  std::uniform_int_distribution<int> jitter_dist(0, std::max(0, null_skip_base - 1));

  arma::vec y_t = y_mat.col(static_cast<arma::uword>(t));
  arma::rowvec b_t(m, arma::fill::zeros);
  arma::Row<int> d_t(m, arma::fill::zeros);

  for (int j = 0; j < m; ++j) {
   b_t(static_cast<arma::uword>(j)) = b_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(j)];
   d_t(static_cast<arma::uword>(j)) = b_t(static_cast<arma::uword>(j)) != 0.0 ? 1 : 0;
  }

  arma::vec xb_t = br_xb(G, marker_maps, marker_order, b_t);
  arma::vec e_t = y_t - xb_t;

  double vb_t = B(t, t);
  double ve_t = E(t, t);
  double vg_t = br_computeG(y_t, e_t);
  double vei_t = ve_t + adjE * vg_t;

  std::vector<double> pi_t = pi_init;
  {
   double psum = 0.0;
   for (int k = 0; k < K; ++k) psum += pi_t[static_cast<std::size_t>(k)];
   if (!std::isfinite(psum) || psum <= 0.0)
    throw std::runtime_error("invalid initial pi: sum <= 0.");
   for (int k = 0; k < K; ++k) {
    pi_t[static_cast<std::size_t>(k)] /= psum;
    if (!std::isfinite(pi_t[static_cast<std::size_t>(k)]) ||
        pi_t[static_cast<std::size_t>(k)] < 0.0)
     throw std::runtime_error("invalid initial pi after normalization.");
   }
  }

  arma::rowvec bm_t(m, arma::fill::zeros);
  arma::rowvec dm_t(m, arma::fill::zeros);
  arma::rowvec component_mean_t(m, arma::fill::zeros);
  arma::mat    pip_k_t(static_cast<arma::uword>(K), static_cast<arma::uword>(m), arma::fill::zeros);
  arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
  arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);
  arma::vec log_inv_cpo_t(G.n, arma::fill::value(-std::numeric_limits<double>::infinity()));

  std::vector<double> mean_pi_acc(static_cast<std::size_t>(K), 0.0);
  int npi_acc = 0;

  const int total_it = nit + nburn;
  const int bucket_count = total_it +
   std::max(2, null_skip_max > 0 ? null_skip_max : 4 * null_skip_base) +
   null_skip_base + 10;

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

  auto add_candidate = [&](int marker) {
   is_candidate[static_cast<std::size_t>(marker)] = 1u;
   if (in_candidate_list[static_cast<std::size_t>(marker)] == 0u) {
    candidate_list.push_back(marker);
    in_candidate_list[static_cast<std::size_t>(marker)] = 1u;
   }
  };

  auto add_active = [&](int marker) {
   if (in_active_list[static_cast<std::size_t>(marker)] == 0u) {
    active_list.push_back(marker);
    in_active_list[static_cast<std::size_t>(marker)] = 1u;
   }
  };

  auto schedule_marker = [&](int marker, int target_it) {
   if (target_it >= bucket_count) target_it = bucket_count - 1;
   if (target_it < 0) target_it = 0;
   scheduled_at[static_cast<std::size_t>(marker)] = target_it;
   scheduled[static_cast<std::size_t>(target_it)].push_back(marker);
  };

  for (int j = 0; j < m; ++j) {
   if (d_t(static_cast<arma::uword>(j)) > 0) {
    add_active(j);
    add_candidate(j);
    last_interesting[static_cast<std::size_t>(j)] = 0;
   } else {
    const int skip = null_skip_base + (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
    schedule_marker(j, skip);
   }
  }

  auto update_one_marker = [&](int marker, int it) {
   if (marker < 0 || marker >= m) return;
   if (last_updated[static_cast<std::size_t>(marker)] == it) return;
   last_updated[static_cast<std::size_t>(marker)] = it;

   const arma::uword ju = static_cast<arma::uword>(marker);
   double bj = b_t(ju);
   int dj = d_t(ju);

   const double p_nonnull = sample_marker_bayesr(
    G,
    marker,
    marker_maps[static_cast<std::size_t>(marker)],
    pi_t,
    c,
    vb_t,
    vei_t,
    e_t,
    bj,
    dj,
    gen_t
   );

   b_t(ju) = bj;
   d_t(ju) = dj;

   if (dj > 0) {
    add_active(marker);
    add_candidate(marker);
    last_interesting[static_cast<std::size_t>(marker)] = it;
    scheduled_at[static_cast<std::size_t>(marker)] = -1;
    return;
   }

   if (p_nonnull >= candidate_threshold) {
    add_candidate(marker);
    last_interesting[static_cast<std::size_t>(marker)] = it;
    scheduled_at[static_cast<std::size_t>(marker)] = -1;
    return;
   }

   if (is_candidate[static_cast<std::size_t>(marker)] != 0u &&
       it - last_interesting[static_cast<std::size_t>(marker)] > candidate_lifetime) {
    is_candidate[static_cast<std::size_t>(marker)] = 0u;
   }

   if (is_candidate[static_cast<std::size_t>(marker)] == 0u) {
    const int skip = br_adaptive_skip(p_nonnull, null_skip_base, null_skip_max) +
     (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
    schedule_marker(marker, it + skip);
   }
  };

  double nsamples_t = 0.0;

  for (int it = 0; it < total_it; ++it) {
   if (progress_every > 0 &&
       (it == 0 || ((it + 1) % progress_every == 0) || it + 1 == total_it)) {
    double n_included_progress = 0.0;
    for (arma::uword jj = 0; jj < d_t.n_elem; ++jj)
     if (d_t(jj) > 0) n_included_progress += 1.0;

    const double vle_progress = br_computeLE(b_t, marker_maps, G.n);
    const double vld_progress = vg_t - vle_progress;
    const double pi_nonnull = 1.0 - pi_t[0];

#ifdef _OPENMP
#pragma omp critical
#endif
{
 Rcpp::Rcout
 << "progress chain " << chain
 << ", trait " << t
 << ": iter " << (it + 1)
 << "/" << total_it
 << ", vb=" << vb_t
 << ", ve=" << ve_t
 << ", vg=" << vg_t
 << ", vle=" << vle_progress
 << ", vld=" << vld_progress
 << ", vei=" << vei_t
 << ", pi_nonnull=" << pi_nonnull
 << ", K=" << K
 << ", n_included=" << n_included_progress
 << ", active=" << active_list.size()
 << ", candidates=" << candidate_list.size()
 << "\n";
}
   }

   const bool skipping_allowed =
    null_skip_base > 1 &&
    (!skip_nulls_burnin_only || it < nburn);

   const bool full_sweep =
    !skipping_allowed ||
    full_sweep_every <= 0 ||
    ((it % full_sweep_every) == 0);

   if (full_sweep) {
    for (int marker : marker_order)
     update_one_marker(marker, it);
   } else {
    for (int marker : active_list)
     if (d_t(static_cast<arma::uword>(marker)) > 0) update_one_marker(marker, it);

    for (int marker : candidate_list)
     if (is_candidate[static_cast<std::size_t>(marker)] != 0u) update_one_marker(marker, it);

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

   if (updateE && rebuild_every > 0 && ((it + 1) % rebuild_every == 0)) {
    xb_t = br_xb(G, marker_maps, marker_order, b_t);
    e_t = y_t - xb_t;
   }

   if (updateB) {
    sampleB_bayesr(m, c, nub, vb_t, b_t, d_t, ssb_prior_mat(t, t), gen_t);
    if (!std::isfinite(vb_t) || vb_t <= 0.0) throw std::runtime_error("invalid vb.");
   }

   if (updateE) {
    sampleE_bayesr(nue, ve_t, e_t, sse_prior_mat(t, t), gen_t);
    if (!std::isfinite(ve_t) || ve_t <= 0.0) throw std::runtime_error("invalid ve.");
   }

   if (updatePi && full_sweep) {
    samplePi_bayesr(d_t, pi_t, alpha, gen_t);
    double psum_check = 0.0;
    for (int k = 0; k < K; ++k) psum_check += pi_t[static_cast<std::size_t>(k)];
    if (!std::isfinite(psum_check) || psum_check <= 0.0)
     throw std::runtime_error("invalid pi after sampling.");
   }

   vg_t = br_computeG(y_t, e_t);
   const double vle_t = br_computeLE(b_t, marker_maps, G.n);
   const double vld_t = vg_t - vle_t;
   vei_t = ve_t + adjE * vg_t;
   if (!std::isfinite(vg_t)) throw std::runtime_error("invalid vg.");
   if (!std::isfinite(vle_t)) throw std::runtime_error("invalid vle.");
   if (!std::isfinite(vld_t)) throw std::runtime_error("invalid vld.");
   if (!std::isfinite(vei_t) || vei_t <= 0.0) throw std::runtime_error("invalid adjusted residual variance.");

   vbs_t(static_cast<arma::uword>(it)) = vb_t;
   ves_t(static_cast<arma::uword>(it)) = ve_t;
   vgs_t(static_cast<arma::uword>(it)) = vg_t;
   vles_t(static_cast<arma::uword>(it)) = vle_t;
   vlds_t(static_cast<arma::uword>(it)) = vld_t;

   // Accumulate pi for mean_pi (post-burnin only, no thinning, mirrors original)
   if (it >= nburn) {
    for (int k = 0; k < K; ++k)
     mean_pi_acc[static_cast<std::size_t>(k)] += pi_t[static_cast<std::size_t>(k)];
    ++npi_acc;
   }

   if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
    nsamples_t += 1.0;
    br_update_log_inv_cpo(e_t, ve_t, log_inv_cpo_t);

    for (int j = 0; j < m; ++j) {
     const arma::uword ju = static_cast<arma::uword>(j);
     bm_t(ju) += b_t(ju);
     dm_t(ju) += d_t(ju) > 0 ? 1.0 : 0.0;
     component_mean_t(ju) += static_cast<double>(d_t(ju));
     pip_k_t(static_cast<arma::uword>(d_t(ju)), ju) += 1.0;
    }
   }
  }

  if (nsamples_t <= 0.0) nsamples_t = 1.0;

  const double log_cpo_t = br_compute_total_log_cpo(log_inv_cpo_t, nsamples_t);
  const double mean_log_cpo_t = log_cpo_t / static_cast<double>(G.n);

  bm_t /= nsamples_t;
  dm_t /= nsamples_t;
  component_mean_t /= nsamples_t;
  pip_k_t /= nsamples_t;

  out.bm = bm_t;
  out.dm = dm_t;
  out.component_mean = component_mean_t;
  out.b = b_t;
  for (int j = 0; j < m; ++j)
   out.d_as_double(static_cast<arma::uword>(j)) = static_cast<double>(d_t(static_cast<arma::uword>(j)));
  out.vbs = vbs_t;
  out.vgs = vgs_t;
  out.ves = ves_t;
  out.pip_k = pip_k_t;
  out.vles = vles_t;
  out.vlds = vlds_t;
  out.final_vb = vb_t;
  out.final_ve = ve_t;
  out.final_vg = vg_t;
  out.final_vle = br_computeLE(b_t, marker_maps, G.n);
  out.final_vld = out.final_vg - out.final_vle;
  out.log_cpo = log_cpo_t;
  out.mean_log_cpo = mean_log_cpo_t;
  out.nsamples = nsamples_t;

  for (int k = 0; k < K; ++k) {
   out.final_pi[static_cast<std::size_t>(k)] = pi_t[static_cast<std::size_t>(k)];
   out.mean_pi[static_cast<std::size_t>(k)] =
    npi_acc > 0
     ? mean_pi_acc[static_cast<std::size_t>(k)] / static_cast<double>(npi_acc)
     : pi_t[static_cast<std::size_t>(k)];
  }

#ifdef _OPENMP
  out.seconds = omp_get_wtime() - wall_start;
#else
  out.seconds = 0.0;
#endif

 } catch (const std::exception& ex) {
  out.failed = 1;
  out.error = ex.what();
#ifdef _OPENMP
  out.seconds = omp_get_wtime() - wall_start;
#endif
 } catch (...) {
  out.failed = 1;
  out.error = "unknown error";
#ifdef _OPENMP
  out.seconds = omp_get_wtime() - wall_start;
#endif
 }

 return out;
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_bed_marker_scheduled_chains_bayesr(
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
  std::vector<double> pi,      // length K: initial mixture probabilities
  std::vector<double> c,       // length K: c[0]=0, c[1..K-1] are non-null variance scalars
  std::vector<double> alpha,   // length K: Dirichlet prior concentrations
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
  int nchains,
  int ncores,
  int seed
) {
 if (nit <= 0 || nburn < 0)
  throw std::runtime_error("stblr_cpg_omp_bed_marker_scheduled_chains_bayesr: nit must be positive and nburn non-negative.");
 if (nthin <= 0) throw std::runtime_error("nthin must be positive.");
 if (rebuild_every < 0) throw std::runtime_error("rebuild_every must be >= 0.");
 if (full_sweep_every < 0) throw std::runtime_error("full_sweep_every must be >= 0.");
 if (null_skip_base <= 0) throw std::runtime_error("null_skip_base must be positive.");
 if (null_skip_max < 0) throw std::runtime_error("null_skip_max must be >= 0.");
 if (nchains <= 0) throw std::runtime_error("nchains must be positive.");
 if (read_block_size <= 0) throw std::runtime_error("read_block_size must be positive.");
 if (progress_every < 0) throw std::runtime_error("progress_every must be >= 0.");
 if (!std::isfinite(candidate_threshold) || candidate_threshold < 0.0 || candidate_threshold > 1.0)
  throw std::runtime_error("candidate_threshold must be in [0,1].");
 if (candidate_lifetime < 0) throw std::runtime_error("candidate_lifetime must be >= 0.");

 const int K = static_cast<int>(c.size());
 if (K < 2)
  throw std::runtime_error("c must have at least 2 elements (null + one non-null component).");
 if (static_cast<int>(pi.size()) != K)
  throw std::runtime_error("pi must have the same length as c (K).");
 if (static_cast<int>(alpha.size()) != K)
  throw std::runtime_error("alpha must have the same length as c (K).");
 if (c[0] != 0.0)
  throw std::runtime_error("c[0] must be 0 (null component).");
 for (int k = 1; k < K; ++k) {
  if (!std::isfinite(c[static_cast<std::size_t>(k)]) || c[static_cast<std::size_t>(k)] <= 0.0)
   throw std::runtime_error("c[k] for k>=1 must be finite and positive.");
 }
 for (int k = 0; k < K; ++k) {
  if (!std::isfinite(alpha[static_cast<std::size_t>(k)]) || alpha[static_cast<std::size_t>(k)] <= 0.0)
   throw std::runtime_error("alpha[k] must be finite and positive.");
 }

 std::vector<std::string> bed_files_cpp = br_copy_bed_files(bed_files);
 std::vector<std::vector<int>> cls_by_file = br_copy_int_list(cls);
 std::vector<int> rows0 = br_copy_rows0(rows, n);

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

 FastPackedBedMatrixBR G = br_read_bed_blocked(
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
 if (y.nrow() != n_used)
  throw std::runtime_error("y rows must equal the number of samples used after rows filtering.");
 if (static_cast<int>(sets.size()) != m)
  throw std::runtime_error("sets length must equal number of markers.");
 if (static_cast<int>(b_init.size()) != nt)
  throw std::runtime_error("b_init must have length nt.");
 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(b_init[static_cast<std::size_t>(t)].size()) != m)
   throw std::runtime_error("each b_init[t] must have length m.");
 }
 if (static_cast<int>(B.n_rows) != nt || static_cast<int>(B.n_cols) != nt)
  throw std::runtime_error("B must be nt x nt.");
 if (static_cast<int>(E.n_rows) != nt || static_cast<int>(E.n_cols) != nt)
  throw std::runtime_error("E must be nt x nt.");
 if (static_cast<int>(ssb_prior.size()) != nt || static_cast<int>(sse_prior.size()) != nt)
  throw std::runtime_error("prior lists must have length nt.");

 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(ssb_prior[t].size()) != nt || static_cast<int>(sse_prior[t].size()) != nt)
   throw std::runtime_error("priors must be nt x nt.");
  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(t, t2) = ssb_prior[t][t2];
   sse_prior_mat(t, t2) = sse_prior[t][t2];
  }
 }

 std::vector<double> af_cpp = br_flatten_af_list(af);
 const bool af_computed = af_cpp.empty();
 if (af_computed) af_cpp = br_compute_af(G);
 if (static_cast<int>(af_cpp.size()) != m)
  throw std::runtime_error("af must have one value per marker.");

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
 << "Building BayesR scheduled packed BED chains: n=" << n_used
 << ", m=" << m
 << ", K=" << K
 << ", nt=" << nt
 << ", nchains=" << nchains
 << ", jobs=" << njobs
 << ", scale=" << scale
 << ", af_computed=" << af_computed
 << ", read_block_size=" << read_block_size
 << ", progress_every=" << progress_every
 << "\n";

 std::vector<MarkerMapBayesR> marker_maps = br_build_marker_maps(G, af_cpp, scale, nthreads);
 std::vector<int> marker_order = br_make_marker_order(sets, m);

 Rcpp::Rcout
 << "BayesR scheduled chains: full_sweep_every=" << full_sweep_every
 << ", null_skip_base=" << null_skip_base
 << ", null_skip_max=" << null_skip_max
 << ", candidate_threshold=" << candidate_threshold
 << ", candidate_lifetime=" << candidate_lifetime
 << ", skip_nulls_burnin_only=" << skip_nulls_burnin_only
 << ", return_wy=" << return_wy
 << ", return_r=" << return_r
 << "\n";

#ifdef _OPENMP
 Rcpp::Rcout
 << "BayesR scheduled packed BED chains OpenMP threads = " << nthreads
 << ", max threads = " << omp_get_max_threads()
 << ", num procs = " << omp_get_num_procs()
 << "\n";
#else
 Rcpp::Rcout << "BayesR scheduled packed BED chains compiled without OpenMP; using one thread.\n";
#endif

 arma::mat y_mat(n_used, nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t)
  for (int i = 0; i < n_used; ++i) y_mat(i, t) = y(i, t);

 std::vector<ChainResultBayesR> job_results(static_cast<std::size_t>(njobs));

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int job = 0; job < njobs; ++job) {
  const int ch = job / nt;
  const int t  = job % nt;

  job_results[static_cast<std::size_t>(job)] = run_one_bayesr_chain(
   G,
   marker_maps,
   marker_order,
   y_mat,
   b_init,
   B,
   E,
   ssb_prior_mat,
   sse_prior_mat,
   pi,
   c,
   alpha,
   nub,
   nue,
   updateB,
   updateE,
   updatePi,
   adjE,
   nit,
   nburn,
   nthin,
   rebuild_every,
   full_sweep_every,
   null_skip_base,
   null_skip_max,
   candidate_threshold,
   candidate_lifetime,
   skip_nulls_burnin_only,
   t,
   ch,
   seed,
   progress_every
  );
 }

 int failed_total = 0;
 for (int job = 0; job < njobs; ++job) {
  const ChainResultBayesR& r = job_results[static_cast<std::size_t>(job)];
  const int ch = job / nt;
  const int t  = job % nt;
  Rcpp::Rcout
  << "chain " << ch
  << ", trait " << t
  << ", failed=" << r.failed
  << ", seconds=" << r.seconds
  << ", nsamples=" << r.nsamples
  << ", log_cpo=" << r.log_cpo
  << ", mean_log_cpo=" << r.mean_log_cpo
  << "\n";
  if (r.failed) {
   ++failed_total;
   Rcpp::Rcout << "  error: " << r.error << "\n";
  }
 }

 if (failed_total > 0) {
  for (int job = 0; job < njobs; ++job) {
   const ChainResultBayesR& r = job_results[static_cast<std::size_t>(job)];
   if (r.failed) {
    const int ch = job / nt;
    const int t  = job % nt;
    throw std::runtime_error(
      "stblr_cpg_omp_bed_marker_scheduled_chains_bayesr failed for chain " +
       std::to_string(ch) + ", trait " + std::to_string(t) + ": " + r.error
    );
   }
  }
 }

 // Aggregate across chains
 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);
 arma::mat component_mean_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat d_mat_double(nt, m, arma::fill::zeros);
 arma::mat bm_sd_mat(nt, m, arma::fill::zeros);
 arma::mat dm_sd_mat(nt, m, arma::fill::zeros);
 arma::mat bm_min_mat(nt, m, arma::fill::zeros);
 arma::mat dm_min_mat(nt, m, arma::fill::zeros);
 arma::mat bm_max_mat(nt, m, arma::fill::zeros);
 arma::mat dm_max_mat(nt, m, arma::fill::zeros);
 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vles_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vlds_mat(nt, nit + nburn, arma::fill::zeros);
 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_vle(nt, arma::fill::zeros);
 arma::vec final_vld(nt, arma::fill::zeros);
 arma::mat final_pi_mat(static_cast<arma::uword>(K), static_cast<arma::uword>(nt), arma::fill::zeros);
 arma::mat mean_pi_mat(static_cast<arma::uword>(K), static_cast<arma::uword>(nt), arma::fill::zeros);
 arma::vec mean_total_log_cpo(nt, arma::fill::zeros);
 arma::vec mean_log_cpo(nt, arma::fill::zeros);
 arma::vec mean_nsamples(nt, arma::fill::zeros);
 arma::vec mean_seconds(nt, arma::fill::zeros);
 arma::vec max_seconds(nt, arma::fill::zeros);

 // pip_k aggregation: one K x m matrix per trait
 std::vector<arma::mat> pip_k_mat(
  static_cast<std::size_t>(nt),
  arma::mat(static_cast<arma::uword>(K), static_cast<arma::uword>(m), arma::fill::zeros)
 );

 for (int ch = 0; ch < nchains; ++ch) {
  for (int t = 0; t < nt; ++t) {
   const int job = ch * nt + t;
   const ChainResultBayesR& r = job_results[static_cast<std::size_t>(job)];
   const arma::uword tu = static_cast<arma::uword>(t);

   bm_mat.row(tu) += r.bm;
   dm_mat.row(tu) += r.dm;
   component_mean_mat.row(tu) += r.component_mean;
   b_mat.row(tu) += r.b;
   d_mat_double.row(tu) += r.d_as_double;
   vbs_mat.row(tu) += r.vbs;
   vgs_mat.row(tu) += r.vgs;
   ves_mat.row(tu) += r.ves;
   vles_mat.row(tu) += r.vles;
   vlds_mat.row(tu) += r.vlds;
   pip_k_mat[static_cast<std::size_t>(t)] += r.pip_k;
   final_vb(tu) += r.final_vb;
   final_vg(tu) += r.final_vg;
   final_ve(tu) += r.final_ve;
   final_vle(tu) += r.final_vle;
   final_vld(tu) += r.final_vld;
   mean_total_log_cpo(tu) += r.log_cpo;
   mean_log_cpo(tu) += r.mean_log_cpo;
   mean_nsamples(tu) += r.nsamples;
   mean_seconds(tu) += r.seconds;
   max_seconds(tu) = std::max(max_seconds(tu), r.seconds);

   for (int k = 0; k < K; ++k) {
    final_pi_mat(static_cast<arma::uword>(k), tu) += r.final_pi[static_cast<std::size_t>(k)];
    mean_pi_mat(static_cast<arma::uword>(k), tu)  += r.mean_pi[static_cast<std::size_t>(k)];
   }
  }
 }

 const double inv_chains = 1.0 / static_cast<double>(nchains);
 bm_mat *= inv_chains;
 dm_mat *= inv_chains;
 component_mean_mat *= inv_chains;
 b_mat *= inv_chains;
 d_mat_double *= inv_chains;
 vbs_mat *= inv_chains;
 vgs_mat *= inv_chains;
 ves_mat *= inv_chains;
 vles_mat *= inv_chains;
 vlds_mat *= inv_chains;
 for (auto& mat : pip_k_mat) mat *= inv_chains;
 final_vb *= inv_chains;
 final_vg *= inv_chains;
 final_ve *= inv_chains;
 final_vle *= inv_chains;
 final_vld *= inv_chains;
 final_pi_mat *= inv_chains;
 mean_pi_mat *= inv_chains;
 mean_total_log_cpo *= inv_chains;
 mean_log_cpo *= inv_chains;
 mean_nsamples *= inv_chains;
 mean_seconds *= inv_chains;

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  bm_min_mat.row(tu).fill(std::numeric_limits<double>::infinity());
  dm_min_mat.row(tu).fill(std::numeric_limits<double>::infinity());
  bm_max_mat.row(tu).fill(-std::numeric_limits<double>::infinity());
  dm_max_mat.row(tu).fill(-std::numeric_limits<double>::infinity());

  for (int ch = 0; ch < nchains; ++ch) {
   const int job = ch * nt + t;
   const ChainResultBayesR& r = job_results[static_cast<std::size_t>(job)];

   for (int j = 0; j < m; ++j) {
    const arma::uword ju = static_cast<arma::uword>(j);
    bm_min_mat(tu, ju) = std::min(bm_min_mat(tu, ju), r.bm(ju));
    dm_min_mat(tu, ju) = std::min(dm_min_mat(tu, ju), r.dm(ju));
    bm_max_mat(tu, ju) = std::max(bm_max_mat(tu, ju), r.bm(ju));
    dm_max_mat(tu, ju) = std::max(dm_max_mat(tu, ju), r.dm(ju));
   }

   if (nchains > 1) {
    const arma::rowvec bm_diff = r.bm - bm_mat.row(tu);
    const arma::rowvec dm_diff = r.dm - dm_mat.row(tu);
    bm_sd_mat.row(tu) += bm_diff % bm_diff;
    dm_sd_mat.row(tu) += dm_diff % dm_diff;
   }
  }

  if (nchains > 1) {
   bm_sd_mat.row(tu) = arma::sqrt(bm_sd_mat.row(tu) / static_cast<double>(nchains - 1));
   dm_sd_mat.row(tu) = arma::sqrt(dm_sd_mat.row(tu) / static_cast<double>(nchains - 1));
  }
 }

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);

 if (return_wy || return_r) {
  for (int t = 0; t < nt; ++t) {
   arma::vec y_t = y_mat.col(static_cast<arma::uword>(t));
   arma::rowvec b_t = b_mat.row(static_cast<arma::uword>(t));
   arma::vec xb_t = br_xb(G, marker_maps, marker_order, b_t);
   arma::vec e_t = y_t - xb_t;

   for (int j = 0; j < m; ++j) {
    if (return_wy)
     wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(j)) =
      br_dot_residual(G, j, marker_maps[static_cast<std::size_t>(j)], y_t.memptr());
    if (return_r)
     r_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(j)) =
      br_dot_residual(G, j, marker_maps[static_cast<std::size_t>(j)], e_t.memptr());
   }
  }
 }

 // --------------------------------------------------------------------------
 // Build named raw schema v1 (same schema as the migrated CSR backends;
 // values below are numerically identical to the previous positional
 // result[0..29] slots).
 // --------------------------------------------------------------------------

 const int n_trace = nit + nburn;

 auto marker_matrix = [&](const arma::mat& x) {
  Rcpp::NumericMatrix out(m, nt);
  for (int t = 0; t < nt; ++t) {
   for (int i = 0; i < m; ++i) {
    out(i, t) = x(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   }
  }
  return out;
 };

 auto trace_matrix = [&](const arma::mat& x) {
  Rcpp::NumericMatrix out(n_trace, nt);
  for (int t = 0; t < nt; ++t) {
   for (int it = 0; it < n_trace; ++it) {
    out(it, t) = x(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
   }
  }
  return out;
 };

 auto diagonal_matrix = [&](const arma::vec& x) {
  Rcpp::NumericMatrix out(nt, nt);
  for (int t = 0; t < nt; ++t) {
   out(t, t) = x(static_cast<arma::uword>(t));
  }
  return out;
 };

 Rcpp::CharacterVector component_names(K);
 for (int k = 0; k < K; ++k) {
  component_names[k] = "component_" + std::to_string(k);
 }

 arma::vec mixture_var_vec(static_cast<arma::uword>(K));
 for (int k = 0; k < K; ++k) mixture_var_vec(static_cast<arma::uword>(k)) = c[static_cast<std::size_t>(k)];

 arma::mat final_pi(nt, K, arma::fill::zeros);
 arma::mat mean_pi(nt, K, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  for (int k = 0; k < K; ++k) {
   final_pi(static_cast<arma::uword>(t), static_cast<arma::uword>(k)) =
    final_pi_mat(static_cast<arma::uword>(k), static_cast<arma::uword>(t));
   mean_pi(static_cast<arma::uword>(t), static_cast<arma::uword>(k)) =
    mean_pi_mat(static_cast<arma::uword>(k), static_cast<arma::uword>(t));
  }
 }

 Rcpp::List comp_prob_out(nt);
 for (int t = 0; t < nt; ++t) {
  comp_prob_out[t] = arma::mat(pip_k_mat[static_cast<std::size_t>(t)].t());
 }

 Rcpp::NumericVector nsamples_out(nt);
 Rcpp::IntegerVector n_used_out(nt);
 Rcpp::NumericVector log_cpo_out(nt);
 Rcpp::NumericVector mean_log_cpo_out(nt);
 Rcpp::NumericVector seconds_mean_out(nt);
 Rcpp::NumericVector seconds_max_out(nt);

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  nsamples_out[t] = mean_nsamples(tu);
  n_used_out[t] = n_used;
  log_cpo_out[t] = mean_total_log_cpo(tu);
  mean_log_cpo_out[t] = mean_log_cpo(tu);
  seconds_mean_out[t] = mean_seconds(tu);
  seconds_max_out[t] = max_seconds(tu);
 }

 Rcpp::List marker = Rcpp::List::create(
  Rcpp::Named("bm") = marker_matrix(bm_mat),
  Rcpp::Named("dm") = marker_matrix(dm_mat),
  Rcpp::Named("wy") = marker_matrix(wy_mat),
  Rcpp::Named("r") = marker_matrix(r_mat),
  Rcpp::Named("b") = marker_matrix(b_mat),
  Rcpp::Named("state") = marker_matrix(d_mat_double),
  Rcpp::Named("bm_sd") = marker_matrix(bm_sd_mat),
  Rcpp::Named("bm_min") = marker_matrix(bm_min_mat),
  Rcpp::Named("bm_max") = marker_matrix(bm_max_mat),
  Rcpp::Named("dm_sd") = marker_matrix(dm_sd_mat),
  Rcpp::Named("dm_min") = marker_matrix(dm_min_mat),
  Rcpp::Named("dm_max") = marker_matrix(dm_max_mat)
 );

 Rcpp::List trace = Rcpp::List::create(
  Rcpp::Named("vbs") = trace_matrix(vbs_mat),
  Rcpp::Named("vgs") = trace_matrix(vgs_mat),
  Rcpp::Named("ves") = trace_matrix(ves_mat),
  Rcpp::Named("vle") = trace_matrix(vles_mat),
  Rcpp::Named("vld") = trace_matrix(vlds_mat),
  Rcpp::Named("pis") = R_NilValue
 );

 Rcpp::List variance = Rcpp::List::create(
  Rcpp::Named("covb") = diagonal_matrix(final_vb),
  Rcpp::Named("covg") = diagonal_matrix(final_vg),
  Rcpp::Named("cove") = diagonal_matrix(final_ve),
  Rcpp::Named("vb") = diagonal_matrix(final_vb),
  Rcpp::Named("vg") = diagonal_matrix(final_vg),
  Rcpp::Named("ve") = diagonal_matrix(final_ve)
 );

 Rcpp::List diagnostics = Rcpp::List::create(
  Rcpp::Named("nsamples") = nsamples_out,
  Rcpp::Named("n_used") = n_used_out,
  Rcpp::Named("log_cpo") = log_cpo_out,
  Rcpp::Named("mean_log_cpo") = mean_log_cpo_out,
  Rcpp::Named("seconds_mean") = seconds_mean_out,
  Rcpp::Named("seconds_max") = seconds_max_out,
  Rcpp::Named("ld_swap") = R_NilValue
 );

 Rcpp::List selection = Rcpp::List::create(
  Rcpp::Named("enabled") = false,
  Rcpp::Named("fixed") = false,
  Rcpp::Named("scale") = "standardized_genotype_effect",
  Rcpp::Named("trace") = R_NilValue,
  Rcpp::Named("mean") = R_NilValue,
  Rcpp::Named("sd") = R_NilValue,
  Rcpp::Named("min") = R_NilValue,
  Rcpp::Named("max") = R_NilValue,
  Rcpp::Named("acceptance") = R_NilValue
 );

 Rcpp::List raw = Rcpp::List::create(
  Rcpp::Named("schema") = Rcpp::List::create(
   Rcpp::Named("class") = "stblr_raw",
   Rcpp::Named("version") = 1
  ),
  Rcpp::Named("meta") = Rcpp::List::create(
   Rcpp::Named("model") = "bayesr",
   Rcpp::Named("backend") = "bed_scheduled_chains_bayesr",
   Rcpp::Named("data_level") = "individual",
   Rcpp::Named("prior_type") = "component",
   Rcpp::Named("m") = m,
   Rcpp::Named("nt") = nt,
   Rcpp::Named("n_trace") = n_trace,
   Rcpp::Named("nit") = nit,
   Rcpp::Named("nburn") = nburn,
   Rcpp::Named("nthin") = nthin,
   Rcpp::Named("nchains") = nchains,
   Rcpp::Named("keep_chains") = false,
   Rcpp::Named("n_components") = K,
   Rcpp::Named("n_annotations") = 0,
   Rcpp::Named("n_groups") = 0
  ),
  Rcpp::Named("marker") = marker,
  Rcpp::Named("trace") = trace,
  Rcpp::Named("variance") = variance,
  Rcpp::Named("pi") = Rcpp::List::create(
   Rcpp::Named("final") = final_pi,
   Rcpp::Named("mean") = mean_pi,
   Rcpp::Named("names") = component_names
  ),
  Rcpp::Named("diagnostics") = diagnostics,
  Rcpp::Named("chains") = R_NilValue,
  Rcpp::Named("prior") = Rcpp::List::create(),
  Rcpp::Named("group") = Rcpp::List::create(),
  Rcpp::Named("annotation") = Rcpp::List::create(),
  Rcpp::Named("component") = Rcpp::List::create(
   Rcpp::Named("names") = component_names,
   Rcpp::Named("mixture_var") = mixture_var_vec,
   Rcpp::Named("prob") = comp_prob_out,
   Rcpp::Named("ncomp") = R_NilValue,
   Rcpp::Named("dm_component_mean") = marker_matrix(component_mean_mat)
  ),
  Rcpp::Named("selection") = selection
 );
 raw.attr("class") = Rcpp::CharacterVector::create("stblr_raw_v1", "stblr_raw", "list");
 return raw;
}
