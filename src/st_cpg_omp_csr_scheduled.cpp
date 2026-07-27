// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "st_csr_common.h"
#include "st_chain_utils.h"
#include "blr_scheduled_execution_types.h"
#include "blr_csr_scheduled_bayesc_types.h"


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
  sblr::core::ScheduledChainRng& rng
) {
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

 const int di = (rng.uniform(rng.engine) < p1) ? 1 : 0;

 double b_new = 0.0;

 if (di == 1) {
  const double lhs = wi + vei_safe / vb;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * rng.normal(rng.engine);
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

#include "blr_csr_scheduled_bayesc_core_impl.h"

struct STScheduledBayesCBindingMetadata {
 const arma::mat& wy;
 const std::vector<int>& sample_sizes;
 int nit,nburn,nthin,nchains;
 bool keep_chains;
};

static Rcpp::List stblr_csr_scheduled_bayesc_result_to_raw(
 const sblr::core::CsrScheduledBayesCExecutionResult& result,
 const STScheduledBayesCBindingMetadata& binding
) {
 const int m=result.marker_count,nt=result.trait_count;
 const int n_trace=binding.nit+binding.nburn;
 auto marker_matrix=[&](const arma::mat& x) {
  Rcpp::NumericMatrix out(m,nt);
  for(int t=0;t<nt;++t) for(int i=0;i<m;++i)
   out(i,t)=x(static_cast<arma::uword>(t),static_cast<arma::uword>(i));
  return out;
 };
 auto trace_matrix=[&](const arma::mat& x) {
  Rcpp::NumericMatrix out(n_trace,nt);
  for(int t=0;t<nt;++t) for(int it=0;it<n_trace;++it)
   out(it,t)=x(static_cast<arma::uword>(t),static_cast<arma::uword>(it));
  return out;
 };
 auto diagonal_matrix=[&](const arma::vec& x) {
  Rcpp::NumericMatrix out(nt,nt);
  for(int t=0;t<nt;++t) out(t,t)=x(static_cast<arma::uword>(t));
  return out;
 };
 Rcpp::NumericMatrix pi_final(nt,2),pi_mean(nt,2);
 for(int t=0;t<nt;++t) {
  const arma::uword tu=static_cast<arma::uword>(t);
  pi_final(t,0)=1.0-result.final_pi(tu); pi_final(t,1)=result.final_pi(tu);
  const double mean_pi=result.mean_pi(tu);
  pi_mean(t,0)=1.0-mean_pi; pi_mean(t,1)=mean_pi;
 }
 Rcpp::NumericVector nsamples(nt),seconds_mean(nt),seconds_max(nt);
 Rcpp::IntegerVector n_used(nt);
 for(int t=0;t<nt;++t) {
  nsamples[t]=result.nsamples(static_cast<arma::uword>(t));
  n_used[t]=binding.sample_sizes[static_cast<std::size_t>(t)];
  double sec_sum=0.0,sec_max=0.0;
  for(int chain=0;chain<binding.nchains;++chain) {
   const double sec=result.task_seconds[static_cast<std::size_t>(t*binding.nchains+chain)];
   sec_sum+=sec; sec_max=std::max(sec_max,sec);
  }
  seconds_mean[t]=sec_sum/static_cast<double>(binding.nchains); seconds_max[t]=sec_max;
 }
 Rcpp::List marker=Rcpp::List::create(
  Rcpp::Named("bm")=marker_matrix(result.bm),Rcpp::Named("dm")=marker_matrix(result.dm),
  Rcpp::Named("wy")=marker_matrix(binding.wy),Rcpp::Named("r")=marker_matrix(result.r),
  Rcpp::Named("b")=marker_matrix(result.b),Rcpp::Named("state")=marker_matrix(result.state));
 if(binding.nchains>1) {
  marker["bm_sd"]=marker_matrix(result.bm_sd); marker["bm_min"]=marker_matrix(result.bm_min);
  marker["bm_max"]=marker_matrix(result.bm_max); marker["dm_sd"]=marker_matrix(result.dm_sd);
  marker["dm_min"]=marker_matrix(result.dm_min); marker["dm_max"]=marker_matrix(result.dm_max);
 }
 Rcpp::List trace=Rcpp::List::create(
  Rcpp::Named("vbs")=trace_matrix(result.vbs),Rcpp::Named("vgs")=trace_matrix(result.vgs),
  Rcpp::Named("ves")=trace_matrix(result.ves),Rcpp::Named("vle")=trace_matrix(result.vle),
  Rcpp::Named("vld")=trace_matrix(result.vld),Rcpp::Named("pis")=trace_matrix(result.pis));
 Rcpp::List variance=Rcpp::List::create(
  Rcpp::Named("covb")=diagonal_matrix(result.final_vb),
  Rcpp::Named("covg")=diagonal_matrix(result.final_vg),
  Rcpp::Named("cove")=diagonal_matrix(result.final_ve),
  Rcpp::Named("vb")=diagonal_matrix(result.final_vb),
  Rcpp::Named("vg")=diagonal_matrix(result.final_vg),
  Rcpp::Named("ve")=diagonal_matrix(result.final_ve));
 Rcpp::List diagnostics=Rcpp::List::create(
  Rcpp::Named("nsamples")=nsamples,Rcpp::Named("n_used")=n_used,
  Rcpp::Named("log_cpo")=Rcpp::NumericVector(nt),
  Rcpp::Named("mean_log_cpo")=Rcpp::NumericVector(nt),
  Rcpp::Named("seconds_mean")=seconds_mean,Rcpp::Named("seconds_max")=seconds_max,
  Rcpp::Named("ld_swap")=R_NilValue);
 Rcpp::List selection=Rcpp::List::create(
  Rcpp::Named("enabled")=false,Rcpp::Named("fixed")=false,
  Rcpp::Named("scale")="standardized_genotype_effect",Rcpp::Named("trace")=R_NilValue,
  Rcpp::Named("mean")=R_NilValue,Rcpp::Named("sd")=R_NilValue,
  Rcpp::Named("min")=R_NilValue,Rcpp::Named("max")=R_NilValue,
  Rcpp::Named("acceptance")=R_NilValue);
 Rcpp::List chains=R_NilValue;
 if (binding.keep_chains) {
  chains=Rcpp::List(nt);
  for (int t=0;t<nt;++t) {
   Rcpp::List trait_chains(binding.nchains);
   for (int chain=0;chain<binding.nchains;++chain) {
    const int task=t*binding.nchains+chain;
    const arma::uword task_u=static_cast<arma::uword>(task);
    trait_chains[chain]=Rcpp::List::create(
     Rcpp::Named("marker")=Rcpp::List::create(
      Rcpp::Named("bm")=result.task_bm.row(task_u).t(),
      Rcpp::Named("dm")=result.task_dm.row(task_u).t(),
      Rcpp::Named("state")=result.task_state.row(task_u).t()),
     Rcpp::Named("trace")=Rcpp::List::create(
      Rcpp::Named("vbs")=result.task_vbs.row(task_u).t(),
      Rcpp::Named("vgs")=result.task_vgs.row(task_u).t(),
      Rcpp::Named("ves")=result.task_ves.row(task_u).t(),
      Rcpp::Named("vle")=result.task_vle.row(task_u).t(),
      Rcpp::Named("vld")=result.task_vld.row(task_u).t(),
      Rcpp::Named("pis")=result.task_pis.row(task_u).t()),
     Rcpp::Named("pi")=Rcpp::List::create(
      Rcpp::Named("final")=Rcpp::NumericVector::create(
       1.0-result.task_final_pi(task_u),result.task_final_pi(task_u)),
      Rcpp::Named("mean")=Rcpp::NumericVector::create(
       1.0-result.task_mean_pi(task_u),result.task_mean_pi(task_u))),
     Rcpp::Named("retained_draw_count")=result.task_nsamples(task_u),
     Rcpp::Named("diagnostics")=Rcpp::List::create(
      Rcpp::Named("seconds")=result.task_seconds[static_cast<std::size_t>(task)]));
   }
   chains[t]=trait_chains;
  }
 }
 Rcpp::List raw=Rcpp::List::create(
  Rcpp::Named("schema")=Rcpp::List::create(Rcpp::Named("class")="stblr_raw",Rcpp::Named("version")=1),
  Rcpp::Named("meta")=Rcpp::List::create(
   Rcpp::Named("model")="bayesc",Rcpp::Named("backend")="csr_scheduled_bayesc",
   Rcpp::Named("data_level")="summary",Rcpp::Named("prior_type")="global",
   Rcpp::Named("m")=m,Rcpp::Named("nt")=nt,Rcpp::Named("n_trace")=n_trace,
   Rcpp::Named("nit")=binding.nit,Rcpp::Named("nburn")=binding.nburn,
   Rcpp::Named("nthin")=binding.nthin,Rcpp::Named("nchains")=binding.nchains,
   Rcpp::Named("keep_chains")=binding.keep_chains,Rcpp::Named("n_components")=2,
   Rcpp::Named("n_annotations")=0,Rcpp::Named("n_groups")=0),
  Rcpp::Named("marker")=marker,Rcpp::Named("trace")=trace,Rcpp::Named("variance")=variance,
  Rcpp::Named("pi")=Rcpp::List::create(Rcpp::Named("final")=pi_final,
   Rcpp::Named("mean")=pi_mean,Rcpp::Named("names")=Rcpp::CharacterVector::create("pi0","pi1")),
  Rcpp::Named("diagnostics")=diagnostics,Rcpp::Named("chains")=chains,
  Rcpp::Named("prior")=Rcpp::List::create(),Rcpp::Named("group")=Rcpp::List::create(),
  Rcpp::Named("annotation")=Rcpp::List::create(),Rcpp::Named("component")=Rcpp::List::create(),
  Rcpp::Named("selection")=selection);
 raw.attr("class")=Rcpp::CharacterVector::create("stblr_raw_v1","stblr_raw","list");
 return raw;
}

// =============================================================================
// Main exported scheduled CSR function
// =============================================================================

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_scheduled(
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
  int seed,
  int nchains,
  bool keep_chains,
  std::vector<int> chain_seeds
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nt must be positive.");
 const int m = static_cast<int>(wy[0].size());
 if (m <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: m must be positive.");
 if (nit <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nit must be positive.");
 if (nburn < 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nburn must be non-negative.");
 if (nthin <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nthin must be positive.");
 if (nchains <= 0) throw std::runtime_error("stblr_cpg_omp_csr_scheduled: nchains must be positive.");
 if (!chain_seeds.empty() && static_cast<int>(chain_seeds.size()) != nchains) {
  throw std::runtime_error("stblr_cpg_omp_csr_scheduled: chain_seeds must have length nchains.");
 }
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
 sblr::core::ScheduledExecutionControl scheduled_control;
 scheduled_control.marker_count=static_cast<std::size_t>(m);
 scheduled_control.trait_count=static_cast<std::size_t>(nt);
 scheduled_control.iterations=nit;
 scheduled_control.burnin=nburn;
 scheduled_control.thinning=nthin;
 scheduled_control.chains=nchains;
 scheduled_control.cores=ncores;
 scheduled_control.seed=seed;
 scheduled_control.chain_seeds=chain_seeds;
 scheduled_control.keep_chains=keep_chains;
 scheduled_control.sweep.full_sweep_every=full_sweep_every;
 scheduled_control.sweep.iteration_zero_is_full=true;
 scheduled_control.skip.base_interval=null_skip_base;
 scheduled_control.skip.maximum_interval=null_skip_max;
 scheduled_control.skip.burnin_only=skip_nulls_burnin_only;
 scheduled_control.candidate.probability_threshold=candidate_threshold;
 scheduled_control.candidate.lifetime=candidate_lifetime;
 scheduled_control.neighbor.enabled=wakeup_ld_neighbors;
 scheduled_control.neighbor.effect_difference_threshold=wakeup_diff_threshold;
 scheduled_control.neighbor.maximum_neighbors=wakeup_max_neighbors;
 scheduled_control.neighbor.friend_data=static_cast<const void*>(&ld);
 scheduled_control.neighbor.friend_marker_count=static_cast<std::size_t>(m);

 const sblr::core::CsrScheduledBayesCExecutionContext<STLDCSR> execution_context{
  ld, wy_mat, ww_mat, b_mat, yy_vec, ssb_prior_mat, sse_prior_mat,
  B, E, pi, n, order, d_init, r_init, scheduled_control,
  m, nt, nub, nue, adjE, pi_prior_a, pi_prior_b,
  use_d_init, use_r_init, rebuild_r_before_updateE,
  updateB, updateE, updatePi
 };
 auto execution_result=sblr::core::run_csr_scheduled_bayesc(execution_context);
 const STScheduledBayesCBindingMetadata binding_metadata{
  wy_mat,n,nit,nburn,nthin,nchains,keep_chains
 };
 return stblr_csr_scheduled_bayesc_result_to_raw(
  execution_result,binding_metadata
 );
}
