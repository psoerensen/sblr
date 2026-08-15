// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "packed_bed.h"
#include "st_bed_bayesr_common.h"
#include "blr_bed_bayesr_types.h"

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

#include "blr_bed_bayesr_core_impl.h"
#include "blr_bed_bayesr_aggregate_impl.h"

struct BedBayesRBindingMetadata {
 int marker_count, trait_count, trace_length, iterations, burnin, thinning;
 int chain_count, component_count, sample_count;
 const std::vector<double>& component_scales;
 const std::vector<int>& convergence_markers;
};

static Rcpp::List stblr_bed_bayesr_result_to_raw(
 const sblr::core::BedBayesRExecutionResult& result,
 const BedBayesRBindingMetadata& metadata,
 const std::vector<sblr::core::BedBayesRChainExecutionResult>& chain_results
);
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
  int seed,
  Rcpp::IntegerVector chain_seeds=Rcpp::IntegerVector::create(),
  Rcpp::IntegerVector convergence_markers=Rcpp::IntegerVector::create(),
  bool convergence_probability=false,
  bool convergence_b=false,
  bool convergence_d=false,
  bool convergence_component=false,
  Rcpp::Nullable<Rcpp::List> execution_contract=R_NilValue
) {
 if (nit <= 0 || nburn < 0)
  throw std::runtime_error("stblr_cpg_omp_bed_marker_scheduled_chains_bayesr: nit must be positive and nburn non-negative.");
 if (nthin <= 0) throw std::runtime_error("nthin must be positive.");
 if (rebuild_every < 0) throw std::runtime_error("rebuild_every must be >= 0.");
 if (full_sweep_every < 0) throw std::runtime_error("full_sweep_every must be >= 0.");
 if (null_skip_base <= 0) throw std::runtime_error("null_skip_base must be positive.");
 if (null_skip_max < 0) throw std::runtime_error("null_skip_max must be >= 0.");
 if (nchains <= 0) throw std::runtime_error("nchains must be positive.");
 if (chain_seeds.size()>0 && chain_seeds.size()!=nchains)
  throw std::runtime_error("chain_seeds must be empty or have length nchains.");
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
 const std::vector<int> convergence_markers_cpp=Rcpp::as<std::vector<int>>(convergence_markers);

 const int* rows0_ptr = rows0.empty() ? nullptr : rows0.data();
 const int n_rows = rows0.empty() ? 0 : static_cast<int>(rows0.size());

 int read_nthreads = 1;
#ifdef _OPENMP
 read_nthreads = ncores > 0 ? std::max(1, ncores) : std::max(1, omp_get_max_threads());
#endif


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

 for (int marker : convergence_markers_cpp) {
  if (marker < 0 || marker >= m) {
   throw std::runtime_error(
    "stblr_cpg_omp_bed_marker_scheduled_chains_bayesr: convergence marker index is out of range.");
  }
 }

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
 const BlrPhase3ExecutionContract phase3 =
  parse_blr_phase3_execution_contract(execution_contract, njobs, nit);
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


 std::vector<MarkerMapBayesR> marker_maps = br_build_marker_maps(G, af_cpp, scale, nthreads);
 std::vector<int> marker_order = br_make_marker_order(sets, m);


#ifdef _OPENMP
#else
#endif

 arma::mat y_mat(n_used, nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t)
  for (int i = 0; i < n_used; ++i) y_mat(i, t) = y(i, t);

 std::vector<sblr::core::BedBayesRChainExecutionResult> job_results(static_cast<std::size_t>(njobs));
 std::vector<int> worker_ids(static_cast<std::size_t>(njobs),0);
 std::vector<int> team_sizes(static_cast<std::size_t>(njobs),1);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int job = 0; job < njobs; ++job) {
  const auto task=sblr::core::make_bed_family_task_index(job,nt);
  const int ch=task.chain;
  const int t=task.trait;
  const int canonical_task=t*nchains+ch;
#ifdef _OPENMP
  worker_ids[static_cast<std::size_t>(canonical_task)]=omp_get_thread_num();
  team_sizes[static_cast<std::size_t>(canonical_task)]=omp_get_num_threads();
#endif

  const sblr::core::BedBayesRPackedGenotypeView<FastPackedBedMatrixBR> genotype{
   G, G.data.data(), G.data.size(), static_cast<std::size_t>(G.m),
   static_cast<std::size_t>(G.n), G.nbytes, G.stride
  };
  const sblr::core::BedBayesRComponentSpec components{0u, c, pi, alpha};
  const sblr::core::BedBayesRSchedulerControl scheduler{
   sblr::core::ScheduledSweepControl{full_sweep_every},
   sblr::core::NullSkipControl{null_skip_base, null_skip_max,
                              skip_nulls_burnin_only, "probability_adaptive"},
   sblr::core::CandidateControl{candidate_threshold, candidate_lifetime}
  };
  const std::uint64_t legacy_seed=chain_seeds.size()==0 ?
   sblr::core::resolve_bed_family_logical_chain_seed(seed,t,ch) :
   static_cast<unsigned int>(chain_seeds[static_cast<std::size_t>(ch)]+
                             1000003*(t+1));
  const std::uint64_t chain_seed=blr_phase3_task_seed(
   phase3,canonical_task,static_cast<std::uint32_t>(legacy_seed));
  const sblr::core::BedBayesRChainExecutionContext<FastPackedBedMatrixBR,MarkerMapBayesR> context{
   genotype, marker_maps, marker_order, y_mat, b_init, B, E,
   ssb_prior_mat, sse_prior_mat, components, scheduler, nub, nue, adjE,
   nit, nburn, nthin, rebuild_every, progress_every, chain_seed, t, ch,
   updateB, updateE, updatePi, convergence_markers_cpp,
   convergence_probability, convergence_b, convergence_d,
   convergence_component,&phase3
  };
  job_results[static_cast<std::size_t>(job)] = sblr::core::run_bed_bayesr_chain(context);
 }

 int failed_total = 0;
 for (int job = 0; job < njobs; ++job) {
  const sblr::core::BedBayesRChainExecutionResult& r = job_results[static_cast<std::size_t>(job)];
  const auto task=sblr::core::make_bed_family_task_index(job,nt);
  const int ch=task.chain;
  const int t=task.trait;
  (void)ch;
  (void)t;
  for (const sblr::core::BedBayesRProgressEvent& event : r.progress_events) {
   (void)event;
  }
  if (r.failed) {
   ++failed_total;
  }
 }

 if (failed_total > 0) {
  for (int job = 0; job < njobs; ++job) {
   const sblr::core::BedBayesRChainExecutionResult& r = job_results[static_cast<std::size_t>(job)];
   if (r.failed) {
    const auto task=sblr::core::make_bed_family_task_index(job,nt);
    const int ch=task.chain;
    const int t=task.trait;
    throw std::runtime_error(
      "stblr_cpg_omp_bed_marker_scheduled_chains_bayesr failed for chain " +
       std::to_string(ch) + ", trait " + std::to_string(t) + ": " + r.error
    );
   }
  }
 }

 const sblr::core::BedBayesRAggregationContext aggregation_context{
  static_cast<std::size_t>(m), static_cast<std::size_t>(nt),
  static_cast<std::size_t>(nchains), static_cast<std::size_t>(nit+nburn),
  static_cast<std::size_t>(K), 0u
 };
 sblr::core::BedBayesRExecutionResult result=
  sblr::core::aggregate_bed_bayesr_results(job_results,aggregation_context);

 if (return_wy || return_r) {
  for (int t = 0; t < nt; ++t) {
   arma::vec y_t = y_mat.col(static_cast<arma::uword>(t));
   arma::rowvec b_t = result.final_effects.row(static_cast<arma::uword>(t));
   arma::vec xb_t = br_xb(G, marker_maps, marker_order, b_t);
   arma::vec e_t = y_t - xb_t;

   for (int j = 0; j < m; ++j) {
    if (return_wy)
     result.wy(static_cast<arma::uword>(t), static_cast<arma::uword>(j)) =
      br_dot_residual(G, j, marker_maps[static_cast<std::size_t>(j)], y_t.memptr());
    if (return_r)
     result.residual_score(static_cast<arma::uword>(t), static_cast<arma::uword>(j)) =
      br_dot_residual(G, j, marker_maps[static_cast<std::size_t>(j)], e_t.memptr());
   }
  }
 }

 const BedBayesRBindingMetadata binding_metadata{
  m,nt,nit+nburn,nit,nburn,nthin,nchains,K,n_used,c,
  convergence_markers_cpp
 };
 Rcpp::List raw=stblr_bed_bayesr_result_to_raw(
  result,binding_metadata,job_results);
 if (phase3.active()) {
  Rcpp::List diagnostics=raw["diagnostics"];
  diagnostics.push_back(blr_phase3_worker_diagnostics(
   phase3,ncores,nthreads,worker_ids,team_sizes),"workers");
  raw["diagnostics"]=diagnostics;
 }
 return raw;
}

static Rcpp::List stblr_bed_bayesr_result_to_raw(
 const sblr::core::BedBayesRExecutionResult& result,
 const BedBayesRBindingMetadata& metadata,
 const std::vector<sblr::core::BedBayesRChainExecutionResult>& chain_results
) {
 const int m=metadata.marker_count, nt=metadata.trait_count;
 const int n_trace=metadata.trace_length, nit=metadata.iterations;
 const int nburn=metadata.burnin, nthin=metadata.thinning;
 const int nchains=metadata.chain_count, K=metadata.component_count;
 const int n_used=metadata.sample_count;
 const std::vector<double>& c=metadata.component_scales;
 const std::vector<int>& convergence_markers=metadata.convergence_markers;


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
    result.final_pi(static_cast<arma::uword>(k), static_cast<arma::uword>(t));
   mean_pi(static_cast<arma::uword>(t), static_cast<arma::uword>(k)) =
    result.mean_pi(static_cast<arma::uword>(k), static_cast<arma::uword>(t));
  }
 }

 Rcpp::List comp_prob_out(nt);
 for (int t = 0; t < nt; ++t) {
  comp_prob_out[t] = arma::mat(result.component_probability[static_cast<std::size_t>(t)].t());
 }

 Rcpp::NumericVector nsamples_out(nt);
 Rcpp::IntegerVector n_used_out(nt);
 Rcpp::NumericVector log_cpo_out(nt);
 Rcpp::NumericVector mean_log_cpo_out(nt);
 Rcpp::NumericVector seconds_mean_out(nt);
 Rcpp::NumericVector seconds_max_out(nt);

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  nsamples_out[t] = result.retained_samples(tu);
  n_used_out[t] = n_used;
  log_cpo_out[t] = result.log_cpo(tu);
  mean_log_cpo_out[t] = result.mean_log_cpo(tu);
  seconds_mean_out[t] = result.seconds_mean(tu);
  seconds_max_out[t] = result.seconds_max(tu);
 }

 Rcpp::List marker = Rcpp::List::create(
  Rcpp::Named("bm") = marker_matrix(result.bm),
  Rcpp::Named("dm") = marker_matrix(result.dm),
  Rcpp::Named("wy") = marker_matrix(result.wy),
  Rcpp::Named("r") = marker_matrix(result.residual_score),
  Rcpp::Named("b") = marker_matrix(result.final_effects),
  Rcpp::Named("state") = marker_matrix(result.final_states),
  Rcpp::Named("bm_sd") = marker_matrix(result.bm_sd),
  Rcpp::Named("bm_min") = marker_matrix(result.bm_min),
  Rcpp::Named("bm_max") = marker_matrix(result.bm_max),
  Rcpp::Named("dm_sd") = marker_matrix(result.dm_sd),
  Rcpp::Named("dm_min") = marker_matrix(result.dm_min),
  Rcpp::Named("dm_max") = marker_matrix(result.dm_max)
 );

 Rcpp::List trace = Rcpp::List::create(
  Rcpp::Named("vbs") = trace_matrix(result.vbs),
  Rcpp::Named("vgs") = trace_matrix(result.vgs),
  Rcpp::Named("ves") = trace_matrix(result.ves),
  Rcpp::Named("vle") = trace_matrix(result.vles),
  Rcpp::Named("vld") = trace_matrix(result.vlds),
  Rcpp::Named("pis") = R_NilValue
 );

 Rcpp::List variance = Rcpp::List::create(
  Rcpp::Named("covb") = diagonal_matrix(result.final_vb),
  Rcpp::Named("covg") = diagonal_matrix(result.final_vg),
  Rcpp::Named("cove") = diagonal_matrix(result.final_ve),
  Rcpp::Named("vb") = diagonal_matrix(result.final_vb),
  Rcpp::Named("vg") = diagonal_matrix(result.final_vg),
  Rcpp::Named("ve") = diagonal_matrix(result.final_ve)
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
 Rcpp::List chains(nt);
 for (int t=0; t<nt; ++t) {
  Rcpp::List trait_chains(nchains);
  for (int chain=0; chain<nchains; ++chain) {
   const auto& value=chain_results[
    static_cast<std::size_t>(chain*nt+t)];
   trait_chains[chain]=Rcpp::List::create(
    Rcpp::Named("marker")=Rcpp::List::create(
     Rcpp::Named("bm")=value.bm.t(),Rcpp::Named("dm")=value.dm.t(),
     Rcpp::Named("state")=value.d_as_double.t()),
    Rcpp::Named("trace")=Rcpp::List::create(
     Rcpp::Named("vbs")=value.vbs.t(),Rcpp::Named("vgs")=value.vgs.t(),
     Rcpp::Named("ves")=value.ves.t(),Rcpp::Named("vle")=value.vles.t(),
     Rcpp::Named("vld")=value.vlds.t()),
    Rcpp::Named("component")=Rcpp::List::create(
     Rcpp::Named("prob")=value.pip_k.t(),
     Rcpp::Named("dm_component_mean")=value.component_mean.t()),
    Rcpp::Named("pi")=Rcpp::List::create(
     Rcpp::Named("final")=value.final_pi,
     Rcpp::Named("mean")=value.mean_pi),
    Rcpp::Named("convergence_trace")=Rcpp::List::create(
     Rcpp::Named("component_pi")=value.convergence_pi,
     Rcpp::Named("b")=value.convergence_b,
     Rcpp::Named("d")=value.convergence_d,
     Rcpp::Named("component")=value.convergence_component,
     Rcpp::Named("component_count")=value.convergence_aggregate.component_count,
     Rcpp::Named("realized_active_count")=value.convergence_aggregate.realized_active_count,
     Rcpp::Named("stick_eligible_count")=value.convergence_aggregate.stick_eligible_count,
     Rcpp::Named("stick_continue_count")=value.convergence_aggregate.stick_continue_count,
     Rcpp::Named("stick_stop_count")=value.convergence_aggregate.stick_stop_count,
     Rcpp::Named("marker_index")=convergence_markers),
    Rcpp::Named("diagnostics")=Rcpp::List::create(
     Rcpp::Named("seconds")=value.seconds));
  }
  chains[t]=trait_chains;
 }

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
   Rcpp::Named("keep_chains") = true,
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
  Rcpp::Named("chains") = chains,
  Rcpp::Named("prior") = Rcpp::List::create(),
  Rcpp::Named("group") = Rcpp::List::create(),
  Rcpp::Named("annotation") = Rcpp::List::create(),
  Rcpp::Named("component") = Rcpp::List::create(
   Rcpp::Named("names") = component_names,
   Rcpp::Named("mixture_var") = mixture_var_vec,
   Rcpp::Named("prob") = comp_prob_out,
   Rcpp::Named("ncomp") = R_NilValue,
   Rcpp::Named("dm_component_mean") = marker_matrix(result.component_mean)
  ),
  Rcpp::Named("selection") = selection
 );
 raw.attr("class") = Rcpp::CharacterVector::create("stblr_raw_v1", "stblr_raw", "list");
 return raw;
}
