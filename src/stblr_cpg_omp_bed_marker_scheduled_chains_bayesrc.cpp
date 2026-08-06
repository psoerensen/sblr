// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "st_bed_bayesr_common.h"
#include "st_bayesrc_annotation_prior.h"
#include "blr_bed_bayesrc_core_impl.h"
#include "blr_bed_bayesrc_aggregate_impl.h"

struct BedBayesRCBindingMetadata {
 const std::vector<double>& gamma;
 int m, nt, ntrace, nit, nburn, nthin, nchains, n_used;
 double n_annotations;
 int annotation_updates_per_chain;
 bool keep_chains, return_wy, return_r;
 const std::vector<int>& convergence_markers;
};

static Rcpp::List stblr_bed_bayesrc_result_to_raw(
 const sblr::core::BedBayesRCExecutionResult& result,
 const BedBayesRCBindingMetadata& metadata
) {
 const int K=static_cast<int>(metadata.gamma.size());
 Rcpp::List cp_out(metadata.nt), marker_prior_final_out(metadata.nt);
 Rcpp::List alpha_out(metadata.nt), alpha_final_out(metadata.nt);
 for (int t=0;t<metadata.nt;++t) {
  cp_out[t]=result.comp_prob[static_cast<std::size_t>(t)];
  marker_prior_final_out[t]=result.marker_prior_final[static_cast<std::size_t>(t)];
  alpha_out[t]=result.alpha_mean[static_cast<std::size_t>(t)];
  alpha_final_out[t]=result.alpha_final[static_cast<std::size_t>(t)];
 }
 Rcpp::CharacterVector component_names(K);
 for (int k=0;k<K;++k) {
  char buf[64];
  std::snprintf(buf,sizeof(buf),"gamma_%.2f",metadata.gamma[static_cast<std::size_t>(k)]);
  component_names[k]=std::string(buf);
 }
 component_names[0]="gamma_0.00";
 Rcpp::RObject chains_out=R_NilValue;
 if (metadata.keep_chains) {
  Rcpp::List retained_chains(metadata.nt);
  for (int t=0;t<metadata.nt;++t) {
   Rcpp::List trait_chains(metadata.nchains);
   for (int ch=0;ch<metadata.nchains;++ch) {
    const auto& z=result.retained_chains[static_cast<std::size_t>(ch*metadata.nt+t)];
    Rcpp::List chain_out=Rcpp::List::create(
     Rcpp::Named("bm")=z.bm,Rcpp::Named("dm")=z.dm,
     Rcpp::Named("b")=z.b,Rcpp::Named("state")=z.state,
     Rcpp::Named("comp_prob")=z.comp_prob,
     Rcpp::Named("alpha")=z.annot_alpha_mean,
     Rcpp::Named("sigmaSqAlpha")=z.annot_sigma_mean,
     Rcpp::Named("convergence_trace")=Rcpp::List::create(
      Rcpp::Named("marker_index")=Rcpp::wrap(metadata.convergence_markers),
      Rcpp::Named("b")=z.convergence_b,Rcpp::Named("d")=z.convergence_d,
      Rcpp::Named("component")=z.convergence_component,
      Rcpp::Named("component_count")=z.convergence_aggregate.component_count,
      Rcpp::Named("realized_active_count")=z.convergence_aggregate.realized_active_count,
      Rcpp::Named("stick_eligible_count")=z.convergence_aggregate.stick_eligible_count,
      Rcpp::Named("stick_continue_count")=z.convergence_aggregate.stick_continue_count,
      Rcpp::Named("stick_stop_count")=z.convergence_aggregate.stick_stop_count,
      Rcpp::Named("alpha")=z.convergence_alpha,
      Rcpp::Named("sigmaSqAlpha")=z.convergence_sigma),
     Rcpp::Named("vbs")=z.vbs,Rcpp::Named("vgs")=z.vgs,
     Rcpp::Named("ves")=z.ves,Rcpp::Named("vle")=z.vles,
     Rcpp::Named("vld")=z.vlds,Rcpp::Named("pis")=z.pis
    );
    if (z.coupling_tempering) {
     chain_out["coupling_tempering"]=Rcpp::List::create(
      Rcpp::Named("lambda")=Rcpp::NumericVector::create(0.0,0.5,1.0),
      Rcpp::Named("replica_identity")=z.coupling_replica_identity,
      Rcpp::Named("active_count")=z.coupling_active_count,
      Rcpp::Named("expected_active_count")=z.coupling_expected_active,
      Rcpp::Named("swap")=z.coupling_swap,
      Rcpp::Named("transition_seconds")=z.coupling_transition_seconds,
      Rcpp::Named("swap_seconds")=z.coupling_swap_seconds);
    }
    trait_chains[ch]=chain_out;
   }
   retained_chains[t]=trait_chains;
  }
  chains_out=retained_chains;
 }
 Rcpp::NumericVector nsamples(result.nsamples.begin(),result.nsamples.end());
 Rcpp::NumericVector log_cpo(result.log_cpo.begin(),result.log_cpo.end());
 Rcpp::NumericVector mean_log_cpo(result.mean_log_cpo.begin(),result.mean_log_cpo.end());
 Rcpp::List raw=Rcpp::List::create(
  Rcpp::Named("schema")=Rcpp::List::create(Rcpp::Named("class")="stblr_raw",Rcpp::Named("version")=1),
  Rcpp::Named("meta")=Rcpp::List::create(
   Rcpp::Named("model")="bayesrc",Rcpp::Named("backend")="bed_bayesrc",
   Rcpp::Named("data_level")="individual",Rcpp::Named("prior_type")="annotation_component",
   Rcpp::Named("m")=metadata.m,Rcpp::Named("nt")=metadata.nt,Rcpp::Named("n_trace")=metadata.ntrace,
   Rcpp::Named("nit")=metadata.nit,Rcpp::Named("nburn")=metadata.nburn,Rcpp::Named("nthin")=metadata.nthin,
   Rcpp::Named("nchains")=metadata.nchains,Rcpp::Named("keep_chains")=metadata.keep_chains,
   Rcpp::Named("n_components")=K,Rcpp::Named("n_annotations")=metadata.n_annotations,
   Rcpp::Named("n_groups")=0,Rcpp::Named("annotations")=true,Rcpp::Named("scheduled")=false),
  Rcpp::Named("marker")=Rcpp::List::create(
   Rcpp::Named("bm")=result.bm,Rcpp::Named("dm")=result.dm,
   Rcpp::Named("wy")=metadata.return_wy ? Rcpp::wrap(result.wy) : R_NilValue,
   Rcpp::Named("r")=metadata.return_r ? Rcpp::wrap(result.residual_marker_score) : R_NilValue,
   Rcpp::Named("b")=result.b,Rcpp::Named("state")=result.state),
  Rcpp::Named("trace")=Rcpp::List::create(
   Rcpp::Named("vbs")=result.vbs,Rcpp::Named("vgs")=result.vgs,Rcpp::Named("ves")=result.ves,
   Rcpp::Named("vle")=result.vle,Rcpp::Named("vld")=result.vld,Rcpp::Named("pis")=result.pis),
  Rcpp::Named("variance")=Rcpp::List::create(
   Rcpp::Named("covb")=arma::diagmat(result.final_vb),Rcpp::Named("covg")=arma::diagmat(result.final_vg),
   Rcpp::Named("cove")=arma::diagmat(result.final_ve),Rcpp::Named("vb")=arma::diagmat(result.final_vb),
   Rcpp::Named("vg")=arma::diagmat(result.final_vg),Rcpp::Named("ve")=arma::diagmat(result.final_ve)),
  Rcpp::Named("pi")=Rcpp::List::create(Rcpp::Named("final")=result.final_prior,
   Rcpp::Named("mean")=result.mean_prior,Rcpp::Named("names")=component_names),
  Rcpp::Named("diagnostics")=Rcpp::List::create(
   Rcpp::Named("nsamples")=nsamples,Rcpp::Named("n_used")=metadata.n_used,
   Rcpp::Named("log_cpo")=log_cpo,Rcpp::Named("mean_log_cpo")=mean_log_cpo,
   Rcpp::Named("full_sweeps")=true,Rcpp::Named("adaptive_skipping")=false,
   Rcpp::Named("annotation_updates_per_chain")=metadata.annotation_updates_per_chain,
   Rcpp::Named("ld_swap")=R_NilValue),
  Rcpp::Named("chains")=chains_out,Rcpp::Named("prior")=Rcpp::List::create(),
  Rcpp::Named("group")=Rcpp::List::create(),
  Rcpp::Named("annotation")=Rcpp::List::create(
   Rcpp::Named("annotation_names")=R_NilValue,Rcpp::Named("alpha_mean")=alpha_out,
   Rcpp::Named("alpha_final")=alpha_final_out,
   Rcpp::Named("sigmaSqAlpha_mean")=result.sigma_mean.t(),
   Rcpp::Named("sigmaSqAlpha_final")=result.sigma_final.t(),
   Rcpp::Named("marker_prior_final")=marker_prior_final_out),
  Rcpp::Named("component")=Rcpp::List::create(
   Rcpp::Named("names")=component_names,Rcpp::Named("mixture_var")=metadata.gamma,
   Rcpp::Named("prob")=cp_out,Rcpp::Named("ncomp")=result.component_counts,
   Rcpp::Named("dm_component_mean")=result.component_mean),
  Rcpp::Named("selection")=Rcpp::List::create(Rcpp::Named("enabled")=false,
   Rcpp::Named("fixed")=false,Rcpp::Named("trace")=R_NilValue)
 );
 raw.attr("class")=Rcpp::CharacterVector::create("stblr_raw_v1","stblr_raw","list");
 return raw;
}

// Development-only hooks for validating the coupling target and exchange
// algebra. They are not exported from the package namespace.
// [[Rcpp::export(name = ".st_bayesrc_tempered_probabilities")]]
arma::mat st_bayesrc_tempered_probabilities_test(
 arma::mat annotation, arma::mat alpha, arma::vec baseline_intercept,
 double coupling, double probability_floor
) {
 return st_bayesrc_compute_tempered_snp_pi(
  annotation,alpha,baseline_intercept,coupling,probability_floor);
}

// [[Rcpp::export(name = ".st_bayesrc_tempered_annotation_update")]]
Rcpp::List st_bayesrc_tempered_annotation_update_test(
 arma::mat annotation, arma::rowvec component_numeric, arma::mat alpha,
 arma::vec sigma_sq_alpha, arma::mat intercept_prior_resolved,
 double sigma_alpha_a, double sigma_alpha_b, arma::vec baseline_intercept,
 double coupling, int seed
) {
 arma::Row<int> component(component_numeric.n_elem);
 for (arma::uword index=0;index<component_numeric.n_elem;++index) {
  const double value=component_numeric(index);
  if (!std::isfinite(value) || value<0.0 || value!=std::floor(value))
   throw std::invalid_argument("component must contain non-negative integers.");
  component(index)=static_cast<int>(value);
 }
 const auto prior=st_bayesrc_parse_intercept_prior(
  intercept_prior_resolved,static_cast<int>(alpha.n_cols));
 std::mt19937 generator(static_cast<std::mt19937::result_type>(seed));
 const auto diagnostics=st_bayesrc_update_tempered_annotation_prior(
  annotation,component,alpha,sigma_sq_alpha,prior,sigma_alpha_a,
  sigma_alpha_b,baseline_intercept,coupling,generator);
 return Rcpp::List::create(
  Rcpp::Named("alpha")=alpha,Rcpp::Named("sigmaSqAlpha")=sigma_sq_alpha,
  Rcpp::Named("eligible")=diagnostics.eligible,
  Rcpp::Named("continuation")=diagnostics.continuation,
  Rcpp::Named("prior_only")=diagnostics.prior_only);
}

// [[Rcpp::export(name = ".st_bayesrc_coupling_swap_log_ratio")]]
double st_bayesrc_coupling_swap_log_ratio_test(
 arma::mat annotation, arma::mat alpha_lower, arma::rowvec component_lower_numeric,
 arma::mat alpha_upper, arma::rowvec component_upper_numeric,
 arma::vec baseline_intercept, double coupling_lower, double coupling_upper,
 double probability_floor
) {
 if (component_lower_numeric.n_elem!=component_upper_numeric.n_elem)
  throw std::invalid_argument("coupling swap component states must be aligned.");
 arma::Row<int> component_lower(component_lower_numeric.n_elem);
 arma::Row<int> component_upper(component_upper_numeric.n_elem);
 for (arma::uword index=0;index<component_lower.n_elem;++index) {
  const double lower=component_lower_numeric(index);
  const double upper=component_upper_numeric(index);
  if (!std::isfinite(lower) || !std::isfinite(upper) || lower<0.0 || upper<0.0 ||
      lower!=std::floor(lower) || upper!=std::floor(upper))
   throw std::invalid_argument("coupling swap components must be non-negative integers.");
  component_lower(index)=static_cast<int>(lower);
  component_upper(index)=static_cast<int>(upper);
 }
 return
  st_bayesrc_tempered_log_allocation_prior(
   annotation,alpha_upper,baseline_intercept,coupling_lower,
   probability_floor,component_upper)+
  st_bayesrc_tempered_log_allocation_prior(
   annotation,alpha_lower,baseline_intercept,coupling_upper,
   probability_floor,component_lower)-
  st_bayesrc_tempered_log_allocation_prior(
   annotation,alpha_lower,baseline_intercept,coupling_lower,
   probability_floor,component_lower)-
  st_bayesrc_tempered_log_allocation_prior(
   annotation,alpha_upper,baseline_intercept,coupling_upper,
   probability_floor,component_upper);
}

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc(
  Rcpp::CharacterVector bed_files, int n, Rcpp::List cls,
  Rcpp::NumericMatrix y, std::vector<std::vector<double>> b_init,
  std::vector<int> sets, Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::Nullable<Rcpp::List> af, bool scale, arma::mat B, arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior, arma::mat A,
  std::vector<double> gamma, arma::mat annot_alpha_init,
  arma::vec annot_sigma_sq_alpha_init, arma::mat intercept_prior_resolved,
  double sigmaSqAlpha_a = 2.0, double sigmaSqAlpha_b = 2.0,
  double pi_floor = 1e-12, double nub = 4.0, double nue = 4.0,
  bool updateAlpha = true, bool updateB = true, bool updateE = true,
  int annot_alpha_update_every = 10, double adjE = 0.9, int nit = 1000,
  int nburn = 100, int nthin = 1, int rebuild_every = 100,
  bool return_wy = true, bool return_r = true, int read_block_size = 256,
  int nchains = 1, bool keep_chains = false, int ncores = 1, int seed = 10,
  Rcpp::IntegerVector chain_seeds=Rcpp::IntegerVector::create(),
  Rcpp::IntegerVector convergence_markers=Rcpp::IntegerVector::create(),
  bool convergence_annotations=false, bool convergence_b=false,
  bool convergence_d=false, bool convergence_component=false
) {
 if (nit <= 0 || nburn < 0 || nthin <= 0 || nchains <= 0 || ncores <= 0)
  throw std::runtime_error("invalid MCMC or chain controls.");
 if (chain_seeds.size()>0 && chain_seeds.size()!=nchains)
  throw std::runtime_error("chain_seeds must be empty or have length nchains.");
 if (annot_alpha_update_every <= 0)
  throw std::runtime_error("annot_alpha_update_every must be positive.");
 if (!std::isfinite(pi_floor) || pi_floor <= 0.0 || pi_floor >= 1.0)
  throw std::runtime_error("pi_floor must be in (0, 1).");
 const int K = static_cast<int>(gamma.size());
 const std::vector<int> convergence_markers_cpp=
  Rcpp::as<std::vector<int>>(convergence_markers);
 if (K < 2 || gamma[0] != 0.0) throw std::runtime_error("gamma must start with exact zero.");
 for (int k = 1; k < K; ++k)
  if (!std::isfinite(gamma[k]) || gamma[k] <= 0.0)
   throw std::runtime_error("non-null gamma values must be positive finite.");
 std::vector<std::string> beds = br_copy_bed_files(bed_files);
 std::vector<std::vector<int>> cls_cpp = br_copy_int_list(cls);
 std::vector<int> rows0 = br_copy_rows0(rows, n);
 FastPackedBedMatrixBR G = br_read_bed_blocked(
  beds, n, rows0.empty() ? nullptr : rows0.data(), static_cast<int>(rows0.size()),
  cls_cpp, read_block_size, std::max(1, ncores)
 );
 const int m = G.m, nt = y.ncol(), n_used = G.n;
 for (int marker : convergence_markers_cpp) {
  if (marker < 0 || marker >= m) {
   throw std::runtime_error(
    "stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc: convergence marker index is out of range.");
  }
 }
 if (y.nrow() != n_used || nt <= 0) throw std::runtime_error("invalid phenotype dimensions.");
 if (A.n_rows != static_cast<arma::uword>(m) || A.n_cols == 0 || !A.is_finite())
  throw std::runtime_error("A must be a finite m by n_annotations matrix aligned to selected BED markers.");
 if (annot_alpha_init.n_rows != A.n_cols || annot_alpha_init.n_cols != static_cast<arma::uword>(K - 1))
  throw std::runtime_error("annot_alpha_init dimensions must be n_annotations by K-1.");
 if (annot_sigma_sq_alpha_init.n_elem != static_cast<arma::uword>(K - 1))
  throw std::runtime_error("annot_sigma_sq_alpha_init must have length K-1.");
 if (!annot_alpha_init.is_finite() || !annot_sigma_sq_alpha_init.is_finite() ||
     arma::any(annot_sigma_sq_alpha_init <= 0.0))
  throw std::runtime_error("annotation initial values must be finite and variances positive.");
 if (static_cast<int>(sets.size()) != m || static_cast<int>(b_init.size()) != nt)
  throw std::runtime_error("sets or b_init dimensions do not match selected markers.");
 for (int t = 0; t < nt; ++t)
  if (static_cast<int>(b_init[t].size()) != m)
   throw std::runtime_error("each b_init trait must have length m.");
 if (B.n_rows != static_cast<arma::uword>(nt) || B.n_cols != static_cast<arma::uword>(nt) ||
     E.n_rows != static_cast<arma::uword>(nt) || E.n_cols != static_cast<arma::uword>(nt))
  throw std::runtime_error("B and E must be nt by nt.");
 arma::mat y_mat(y.begin(), n_used, nt, true);
 arma::mat ssb(nt, nt, arma::fill::zeros), sse(nt, nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) for (int u = 0; u < nt; ++u) {
  ssb(t, u) = ssb_prior[t][u]; sse(t, u) = sse_prior[t][u];
 }
 std::vector<double> af_cpp = br_flatten_af_list(af);
 if (af_cpp.empty()) af_cpp = br_compute_af(G);
 std::vector<MarkerMapBayesR> maps = br_build_marker_maps(G, af_cpp, scale, std::max(1, ncores));
 std::vector<int> order = br_make_marker_order(sets, m);
 const StBayesRCInterceptPrior intercept_prior =
  st_bayesrc_parse_intercept_prior(intercept_prior_resolved, K - 1);
 const int njobs = nt * nchains;
 std::vector<sblr::core::BedBayesRCChainExecutionResult> jobs(static_cast<std::size_t>(njobs));
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores) schedule(static)
#endif
 for (int job = 0; job < njobs; ++job) {
  const auto task=sblr::core::make_bed_family_task_index(job,nt);
  const int trait=task.trait, chain=task.chain;
  const std::uint64_t chain_seed=chain_seeds.size()==0 ?
   sblr::core::resolve_bed_family_logical_chain_seed(seed,trait,chain) :
   static_cast<unsigned int>(chain_seeds[static_cast<std::size_t>(chain)]+
                             1000003*(trait+1));
  const sblr::core::BedBayesRCChainExecutionContext<
   FastPackedBedMatrixBR,arma::mat,MarkerMapBayesR
  > context{
   {G,static_cast<std::size_t>(m),static_cast<std::size_t>(n_used),G.nbytes},
   {A,static_cast<std::size_t>(m),static_cast<std::size_t>(A.n_cols),0u},
   {gamma,0u,static_cast<std::size_t>(K-1)},
   {annot_alpha_init,annot_sigma_sq_alpha_init,intercept_prior,
    sigmaSqAlpha_a,sigmaSqAlpha_b,updateAlpha,annot_alpha_update_every},
   maps,order,y_mat,b_init,B,E,ssb,sse,pi_floor,nub,nue,adjE,
   updateB,updateE,nit,nburn,nthin,rebuild_every,chain_seed,trait,chain,
   convergence_markers_cpp,convergence_annotations,convergence_b,
   convergence_d,convergence_component
  };
  jobs[job]=sblr::core::run_bed_bayesrc_chain(context);
 }
 for (int job = 0; job < njobs; ++job)
  if (jobs[job].failed) throw std::runtime_error(jobs[job].error);
 const int ntrace = nit + nburn;
 const sblr::core::BedBayesRCAggregationContext aggregation_context{
  A,gamma,static_cast<std::size_t>(m),static_cast<std::size_t>(A.n_cols),
  static_cast<std::size_t>(K),static_cast<std::size_t>(nt),
  static_cast<std::size_t>(nchains),static_cast<std::size_t>(ntrace),pi_floor,keep_chains
 };
 sblr::core::BedBayesRCExecutionResult result=
  sblr::core::aggregate_bed_bayesrc_results(jobs,aggregation_context);
 if (!result.failures.empty()) throw std::runtime_error(result.failures.front());
 result.wy.zeros(m,nt); result.residual_marker_score.zeros(m,nt);
 if (return_wy || return_r) {
  for (int t = 0; t < nt; ++t) {
   arma::vec y_t = y_mat.col(t);
   arma::vec residual = y_t - br_xb(G,maps,order,result.b.col(t).t());
   for (int j = 0; j < m; ++j) {
    if (return_wy) result.wy(j,t)=br_dot_residual(G,j,maps[j],y_t.memptr());
    if (return_r) result.residual_marker_score(j,t)=
     br_dot_residual(G,j,maps[j],residual.memptr());
   }
  }
 }
 const BedBayesRCBindingMetadata metadata{
  gamma,m,nt,ntrace,nit,nburn,nthin,nchains,n_used,static_cast<double>(A.n_cols),
  updateAlpha ? ((nit+nburn)/annot_alpha_update_every) *
   intercept_prior.annotation_updates_per_cycle : 0,
  keep_chains,return_wy,return_r,convergence_markers_cpp
 };
 return stblr_bed_bayesrc_result_to_raw(result,metadata);
}
