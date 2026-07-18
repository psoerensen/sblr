// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "st_bed_bayesr_common.h"
#include "st_bayesrc_annotation_prior.h"
#include "blr_bed_bayesrc_core_impl.h"

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_bed_marker_scheduled_chains_bayesrc(
  Rcpp::CharacterVector bed_files, int n, Rcpp::List cls,
  Rcpp::NumericMatrix y, std::vector<std::vector<double>> b_init,
  std::vector<int> sets, Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::Nullable<Rcpp::List> af, bool scale, arma::mat B, arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior, arma::mat A,
  std::vector<double> gamma, arma::mat annot_alpha_init,
  arma::vec annot_sigma_sq_alpha_init, bool intercept_flat = true,
  double sigmaSqAlpha_a = 2.0, double sigmaSqAlpha_b = 2.0,
  double pi_floor = 1e-12, double nub = 4.0, double nue = 4.0,
  bool updateAlpha = true, bool updateB = true, bool updateE = true,
  int annot_alpha_update_every = 10, double adjE = 0.9, int nit = 1000,
  int nburn = 100, int nthin = 1, int rebuild_every = 100,
  bool return_wy = true, bool return_r = true, int read_block_size = 256,
  int nchains = 1, bool keep_chains = false, int ncores = 1, int seed = 10
) {
 if (nit <= 0 || nburn < 0 || nthin <= 0 || nchains <= 0 || ncores <= 0)
  throw std::runtime_error("invalid MCMC or chain controls.");
 if (annot_alpha_update_every <= 0)
  throw std::runtime_error("annot_alpha_update_every must be positive.");
 if (!std::isfinite(pi_floor) || pi_floor <= 0.0 || pi_floor >= 1.0)
  throw std::runtime_error("pi_floor must be in (0, 1).");
 const int K = static_cast<int>(gamma.size());
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
 const int njobs = nt * nchains;
 std::vector<sblr::core::BedBayesRCChainExecutionResult> jobs(static_cast<std::size_t>(njobs));
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores) schedule(static)
#endif
 for (int job = 0; job < njobs; ++job) {
  const int trait=job%nt, chain=job/nt;
  const std::uint64_t chain_seed=static_cast<unsigned int>(
   seed+1000003*(trait+1)+9176*(chain+1)
  );
  const sblr::core::BedBayesRCChainExecutionContext<
   FastPackedBedMatrixBR,arma::mat,MarkerMapBayesR
  > context{
   {G,static_cast<std::size_t>(m),static_cast<std::size_t>(n_used),G.nbytes},
   {A,static_cast<std::size_t>(m),static_cast<std::size_t>(A.n_cols),0u},
   {gamma,0u,static_cast<std::size_t>(K-1)},
   {annot_alpha_init,annot_sigma_sq_alpha_init,intercept_flat,
    sigmaSqAlpha_a,sigmaSqAlpha_b,updateAlpha,annot_alpha_update_every},
   maps,order,y_mat,b_init,B,E,ssb,sse,pi_floor,nub,nue,adjE,
   updateB,updateE,nit,nburn,nthin,rebuild_every,chain_seed,trait,chain
  };
  jobs[job]=sblr::core::run_bed_bayesrc_chain(context);
 }
 for (int job = 0; job < njobs; ++job)
  if (jobs[job].failed) throw std::runtime_error(jobs[job].error);
 const int ntrace = nit + nburn;
 arma::mat bm(m, nt, arma::fill::zeros), dm(m, nt, arma::fill::zeros);
 arma::mat b(m, nt, arma::fill::zeros), state(m, nt, arma::fill::zeros);
 arma::mat component_mean(m, nt, arma::fill::zeros);
 arma::mat vbs(ntrace, nt, arma::fill::zeros), vgs(ntrace, nt, arma::fill::zeros);
 arma::mat ves(ntrace, nt, arma::fill::zeros), vle(ntrace, nt, arma::fill::zeros);
 arma::mat vld(ntrace, nt, arma::fill::zeros), pis(ntrace, nt, arma::fill::zeros);
 arma::mat final_prior(nt, K, arma::fill::zeros), mean_prior(nt, K, arma::fill::zeros);
 std::vector<arma::mat> comp_prob(static_cast<std::size_t>(nt));
 std::vector<arma::mat> marker_prior_final(static_cast<std::size_t>(nt));
 std::vector<arma::mat> alpha_mean(static_cast<std::size_t>(nt));
 std::vector<arma::mat> alpha_final(static_cast<std::size_t>(nt));
 arma::mat sigma_mean(nt, K - 1, arma::fill::zeros);
 arma::mat sigma_final(nt, K - 1, arma::fill::zeros);
 arma::vec final_vb(nt, arma::fill::zeros), final_vg(nt, arma::fill::zeros), final_ve(nt, arma::fill::zeros);
 Rcpp::NumericVector log_cpo(nt), mean_log_cpo(nt), nsamples(nt);
 for (int t = 0; t < nt; ++t) {
  comp_prob[t].zeros(m, K); marker_prior_final[t].zeros(m, K);
  alpha_mean[t].zeros(A.n_cols, K - 1);
  alpha_final[t].zeros(A.n_cols, K - 1);
  for (int ch = 0; ch < nchains; ++ch) {
   const sblr::core::BedBayesRCChainExecutionResult& z = jobs[ch * nt + t];
   bm.col(t) += z.bm.t(); dm.col(t) += z.dm.t(); b.col(t) += z.b.t();
   state.col(t) += z.state.t(); component_mean.col(t) += z.component_mean.t();
   vbs.col(t) += z.vbs.t(); vgs.col(t) += z.vgs.t(); ves.col(t) += z.ves.t();
   vle.col(t) += z.vles.t(); vld.col(t) += z.vlds.t(); pis.col(t) += z.pis.t();
   comp_prob[t] += z.comp_prob; alpha_mean[t] += z.annot_alpha_mean;
   alpha_final[t] += z.annot_alpha_final;
   sigma_mean.row(t) += z.annot_sigma_mean.t(); mean_prior.row(t) += z.mean_prior;
   sigma_final.row(t) += z.annot_sigma_final.t();
   marker_prior_final[t] += st_bayesrc_compute_snp_pi(A, z.annot_alpha_final, pi_floor);
   final_prior.row(t) += arma::mean(st_bayesrc_compute_snp_pi(A, z.annot_alpha_final, pi_floor), 0);
   final_vb(t) += z.final_vb; final_vg(t) += z.final_vg; final_ve(t) += z.final_ve;
   log_cpo[t] += z.log_cpo; mean_log_cpo[t] += z.mean_log_cpo; nsamples[t] += z.nsamples;
  }
 }
 const double inv = 1.0 / nchains;
 bm *= inv; dm *= inv; b *= inv; state *= inv; component_mean *= inv;
 vbs *= inv; vgs *= inv; ves *= inv; vle *= inv; vld *= inv; pis *= inv;
 final_prior *= inv; mean_prior *= inv; sigma_mean *= inv; sigma_final *= inv;
 final_vb *= inv; final_vg *= inv; final_ve *= inv;
 Rcpp::List cp_out(nt), marker_prior_final_out(nt), alpha_out(nt), alpha_final_out(nt);
 for (int t = 0; t < nt; ++t) {
  cp_out[t] = comp_prob[t] * inv;
  marker_prior_final_out[t] = marker_prior_final[t] * inv;
  alpha_out[t] = alpha_mean[t] * inv;
  alpha_final_out[t] = alpha_final[t] * inv;
 }
 arma::mat wy(m, nt, arma::fill::zeros), r(m, nt, arma::fill::zeros);
 if (return_wy || return_r) {
  for (int t = 0; t < nt; ++t) {
   arma::vec y_t = y_mat.col(t);
   arma::vec residual = y_t - br_xb(G, maps, order, b.col(t).t());
   for (int j = 0; j < m; ++j) {
    if (return_wy) wy(j, t) = br_dot_residual(G, j, maps[j], y_t.memptr());
    if (return_r) r(j, t) = br_dot_residual(G, j, maps[j], residual.memptr());
   }
  }
 }
 Rcpp::CharacterVector component_names(K);
 for (int k = 0; k < K; ++k) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "gamma_%.2f", gamma[k]);
  component_names[k] = std::string(buf);
 }
 component_names[0] = "gamma_0.00";
 arma::mat ncomp(nt, K, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) ncomp.row(t) = arma::sum(comp_prob[t], 0) * inv;
 Rcpp::RObject chains_out = R_NilValue;
 if (keep_chains) {
  Rcpp::List retained_chains(nt);
  for (int t = 0; t < nt; ++t) {
   Rcpp::List trait_chains(nchains);
   for (int ch = 0; ch < nchains; ++ch) {
    const sblr::core::BedBayesRCChainExecutionResult& z = jobs[ch * nt + t];
    trait_chains[ch] = Rcpp::List::create(
     Rcpp::Named("bm")=z.bm, Rcpp::Named("dm")=z.dm,
     Rcpp::Named("b")=z.b, Rcpp::Named("state")=z.state,
     Rcpp::Named("comp_prob")=z.comp_prob,
     Rcpp::Named("alpha")=z.annot_alpha_mean,
     Rcpp::Named("sigmaSqAlpha")=z.annot_sigma_mean,
     Rcpp::Named("pis")=z.pis
    );
   }
   retained_chains[t] = trait_chains;
  }
  chains_out = retained_chains;
 }
 Rcpp::List raw = Rcpp::List::create(
  Rcpp::Named("schema") = Rcpp::List::create(Rcpp::Named("class")="stblr_raw", Rcpp::Named("version")=1),
  Rcpp::Named("meta") = Rcpp::List::create(
   Rcpp::Named("model")="bayesrc", Rcpp::Named("backend")="bed_bayesrc",
   Rcpp::Named("data_level")="individual", Rcpp::Named("prior_type")="annotation_component",
   Rcpp::Named("m")=m, Rcpp::Named("nt")=nt, Rcpp::Named("n_trace")=ntrace,
   Rcpp::Named("nit")=nit, Rcpp::Named("nburn")=nburn, Rcpp::Named("nthin")=nthin,
   Rcpp::Named("nchains")=nchains, Rcpp::Named("keep_chains")=keep_chains,
   Rcpp::Named("n_components")=K, Rcpp::Named("n_annotations")=A.n_cols,
   Rcpp::Named("n_groups")=0, Rcpp::Named("annotations")=true,
   Rcpp::Named("scheduled")=false),
  Rcpp::Named("marker") = Rcpp::List::create(
   Rcpp::Named("bm")=bm, Rcpp::Named("dm")=dm,
   Rcpp::Named("wy")=return_wy ? Rcpp::wrap(wy) : R_NilValue,
   Rcpp::Named("r")=return_r ? Rcpp::wrap(r) : R_NilValue,
   Rcpp::Named("b")=b, Rcpp::Named("state")=state),
  Rcpp::Named("trace") = Rcpp::List::create(
   Rcpp::Named("vbs")=vbs, Rcpp::Named("vgs")=vgs, Rcpp::Named("ves")=ves,
   Rcpp::Named("vle")=vle, Rcpp::Named("vld")=vld, Rcpp::Named("pis")=pis),
  Rcpp::Named("variance") = Rcpp::List::create(
   Rcpp::Named("covb")=arma::diagmat(final_vb), Rcpp::Named("covg")=arma::diagmat(final_vg),
   Rcpp::Named("cove")=arma::diagmat(final_ve), Rcpp::Named("vb")=arma::diagmat(final_vb),
   Rcpp::Named("vg")=arma::diagmat(final_vg), Rcpp::Named("ve")=arma::diagmat(final_ve)),
  Rcpp::Named("pi") = Rcpp::List::create(Rcpp::Named("final")=final_prior,
   Rcpp::Named("mean")=mean_prior, Rcpp::Named("names")=component_names),
  Rcpp::Named("diagnostics") = Rcpp::List::create(
   Rcpp::Named("nsamples")=nsamples*inv, Rcpp::Named("n_used")=n_used,
   Rcpp::Named("log_cpo")=log_cpo*inv, Rcpp::Named("mean_log_cpo")=mean_log_cpo*inv,
   Rcpp::Named("full_sweeps")=true, Rcpp::Named("adaptive_skipping")=false,
   Rcpp::Named("annotation_updates_per_chain")=
    (updateAlpha ? (nit + nburn) / annot_alpha_update_every : 0),
   Rcpp::Named("ld_swap")=R_NilValue),
  Rcpp::Named("chains")=chains_out, Rcpp::Named("prior")=Rcpp::List::create(),
  Rcpp::Named("group")=Rcpp::List::create(),
  Rcpp::Named("annotation") = Rcpp::List::create(
   Rcpp::Named("annotation_names")=R_NilValue, Rcpp::Named("alpha_mean")=alpha_out,
   Rcpp::Named("alpha_final")=alpha_final_out,
   Rcpp::Named("sigmaSqAlpha_mean")=sigma_mean.t(),
   Rcpp::Named("sigmaSqAlpha_final")=sigma_final.t(),
   Rcpp::Named("marker_prior_final")=marker_prior_final_out),
  Rcpp::Named("component") = Rcpp::List::create(
   Rcpp::Named("names")=component_names, Rcpp::Named("mixture_var")=gamma,
   Rcpp::Named("prob")=cp_out, Rcpp::Named("ncomp")=ncomp,
   Rcpp::Named("dm_component_mean")=component_mean),
  Rcpp::Named("selection") = Rcpp::List::create(Rcpp::Named("enabled")=false,
   Rcpp::Named("fixed")=false, Rcpp::Named("trace")=R_NilValue)
 );
 raw.attr("class") = Rcpp::CharacterVector::create("stblr_raw_v1", "stblr_raw", "list");
 return raw;
}
