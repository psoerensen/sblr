// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "st_bed_bayesr_common.h"
#include "st_bayesrc_annotation_prior.h"

struct ChainResultBayesRC {
 arma::rowvec bm, dm, component_mean, b, state;
 arma::rowvec vbs, vgs, ves, vles, vlds, pis;
 arma::mat comp_prob;
 arma::mat annot_alpha_mean, annot_alpha_final;
 arma::vec annot_sigma_mean, annot_sigma_final;
 arma::rowvec mean_prior;
 arma::vec residual;
 double final_vb = 0.0, final_vg = 0.0, final_ve = 0.0;
 double log_cpo = NA_REAL, mean_log_cpo = NA_REAL, nsamples = 0.0;
 int failed = 0;
 std::string error;
};

static inline void sample_marker_bayesrc(
  const FastPackedBedMatrixBR& G, int marker, const MarkerMapBayesR& map,
  const arma::rowvec& marker_prior, const std::vector<double>& gamma,
  double vb, double vei, arma::vec& residual, double& marker_effect,
  int& component, std::mt19937& gen
) {
 static thread_local std::uniform_real_distribution<double> runif(0.0, 1.0);
 static thread_local std::normal_distribution<double> norm01(0.0, 1.0);
 const int ncomponent = static_cast<int>(gamma.size());
 const double vei_safe = std::max(vei, 1e-300);
 const double score = br_dot_residual(G, marker, map, residual.memptr()) +
  map.xx * marker_effect;
 std::vector<double> logp(static_cast<std::size_t>(ncomponent));
 logp[0] = std::log(std::max(marker_prior(0), 1e-300));
 for (int k = 1; k < ncomponent; ++k) {
  const double vbk = vb * gamma[static_cast<std::size_t>(k)];
  const double denom = std::max(vei_safe + map.xx * vbk, 1e-300);
  logp[static_cast<std::size_t>(k)] =
   std::log(std::max(marker_prior(static_cast<arma::uword>(k)), 1e-300)) +
   0.5 * std::log(vei_safe / denom) +
   0.5 * score * score * vbk / (vei_safe * denom);
 }
 const double mx = *std::max_element(logp.begin(), logp.end());
 std::vector<double> prob(static_cast<std::size_t>(ncomponent));
 double psum = 0.0;
 for (int k = 0; k < ncomponent; ++k) {
  prob[static_cast<std::size_t>(k)] = std::exp(logp[static_cast<std::size_t>(k)] - mx);
  psum += prob[static_cast<std::size_t>(k)];
 }
 int component_new = ncomponent - 1;
 const double u = runif(gen);
 double cumulative = 0.0;
 for (int k = 0; k < ncomponent - 1; ++k) {
  cumulative += prob[static_cast<std::size_t>(k)] / psum;
  if (u < cumulative) { component_new = k; break; }
 }
 double effect_new = 0.0;
 if (component_new > 0) {
  const double vbk = vb * gamma[static_cast<std::size_t>(component_new)];
  const double lhs = map.xx + vei_safe / vbk;
  effect_new = score / lhs + std::sqrt(vei_safe / lhs) * norm01(gen);
 }
 const double diff = effect_new - marker_effect;
 if (diff != 0.0)
  br_update_residual(G, marker, map, residual.memptr(), diff);
 marker_effect = effect_new;
 component = component_new;
}

static ChainResultBayesRC run_one_bayesrc_chain(
  const FastPackedBedMatrixBR& G,
  const std::vector<MarkerMapBayesR>& maps,
  const std::vector<int>& marker_order, const arma::mat& y,
  const std::vector<std::vector<double>>& b_init, const arma::mat& B,
  const arma::mat& E, const arma::mat& ssb_prior, const arma::mat& sse_prior,
  const arma::mat& annotation, const std::vector<double>& gamma,
  const arma::mat& annot_alpha_init, const arma::vec& annot_sigma_init,
  bool intercept_flat, double sigma_a, double sigma_b, double pi_floor,
  double nub, double nue, bool update_alpha, bool update_b, bool update_e,
  int alpha_every, double adjE, int nit, int nburn, int nthin,
  int rebuild_every, int trait, int chain, int seed
) {
 const int m = G.m, K = static_cast<int>(gamma.size());
 const int total_it = nit + nburn;
 ChainResultBayesRC out;
 out.bm.zeros(m); out.dm.zeros(m); out.component_mean.zeros(m);
 out.b.zeros(m); out.state.zeros(m);
 out.vbs.zeros(total_it); out.vgs.zeros(total_it); out.ves.zeros(total_it);
 out.vles.zeros(total_it); out.vlds.zeros(total_it); out.pis.zeros(total_it);
 out.comp_prob.zeros(m, K); out.mean_prior.zeros(K);
 try {
  std::mt19937 gen(static_cast<unsigned int>(
   seed + 1000003 * (trait + 1) + 9176 * (chain + 1)
  ));
  arma::vec y_t = y.col(static_cast<arma::uword>(trait));
  arma::rowvec b_t(m, arma::fill::zeros);
  arma::Row<int> component_t(m, arma::fill::zeros);
  for (int j = 0; j < m; ++j) {
   b_t(static_cast<arma::uword>(j)) =
    b_init[static_cast<std::size_t>(trait)][static_cast<std::size_t>(j)];
   component_t(static_cast<arma::uword>(j)) = b_t(static_cast<arma::uword>(j)) != 0.0 ? 1 : 0;
  }
  arma::vec residual = y_t - br_xb(G, maps, marker_order, b_t);
  double vb = B(trait, trait), ve = E(trait, trait);
  double vg = br_computeG(y_t, residual), vei = ve + adjE * vg;
  arma::mat annot_alpha = annot_alpha_init;
  arma::vec annot_sigma = annot_sigma_init;
  arma::mat snp_pi = st_bayesrc_compute_snp_pi(annotation, annot_alpha, pi_floor);
  arma::mat alpha_acc(annotation.n_cols, K - 1, arma::fill::zeros);
  arma::vec sigma_acc(K - 1, arma::fill::zeros);
  arma::vec log_inv_cpo(G.n, arma::fill::value(-std::numeric_limits<double>::infinity()));
  double nsamples = 0.0;

  for (int it = 0; it < total_it; ++it) {
   // Exact full sweep: every marker is visited once in marker_order every iteration.
   for (int marker : marker_order) {
    const arma::uword ju = static_cast<arma::uword>(marker);
    double effect = b_t(ju);
    int component = component_t(ju);
    sample_marker_bayesrc(
     G, marker, maps[static_cast<std::size_t>(marker)], snp_pi.row(ju), gamma,
     vb, vei, residual, effect, component, gen
    );
    b_t(ju) = effect;
    component_t(ju) = component;
   }
   if (update_e && rebuild_every > 0 && ((it + 1) % rebuild_every == 0))
    residual = y_t - br_xb(G, maps, marker_order, b_t);
   if (update_b)
    sampleB_bayesr(m, gamma, nub, vb, b_t, component_t, ssb_prior(trait, trait), gen);
   if (update_e)
    sampleE_bayesr(nue, ve, residual, sse_prior(trait, trait), gen);
   if (update_alpha && ((it + 1) % alpha_every == 0)) {
    st_bayesrc_update_annotation_prior(
     annotation, component_t, annot_alpha, annot_sigma, intercept_flat,
     sigma_a, sigma_b, gen
    );
    snp_pi = st_bayesrc_compute_snp_pi(annotation, annot_alpha, pi_floor);
   }
   vg = br_computeG(y_t, residual);
   const double vle = br_computeLE(b_t, maps, G.n);
   vei = ve + adjE * vg;
   out.vbs(static_cast<arma::uword>(it)) = vb;
   out.vgs(static_cast<arma::uword>(it)) = vg;
   out.ves(static_cast<arma::uword>(it)) = ve;
   out.vles(static_cast<arma::uword>(it)) = vle;
   out.vlds(static_cast<arma::uword>(it)) = vg - vle;
   out.pis(static_cast<arma::uword>(it)) = 1.0 - arma::mean(snp_pi.col(0));
   if (it >= nburn && ((it - nburn) % nthin == 0)) {
    nsamples += 1.0;
    br_update_log_inv_cpo(residual, ve, log_inv_cpo);
    for (int j = 0; j < m; ++j) {
     const arma::uword ju = static_cast<arma::uword>(j);
     out.bm(ju) += b_t(ju);
     out.dm(ju) += component_t(ju) > 0 ? 1.0 : 0.0;
     out.component_mean(ju) += component_t(ju);
     out.comp_prob(ju, static_cast<arma::uword>(component_t(ju))) += 1.0;
    }
    alpha_acc += annot_alpha;
    sigma_acc += annot_sigma;
    out.mean_prior += arma::mean(snp_pi, 0);
   }
  }
  if (nsamples <= 0.0) nsamples = 1.0;
  out.bm /= nsamples; out.dm /= nsamples; out.component_mean /= nsamples;
  out.comp_prob /= nsamples; alpha_acc /= nsamples; sigma_acc /= nsamples;
  out.mean_prior /= nsamples;
  out.b = b_t;
  for (int j = 0; j < m; ++j) out.state(j) = component_t(j);
  out.annot_alpha_mean = alpha_acc; out.annot_alpha_final = annot_alpha;
  out.annot_sigma_mean = sigma_acc; out.annot_sigma_final = annot_sigma;
  out.residual = residual; out.final_vb = vb; out.final_vg = vg; out.final_ve = ve;
  out.nsamples = nsamples;
  out.log_cpo = br_compute_total_log_cpo(log_inv_cpo, nsamples);
  out.mean_log_cpo = out.log_cpo / G.n;
 } catch (const std::exception& ex) { out.failed = 1; out.error = ex.what(); }
 return out;
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
 std::vector<ChainResultBayesRC> jobs(static_cast<std::size_t>(njobs));
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores) schedule(static)
#endif
 for (int job = 0; job < njobs; ++job) {
  jobs[job] = run_one_bayesrc_chain(
   G, maps, order, y_mat, b_init, B, E, ssb, sse, A, gamma,
   annot_alpha_init, annot_sigma_sq_alpha_init, intercept_flat,
   sigmaSqAlpha_a, sigmaSqAlpha_b, pi_floor, nub, nue, updateAlpha,
   updateB, updateE, annot_alpha_update_every, adjE, nit, nburn, nthin,
   rebuild_every, job % nt, job / nt, seed
  );
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
   const ChainResultBayesRC& z = jobs[ch * nt + t];
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
 Rcpp::List chains_out = R_NilValue;
 if (keep_chains) {
  chains_out = Rcpp::List(nt);
  for (int t = 0; t < nt; ++t) {
   Rcpp::List trait_chains(nchains);
   for (int ch = 0; ch < nchains; ++ch) {
    const ChainResultBayesRC& z = jobs[ch * nt + t];
    trait_chains[ch] = Rcpp::List::create(
     Rcpp::Named("bm")=z.bm, Rcpp::Named("dm")=z.dm,
     Rcpp::Named("b")=z.b, Rcpp::Named("state")=z.state,
     Rcpp::Named("comp_prob")=z.comp_prob,
     Rcpp::Named("alpha")=z.annot_alpha_mean,
     Rcpp::Named("sigmaSqAlpha")=z.annot_sigma_mean,
     Rcpp::Named("pis")=z.pis
    );
   }
   chains_out[t] = trait_chains;
  }
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
