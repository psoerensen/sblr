// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "blr_csr_bayesr_rcpp_adapter.h"
#include "blr_csr_logvar_bayesr_core_impl.h"

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_logvar_bayesr(
 std::vector<std::vector<double>> wy,
 std::vector<std::vector<double>> ww,
 std::vector<double> yy,
 std::vector<std::vector<double>> b_init,
 std::vector<std::vector<double>> comp_init,
 bool use_comp_init,
 std::vector<std::vector<double>> r_init,
 bool use_r_init,
 bool rebuild_r_before_updateE,
 std::string ld_prefix,
 arma::mat B,
 arma::mat E,
 std::vector<std::vector<double>> ssb_prior,
 std::vector<std::vector<double>> sse_prior,
 std::vector<double> pi,
 std::vector<double> mixture_var,
 std::vector<double> alpha,
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
 int ncores,
 int seed,
 int nchains,
 bool keep_chains,
 std::vector<int> chain_seeds,
 int updateE_start,
 int updateE_every,
 bool updateLDswap,
 double ld_swap_prob,
 double ld_swap_r2,
 int ld_swap_max_friends,
 int ld_swap_moves,
 arma::mat annotation,
 arma::mat theta_init,
 double theta_prior_sd,
 bool updateTheta,
 Rcpp::IntegerVector convergence_markers,
 bool convergence_probability,
 bool convergence_b,
 bool convergence_d,
 bool convergence_component
) {
 const int nt = static_cast<int>(wy.size());
 if (nt <= 0) throw std::invalid_argument("BayesR-LV requires at least one trait.");
 const int m = static_cast<int>(wy[0].size());
 if (m <= 0 || nit <= 0 || nburn < 0 || nthin <= 0 || nchains <= 0 ||
     ncores <= 0)
  throw std::invalid_argument("BayesR-LV MCMC dimensions are invalid.");
 if (annotation.n_rows != static_cast<arma::uword>(m) ||
     annotation.n_cols == 0 || !annotation.is_finite())
  throw std::invalid_argument("BayesR-LV annotation must be a finite m by p matrix.");
 const int p = static_cast<int>(annotation.n_cols);
 if (theta_init.n_rows != static_cast<arma::uword>(p) ||
     theta_init.n_cols != static_cast<arma::uword>(nt) ||
     !theta_init.is_finite())
  throw std::invalid_argument("BayesR-LV theta_init must be finite p by nt.");
 if (!std::isfinite(theta_prior_sd) || theta_prior_sd <= 0.0)
  throw std::invalid_argument("BayesR-LV theta_prior_sd must be positive finite.");
 if (!chain_seeds.empty() && static_cast<int>(chain_seeds.size()) != nchains)
  throw std::invalid_argument("BayesR-LV chain_seeds must have length nchains.");

 const int ntasks = nt * nchains;
 std::vector<sblr::logvar::CsrLogvarBayesRChainOutput> policy_outputs(
  static_cast<std::size_t>(ntasks));
 sblr::logvar::CsrLogvarBayesRPolicyInput policy_input;
 policy_input.annotation = &annotation;
 policy_input.theta_initial = &theta_init;
 policy_input.theta_prior_sd = theta_prior_sd;
 policy_input.update_theta = updateTheta;
 policy_input.trace_count = nit + nburn;
 policy_input.outputs = &policy_outputs;
 sblr::logvar::CsrLogvarBayesRPolicyFactory policy_factory(policy_input);

 Rcpp::List raw = stblr_cpg_omp_csr_bayesr_with_policy(
  std::move(wy), std::move(ww), std::move(yy), std::move(b_init),
  std::move(comp_init), use_comp_init, std::move(r_init), use_r_init,
  rebuild_r_before_updateE, std::move(ld_prefix), std::move(B), std::move(E),
  std::move(ssb_prior), std::move(sse_prior), std::move(pi),
  std::move(mixture_var), std::move(alpha), nub, nue, updateB, updateE,
  updatePi, adjE, std::move(n), nit, nburn, nthin, ncores, seed, nchains,
  keep_chains, std::move(chain_seeds), updateE_start, updateE_every,
  updateLDswap, ld_swap_prob, ld_swap_r2, ld_swap_max_friends,
  ld_swap_moves, R_NilValue, false, 0.0,
  Rcpp::NumericVector::create(-3.0, 2.0), 0.35, R_NilValue,
  convergence_markers, convergence_probability, convergence_b,
  convergence_d, convergence_component, &policy_factory);

 Rcpp::List meta = raw["meta"];
 meta["model"] = "sbayesr_logvar";
 meta["backend"] = "csr_logvar_bayesr";
 meta["prior_type"] = "annotation_log_variance";
 meta["n_annotations"] = p;

 Rcpp::NumericMatrix theta_mean(p, nt);
 Rcpp::NumericMatrix q_mean(m, nt);
 Rcpp::NumericMatrix q_chain_mean(m, ntasks);
 Rcpp::NumericVector theta_trace(
  static_cast<R_xlen_t>(nit + nburn) * p * ntasks);
 theta_trace.attr("dim") = Rcpp::IntegerVector::create(nit + nburn, p, ntasks);
 double updates = 0.0, evaluations = 0.0, contractions = 0.0;
 std::size_t max_evaluations = 0, max_contractions = 0;
 double min_log_q = std::numeric_limits<double>::infinity();
 double max_log_q = -std::numeric_limits<double>::infinity();
 for (int trait = 0; trait < nt; ++trait) for (int chain = 0; chain < nchains; ++chain) {
  const int task = trait * nchains + chain;
  const auto& current = policy_outputs[static_cast<std::size_t>(task)];
  const double retained = std::max(current.retained_samples, 1.0);
  for (int column = 0; column < p; ++column) {
   theta_mean(column, trait) += current.theta_sum(column) /
    retained / static_cast<double>(nchains);
   for (int iteration = 0; iteration < nit + nburn; ++iteration)
    theta_trace[iteration + (nit + nburn) * (column + p * task)] =
     current.theta_trace(iteration, column);
  }
  for (int marker = 0; marker < m; ++marker) {
   const double chain_q = current.prior_scale_sum(marker) / retained;
   q_chain_mean(marker, task) = chain_q;
   q_mean(marker, trait) += chain_q / static_cast<double>(nchains);
  }
  updates += current.diagnostics.theta_updates;
  evaluations += current.diagnostics.likelihood_evaluations;
  contractions += current.diagnostics.bracket_contractions;
  max_evaluations = std::max(max_evaluations,
   current.diagnostics.max_likelihood_evaluations);
  max_contractions = std::max(max_contractions,
   current.diagnostics.max_bracket_contractions);
  min_log_q = std::min(min_log_q, current.diagnostics.min_log_q);
  max_log_q = std::max(max_log_q, current.diagnostics.max_log_q);
 }
 Rcpp::NumericMatrix variance_ratio(p, nt);
 for (int trait = 0; trait < nt; ++trait) for (int column = 0; column < p; ++column)
  variance_ratio(column, trait) = std::exp(theta_mean(column, trait));
 raw["annotation"] = Rcpp::List::create(
  Rcpp::Named("theta") = theta_mean,
  Rcpp::Named("theta_trace") = theta_trace,
  Rcpp::Named("variance_ratio") = variance_ratio,
  Rcpp::Named("marker_prior_scale") = q_mean,
  Rcpp::Named("marker_prior_scale_chain") = q_chain_mean,
  Rcpp::Named("theta_prior_sd") = theta_prior_sd,
  Rcpp::Named("update_theta") = updateTheta);
 Rcpp::List diagnostics = raw["diagnostics"];
 diagnostics["logvar"] = Rcpp::List::create(
  Rcpp::Named("theta_updates") = updates,
  Rcpp::Named("mean_likelihood_evaluations_per_update") =
   updates > 0.0 ? evaluations / updates : 0.0,
  Rcpp::Named("max_likelihood_evaluations") = max_evaluations,
  Rcpp::Named("mean_bracket_contractions") =
   updates > 0.0 ? contractions / updates : 0.0,
  Rcpp::Named("max_bracket_contractions") = max_contractions,
  Rcpp::Named("min_log_q") = min_log_q,
  Rcpp::Named("max_log_q") = max_log_q);
 raw["diagnostics"] = diagnostics;
 return raw;
}
