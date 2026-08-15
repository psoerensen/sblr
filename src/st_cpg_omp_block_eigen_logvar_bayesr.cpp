// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "blr_block_eigen_rcpp_adapter.h"
#include "blr_csr_logvar_bayesr_core_impl.h"
#include "st_logvar_annotation_rcpp.h"

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_block_eigen_logvar_bayesr(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> comp_init,
  bool use_comp_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
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
  Rcpp::IntegerVector convergence_markers,
  bool convergence_probability,
  bool convergence_b,
  bool convergence_d,
  bool convergence_component,
  Rcpp::CharacterVector bed_files,
  int n_bed,
  Rcpp::List cls,
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::NumericVector af,
  Rcpp::IntegerVector block_start,
  std::string eigen_filter,
  double eigen_tau,
  double eigen_eta,
  std::string representation,
  double eigen_prop,
  int low_rank_residual_rebuild_every,
  Rcpp::List block_residual_config,
  arma::mat annotation,
  arma::mat theta_init,
  double theta_prior_sd,
  bool updateTheta
) {
  const int nt = static_cast<int>(wy.size());
  if (nt <= 0) {
    throw std::invalid_argument("BayesR-LV block fit requires at least one trait.");
  }
  const int m = static_cast<int>(wy[0].size());
  if (m <= 0 || annotation.n_rows != static_cast<arma::uword>(m) ||
      annotation.n_cols == 0 || !annotation.is_finite()) {
    throw std::invalid_argument(
      "BayesR-LV block annotation must be a finite m by p matrix.");
  }
  const int p = static_cast<int>(annotation.n_cols);
  if (theta_init.n_rows != static_cast<arma::uword>(p) ||
      theta_init.n_cols != static_cast<arma::uword>(nt) ||
      !theta_init.is_finite()) {
    throw std::invalid_argument(
      "BayesR-LV block theta_init must be finite p by nt.");
  }
  if (!std::isfinite(theta_prior_sd) || theta_prior_sd <= 0.0) {
    throw std::invalid_argument(
      "BayesR-LV block theta_prior_sd must be positive finite.");
  }

  const int task_count = nt * nchains;
  std::vector<sblr::logvar::CsrLogvarBayesRChainOutput> outputs(
    static_cast<std::size_t>(task_count));
  sblr::logvar::CsrLogvarBayesRPolicyInput policy_input;
  policy_input.annotation = &annotation;
  policy_input.theta_initial = &theta_init;
  policy_input.theta_prior_sd = theta_prior_sd;
  policy_input.update_theta = updateTheta;
  policy_input.trace_count = nit + nburn;
  policy_input.outputs = &outputs;
  sblr::logvar::CsrLogvarBayesRPolicyFactory policy_factory(policy_input);

  Rcpp::List raw = stblr_cpg_omp_csr_bayesr_block_eigen_with_policy(
    std::move(wy), std::move(ww), std::move(yy), std::move(b_init),
    std::move(comp_init), use_comp_init, std::move(r_init), use_r_init,
    rebuild_r_before_updateE, "", std::move(B), std::move(E),
    std::move(ssb_prior), std::move(sse_prior), std::move(pi),
    std::move(mixture_var), std::move(alpha), nub, nue, updateB, updateE,
    updatePi, adjE, std::move(n), nit, nburn, nthin, ncores, seed, nchains,
    keep_chains, std::move(chain_seeds), updateE_start, updateE_every,
    false, 0.05, 0.8, 50, 1, R_NilValue, false, 0.0,
    Rcpp::NumericVector::create(-3.0, 2.0), 0.35, R_NilValue,
    convergence_markers, convergence_probability, convergence_b,
    convergence_d, convergence_component, bed_files, n_bed, cls, rows, af,
    block_start, std::move(eigen_filter), eigen_tau, eigen_eta,
    std::move(representation), eigen_prop, low_rank_residual_rebuild_every,
    block_residual_config, &policy_factory, BlrPhase3ExecutionContract()
  );
  return decorate_logvar_raw(
    raw, outputs, m, nt, p, nit + nburn, nchains, theta_prior_sd,
    updateTheta, "sbayesr_logvar", "block_eigen_logvar_bayesr"
  );
}
