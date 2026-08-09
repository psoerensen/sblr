// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <random>
#include <string>

#include "st_logvar_annotation_prior.h"

// Internal bindings below exist solely to qualify the shared mathematical
// kernel independently of any CSR or retained-block residual operator.

// [[Rcpp::export]]
Rcpp::List st_logvar_eta_q_internal(
  const arma::mat& annotation,
  const arma::vec& theta
) {
  sblr::logvar::EssUpdateDiagnostics diagnostics;
  const arma::vec eta = sblr::logvar::calculate_eta(annotation, theta);
  const arma::vec q = sblr::logvar::calculate_prior_scale_from_eta(
    eta, &diagnostics);
  return Rcpp::List::create(
    Rcpp::Named("eta") = eta,
    Rcpp::Named("q") = q,
    Rcpp::Named("min_log_q") = diagnostics.min_log_q,
    Rcpp::Named("max_log_q") = diagnostics.max_log_q);
}

// [[Rcpp::export]]
double st_logvar_loglik_bayesc_internal(
  const arma::vec& theta,
  const arma::mat& annotation,
  const arma::vec& effect,
  const arma::ivec& state,
  double marker_variance
) {
  return sblr::logvar::theta_log_likelihood_bayesc(
    theta, annotation, effect, state, marker_variance);
}

// [[Rcpp::export]]
double st_logvar_loglik_bayesr_internal(
  const arma::vec& theta,
  const arma::mat& annotation,
  const arma::vec& effect,
  const arma::ivec& component,
  double marker_variance,
  const arma::vec& gamma
) {
  return sblr::logvar::theta_log_likelihood_bayesr(
    theta, annotation, effect, component, marker_variance, gamma);
}

// [[Rcpp::export]]
Rcpp::List st_logvar_ess_fixture_internal(
  const arma::vec& theta,
  const arma::mat& annotation,
  const arma::vec& effect,
  const arma::ivec& state,
  double marker_variance,
  double theta_prior_sd = 0.7,
  int updates = 1,
  int seed = 1,
  bool empty_active_set = false,
  std::string model = "bayesc",
  Rcpp::Nullable<Rcpp::NumericVector> gamma = R_NilValue
) {
  if (updates <= 0) Rcpp::stop("updates must be positive.");
  std::mt19937 generator(static_cast<std::mt19937::result_type>(seed));
  sblr::logvar::EssUpdateDiagnostics diagnostics;
  arma::vec current = theta;
  arma::vec gamma_value;
  if (model == "bayesr") {
    if (gamma.isNull()) Rcpp::stop("gamma is required for BayesR-LV ESS.");
    gamma_value = Rcpp::as<arma::vec>(gamma);
  } else if (model != "bayesc") {
    Rcpp::stop("model must be 'bayesc' or 'bayesr'.");
  }

  arma::mat draws(static_cast<arma::uword>(updates), theta.n_elem);
  for (int iteration = 0; iteration < updates; ++iteration) {
    if (model == "bayesc") {
      const auto log_likelihood = [&](const arma::vec& value) {
        return sblr::logvar::theta_log_likelihood_bayesc(
          value, annotation, effect, state, marker_variance);
      };
      current = sblr::logvar::elliptical_slice_update(
        current, theta_prior_sd, empty_active_set, log_likelihood,
        generator, diagnostics);
    } else {
      const auto log_likelihood = [&](const arma::vec& value) {
        return sblr::logvar::theta_log_likelihood_bayesr(
          value, annotation, effect, state, marker_variance, gamma_value);
      };
      current = sblr::logvar::elliptical_slice_update(
        current, theta_prior_sd, empty_active_set, log_likelihood,
        generator, diagnostics);
    }
    draws.row(static_cast<arma::uword>(iteration)) = current.t();
  }

  const arma::vec eta = sblr::logvar::calculate_eta(annotation, current);
  sblr::logvar::calculate_prior_scale_from_eta(eta, &diagnostics);
  const double update_count = static_cast<double>(diagnostics.theta_updates);
  return Rcpp::List::create(
    Rcpp::Named("theta") = current,
    Rcpp::Named("draws") = draws,
    Rcpp::Named("theta_updates") = diagnostics.theta_updates,
    Rcpp::Named("mean_likelihood_evaluations_per_update") =
      update_count > 0.0 ? diagnostics.likelihood_evaluations / update_count : 0.0,
    Rcpp::Named("max_likelihood_evaluations") =
      diagnostics.max_likelihood_evaluations,
    Rcpp::Named("mean_bracket_contractions") =
      update_count > 0.0 ? diagnostics.bracket_contractions / update_count : 0.0,
    Rcpp::Named("max_bracket_contractions") =
      diagnostics.max_bracket_contractions,
    Rcpp::Named("min_log_q") = diagnostics.min_log_q,
    Rcpp::Named("max_log_q") = diagnostics.max_log_q);
}
