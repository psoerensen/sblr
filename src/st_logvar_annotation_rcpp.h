#ifndef SBLR_ST_LOGVAR_ANNOTATION_RCPP_H
#define SBLR_ST_LOGVAR_ANNOTATION_RCPP_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

template <class ChainOutput>
inline Rcpp::List decorate_logvar_raw(
  Rcpp::List raw,
  const std::vector<ChainOutput>& outputs,
  int marker_count,
  int trait_count,
  int annotation_count,
  int trace_count,
  int chain_count,
  double theta_prior_sd,
  bool update_theta,
  const std::string& model,
  const std::string& backend
) {
  Rcpp::List meta = raw["meta"];
  meta["model"] = model;
  meta["backend"] = backend;
  meta["prior_type"] = "annotation_log_variance";
  meta["n_annotations"] = annotation_count;

  const int task_count = trait_count * chain_count;
  Rcpp::NumericMatrix theta_mean(annotation_count, trait_count);
  Rcpp::NumericMatrix q_mean(marker_count, trait_count);
  Rcpp::NumericMatrix q_chain_mean(marker_count, task_count);
  Rcpp::NumericVector theta_trace(
    static_cast<R_xlen_t>(trace_count) * annotation_count * task_count
  );
  theta_trace.attr("dim") = Rcpp::IntegerVector::create(
    trace_count, annotation_count, task_count
  );
  double updates = 0.0;
  double evaluations = 0.0;
  double contractions = 0.0;
  std::size_t max_evaluations = 0;
  std::size_t max_contractions = 0;
  double min_log_q = std::numeric_limits<double>::infinity();
  double max_log_q = -std::numeric_limits<double>::infinity();
  for (int trait = 0; trait < trait_count; ++trait) {
    for (int chain = 0; chain < chain_count; ++chain) {
      const int task = trait * chain_count + chain;
      const ChainOutput& current = outputs[static_cast<std::size_t>(task)];
      const double retained = std::max(current.retained_samples, 1.0);
      for (int column = 0; column < annotation_count; ++column) {
        theta_mean(column, trait) += current.theta_sum(column) /
          retained / static_cast<double>(chain_count);
        for (int iteration = 0; iteration < trace_count; ++iteration) {
          theta_trace[iteration + trace_count *
            (column + annotation_count * task)] =
            current.theta_trace(iteration, column);
        }
      }
      for (int marker = 0; marker < marker_count; ++marker) {
        const double chain_q = current.prior_scale_sum(marker) / retained;
        q_chain_mean(marker, task) = chain_q;
        q_mean(marker, trait) += chain_q / static_cast<double>(chain_count);
      }
      updates += current.diagnostics.theta_updates;
      evaluations += current.diagnostics.likelihood_evaluations;
      contractions += current.diagnostics.bracket_contractions;
      max_evaluations = std::max(
        max_evaluations, current.diagnostics.max_likelihood_evaluations);
      max_contractions = std::max(
        max_contractions, current.diagnostics.max_bracket_contractions);
      min_log_q = std::min(min_log_q, current.diagnostics.min_log_q);
      max_log_q = std::max(max_log_q, current.diagnostics.max_log_q);
    }
  }
  Rcpp::NumericMatrix variance_ratio(annotation_count, trait_count);
  for (int trait = 0; trait < trait_count; ++trait) {
    for (int column = 0; column < annotation_count; ++column) {
      variance_ratio(column, trait) = std::exp(theta_mean(column, trait));
    }
  }
  raw["annotation"] = Rcpp::List::create(
    Rcpp::Named("theta") = theta_mean,
    Rcpp::Named("theta_trace") = theta_trace,
    Rcpp::Named("variance_ratio") = variance_ratio,
    Rcpp::Named("marker_prior_scale") = q_mean,
    Rcpp::Named("marker_prior_scale_chain") = q_chain_mean,
    Rcpp::Named("theta_prior_sd") = theta_prior_sd,
    Rcpp::Named("update_theta") = update_theta
  );
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
    Rcpp::Named("max_log_q") = max_log_q
  );
  raw["diagnostics"] = diagnostics;
  return raw;
}

#endif
