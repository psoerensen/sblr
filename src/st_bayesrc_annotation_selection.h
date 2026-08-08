#ifndef ST_BAYESRC_ANNOTATION_SELECTION_H
#define ST_BAYESRC_ANNOTATION_SELECTION_H

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include <armadillo>

#include "st_bayesrc_annotation_prior.h"

// SBayesRC-S model-specific annotation hierarchy.  This header deliberately
// owns no marker likelihood, LD operator, residual state, or production model
// dispatch.  It composes the validated neutral BayesRC probit helpers.

struct StBayesRCSelectionState {
 arma::uvec delta;   // J selectable non-intercept annotations
 arma::mat alpha;    // (J + 1) x K; row zero is the always-included intercept
 double pi_a;
 arma::vec tau2;     // K stick-specific slab variances
};

struct StBayesRCSelectionHyperParameters {
 double a_pi;
 double b_pi;
 double a_tau;
 double b_tau;
};

// Optional policy object consumed by the existing CSR BayesRC genomic engine.
// A disabled instance is inert and must not consume RNG or alter the standard
// SBayesRC path.  The annotation matrix passed to that engine contains the
// intercept in column zero; delta_init therefore has one fewer element.
struct StBayesRCSelectionGenomicConfig {
 bool enabled = false;
 arma::uvec delta_init;
 double pi_a_init = 0.5;
 arma::vec tau2_init;
 StBayesRCSelectionHyperParameters hyper{1.0, 1.0, 3.0, 1.6};
 arma::vec intercept_mean;
 arma::vec intercept_precision;
 bool fixed_delta = false;
 arma::ivec fixed_delta_value;
};

struct StBayesRCSelectionMoments {
 double s;
 double t;
 double log_bf;
 double mean;
 double variance;
};

inline arma::vec st_bayesrc_selection_column_rows(
 const arma::mat& matrix,
 arma::uword column,
 const arma::uvec& rows
) {
 arma::vec result(rows.n_elem);
 for (arma::uword i = 0; i < rows.n_elem; ++i) {
  result(i) = matrix(rows(i), column);
 }
 return result;
}

inline StBayesRCSelectionMoments st_bayesrc_selection_moments(
 const arma::vec& x,
 const arma::vec& residual,
 double tau2
) {
 if (x.n_elem != residual.n_elem || !x.is_finite() ||
     !residual.is_finite() || !std::isfinite(tau2) || tau2 <= 0.0) {
  throw std::invalid_argument("invalid SBayesRC-S conditional inputs");
 }
 const double s = arma::dot(x, x);
 const double t = arma::dot(x, residual);
 const double denominator = 1.0 + tau2 * s;
 const double variance = tau2 / denominator;
 const double mean = variance * t;
 const double log_bf = -0.5 * std::log1p(tau2 * s) +
  0.5 * tau2 * t * t / denominator;
 if (!std::isfinite(log_bf) || !std::isfinite(mean) ||
     !std::isfinite(variance) || variance <= 0.0) {
  throw std::runtime_error("invalid SBayesRC-S conditional moments");
 }
 return StBayesRCSelectionMoments{s, t, log_bf, mean, variance};
}

inline arma::vec st_bayesrc_selection_beta_parameters(
 const arma::uvec& delta,
 double a_pi,
 double b_pi
) {
 if (!delta.is_finite() || arma::any(delta > 1u) ||
     !std::isfinite(a_pi) || a_pi <= 0.0 ||
     !std::isfinite(b_pi) || b_pi <= 0.0) {
  throw std::invalid_argument("invalid SBayesRC-S beta parameters");
 }
 const double included = static_cast<double>(arma::accu(delta));
 return arma::vec({a_pi + included,
  b_pi + static_cast<double>(delta.n_elem) - included});
}

inline arma::mat st_bayesrc_selection_ig_parameters(
 const arma::mat& alpha,
 const arma::uvec& delta,
 double a_tau,
 double b_tau
) {
 if (alpha.n_rows != delta.n_elem + 1u || !alpha.is_finite() ||
     !delta.is_finite() || arma::any(delta > 1u) ||
     !std::isfinite(a_tau) || a_tau <= 0.0 ||
     !std::isfinite(b_tau) || b_tau <= 0.0) {
  throw std::invalid_argument("invalid SBayesRC-S inverse-gamma inputs");
 }
 arma::mat result(alpha.n_cols, 2u, arma::fill::zeros);
 const double shape = a_tau + 0.5 * static_cast<double>(arma::accu(delta));
 for (arma::uword stick = 0; stick < alpha.n_cols; ++stick) {
  double sum_squares = 0.0;
  for (arma::uword annotation = 0; annotation < delta.n_elem; ++annotation) {
   if (delta(annotation) == 1u) {
    const double value = alpha(annotation + 1u, stick);
    sum_squares += value * value;
   }
  }
  result(stick, 0u) = shape;
  result(stick, 1u) = b_tau + 0.5 * sum_squares;
 }
 return result;
}

inline double st_bayesrc_selection_draw_beta(
 double shape_a,
 double shape_b,
 std::mt19937& generator
) {
 std::gamma_distribution<double> gamma_a(shape_a, 1.0);
 std::gamma_distribution<double> gamma_b(shape_b, 1.0);
 const double x = gamma_a(generator);
 const double y = gamma_b(generator);
 const double result = x / (x + y);
 if (!std::isfinite(result) || result <= 0.0 || result >= 1.0) {
  throw std::runtime_error("SBayesRC-S beta update produced invalid pi_A");
 }
 return result;
}

inline double st_bayesrc_selection_draw_ig(
 double shape,
 double scale,
 std::mt19937& generator
) {
 if (!std::isfinite(shape) || shape <= 0.0 ||
     !std::isfinite(scale) || scale <= 0.0) {
  throw std::invalid_argument("invalid SBayesRC-S inverse-gamma parameters");
 }
 std::gamma_distribution<double> precision(shape, 1.0 / scale);
 const double draw = 1.0 / std::max(
  precision(generator), std::numeric_limits<double>::min());
 if (!std::isfinite(draw) || draw <= 0.0) {
  throw std::runtime_error("SBayesRC-S inverse-gamma update failed");
 }
 return draw;
}

inline void st_bayesrc_selection_validate(
 const arma::mat& annotation,
 const std::vector<arma::uvec>& eligible,
 const std::vector<arma::ivec>* outcome,
 const StBayesRCSelectionState& state,
 const StBayesRCSelectionHyperParameters& hyper,
 const arma::vec& intercept_mean,
 const arma::vec& intercept_precision
) {
 const arma::uword annotations = annotation.n_cols;
 const arma::uword sticks = eligible.size();
 if (!annotation.is_finite() || annotations == 0u || sticks == 0u ||
     state.delta.n_elem != annotations || arma::any(state.delta > 1u) ||
     state.alpha.n_rows != annotations + 1u || state.alpha.n_cols != sticks ||
     !state.alpha.is_finite() || state.tau2.n_elem != sticks ||
     !state.tau2.is_finite() || arma::any(state.tau2 <= 0.0) ||
     !std::isfinite(state.pi_a) || state.pi_a <= 0.0 || state.pi_a >= 1.0 ||
     !std::isfinite(hyper.a_pi) || hyper.a_pi <= 0.0 ||
     !std::isfinite(hyper.b_pi) || hyper.b_pi <= 0.0 ||
     !std::isfinite(hyper.a_tau) || hyper.a_tau <= 0.0 ||
     !std::isfinite(hyper.b_tau) || hyper.b_tau <= 0.0 ||
     intercept_mean.n_elem != sticks || !intercept_mean.is_finite() ||
     intercept_precision.n_elem != sticks ||
     !intercept_precision.is_finite() ||
     arma::any(intercept_precision <= 0.0)) {
  throw std::invalid_argument("invalid SBayesRC-S hierarchy dimensions/state");
 }
 if (outcome != nullptr && outcome->size() != sticks) {
  throw std::invalid_argument("SBayesRC-S outcome stick count mismatch");
 }
 for (arma::uword stick = 0; stick < sticks; ++stick) {
  const arma::uvec& rows = eligible[static_cast<std::size_t>(stick)];
  if (rows.n_elem > 0u && arma::any(rows >= annotation.n_rows)) {
   throw std::invalid_argument("SBayesRC-S eligible rows are invalid");
  }
  if (outcome != nullptr) {
   const arma::ivec& d = (*outcome)[static_cast<std::size_t>(stick)];
   if (d.n_elem != rows.n_elem || arma::any(d < 0) || arma::any(d > 1)) {
    throw std::invalid_argument("SBayesRC-S outcome dimensions are invalid");
   }
  }
 }
}

inline arma::mat st_bayesrc_selection_compute_q(
 const arma::mat& annotation,
 const arma::mat& alpha
) {
 if (alpha.n_rows != annotation.n_cols + 1u || !annotation.is_finite() ||
     !alpha.is_finite()) {
  throw std::invalid_argument("invalid SBayesRC-S q inputs");
 }
 arma::mat design(annotation.n_rows, annotation.n_cols + 1u,
                  arma::fill::ones);
 design.cols(1u, annotation.n_cols) = annotation;
 arma::mat q(annotation.n_rows, alpha.n_cols, arma::fill::zeros);
 const arma::mat eta = design * alpha;
 for (arma::uword i = 0; i < eta.n_rows; ++i) {
  for (arma::uword stick = 0; stick < eta.n_cols; ++stick) {
   q(i, stick) = st_bayesrc_safe_pnorm(eta(i, stick));
  }
 }
 return q;
}

inline arma::mat st_bayesrc_selection_compute_pi(
 const arma::mat& annotation,
 const arma::mat& alpha,
 double probability_floor
) {
 arma::mat design(annotation.n_rows, annotation.n_cols + 1u,
                  arma::fill::ones);
 design.cols(1u, annotation.n_cols) = annotation;
 return st_bayesrc_compute_snp_pi(design, alpha, probability_floor);
}

inline void st_bayesrc_selection_build_observed_sticks(
 const arma::Row<int>& component,
 int stick_count,
 std::vector<arma::uvec>& eligible,
 std::vector<arma::ivec>& outcome
) {
 if (stick_count <= 0) throw std::invalid_argument("invalid SBayesRC-S stick count");
 arma::Mat<int> indicators(component.n_elem, static_cast<arma::uword>(stick_count),
                           arma::fill::zeros);
 st_bayesrc_build_step_indicators(component, indicators);
 eligible.resize(static_cast<std::size_t>(stick_count));
 outcome.resize(static_cast<std::size_t>(stick_count));
 for (int stick = 0; stick < stick_count; ++stick) {
  std::vector<arma::uword> rows;
  rows.reserve(component.n_elem);
  for (arma::uword marker = 0; marker < component.n_elem; ++marker) {
   if (stick == 0 || indicators(marker, static_cast<arma::uword>(stick - 1)) > 0)
    rows.push_back(marker);
  }
  eligible[static_cast<std::size_t>(stick)] = arma::conv_to<arma::uvec>::from(rows);
  arma::ivec d(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i)
   d(static_cast<arma::uword>(i)) = indicators(rows[i], static_cast<arma::uword>(stick));
  outcome[static_cast<std::size_t>(stick)] = std::move(d);
 }
}

inline std::vector<arma::vec> st_bayesrc_selection_sample_latent(
 const arma::mat& annotation,
 const std::vector<arma::uvec>& eligible,
 const std::vector<arma::ivec>& outcome,
 const arma::mat& alpha,
 std::mt19937& generator
) {
 std::vector<arma::vec> latent(eligible.size());
 for (std::size_t stick = 0; stick < eligible.size(); ++stick) {
  const arma::uvec& rows = eligible[stick];
  arma::vec eta(rows.n_elem, arma::fill::zeros);
  for (arma::uword ii = 0; ii < rows.n_elem; ++ii) {
   const arma::uword row = rows(ii);
   double value = alpha(0u, static_cast<arma::uword>(stick));
   for (arma::uword j = 0; j < annotation.n_cols; ++j) {
    value += annotation(row, j) * alpha(j + 1u,
     static_cast<arma::uword>(stick));
   }
   eta(ii) = value;
  }
  latent[stick].set_size(rows.n_elem);
  for (arma::uword ii = 0; ii < rows.n_elem; ++ii) {
   latent[stick](ii) = st_bayesrc_sample_truncated_normal_std(
    eta(ii), outcome[stick](ii) == 1, generator);
  }
 }
 return latent;
}

inline void st_bayesrc_selection_delta_sweep(
 const arma::mat& annotation,
 const std::vector<arma::uvec>& eligible,
 const std::vector<arma::vec>& latent,
 StBayesRCSelectionState& state,
 std::mt19937& generator,
 const arma::ivec* fixed_delta = nullptr
) {
 const arma::uword annotation_count = annotation.n_cols;
 const arma::uword stick_count = eligible.size();
 std::uniform_real_distribution<double> uniform(0.0, 1.0);
 std::normal_distribution<double> standard_normal(0.0, 1.0);
 for (arma::uword j = 0; j < annotation_count; ++j) {
  double log_odds = std::log(state.pi_a) - std::log1p(-state.pi_a);
  std::vector<StBayesRCSelectionMoments> moments;
  moments.reserve(stick_count);
  for (arma::uword stick = 0; stick < stick_count; ++stick) {
   const arma::uvec& rows = eligible[static_cast<std::size_t>(stick)];
   arma::vec residual = latent[static_cast<std::size_t>(stick)] -
    state.alpha(0u, stick);
   for (arma::uword other = 0; other < annotation_count; ++other) {
    if (other == j) continue;
    residual -= st_bayesrc_selection_column_rows(annotation, other, rows) *
     state.alpha(other + 1u, stick);
   }
   const auto value = st_bayesrc_selection_moments(
    st_bayesrc_selection_column_rows(annotation, j, rows),
    residual, state.tau2(stick));
   moments.push_back(value);
   log_odds += value.log_bf;
  }
  bool included;
  if (fixed_delta != nullptr) {
   if (fixed_delta->n_elem != annotation_count ||
       (*fixed_delta)(j) < 0 || (*fixed_delta)(j) > 1) {
    throw std::invalid_argument("invalid fixed SBayesRC-S delta");
   }
   included = (*fixed_delta)(j) == 1;
  } else {
   const double probability = log_odds >= 0.0 ?
    1.0 / (1.0 + std::exp(-log_odds)) :
    std::exp(log_odds) / (1.0 + std::exp(log_odds));
   included = uniform(generator) < probability;
  }
  state.delta(j) = included ? 1u : 0u;
  for (arma::uword stick = 0; stick < stick_count; ++stick) {
   state.alpha(j + 1u, stick) = included ?
    moments[static_cast<std::size_t>(stick)].mean +
     std::sqrt(moments[static_cast<std::size_t>(stick)].variance) *
      standard_normal(generator) : 0.0;
  }
 }
}

inline void st_bayesrc_selection_blocked_redraw(
 const arma::mat& annotation,
 const std::vector<arma::uvec>& eligible,
 const std::vector<arma::vec>& latent,
 StBayesRCSelectionState& state,
 const arma::vec& intercept_mean,
 const arma::vec& intercept_precision,
 std::mt19937& generator
) {
 const arma::uvec selected = arma::find(state.delta == 1u);
 std::normal_distribution<double> standard_normal(0.0, 1.0);
 for (arma::uword stick = 0; stick < state.alpha.n_cols; ++stick) {
  const arma::uvec& rows = eligible[static_cast<std::size_t>(stick)];
  arma::mat design(rows.n_elem, selected.n_elem + 1u, arma::fill::ones);
  for (arma::uword column = 0; column < selected.n_elem; ++column) {
   design.col(column + 1u) = st_bayesrc_selection_column_rows(
    annotation, selected(column), rows);
  }
  arma::vec prior_precision(selected.n_elem + 1u, arma::fill::zeros);
  prior_precision(0u) = intercept_precision(stick);
  if (selected.n_elem > 0u) {
   prior_precision.subvec(1u, selected.n_elem).fill(1.0 / state.tau2(stick));
  }
  const arma::mat precision = design.t() * design +
   arma::diagmat(prior_precision);
  arma::mat chol_precision;
  if (!arma::chol(chol_precision, precision)) {
   throw std::runtime_error("SBayesRC-S blocked precision is not positive definite");
  }
  arma::vec rhs = design.t() * latent[static_cast<std::size_t>(stick)];
  rhs(0u) += intercept_precision(stick) * intercept_mean(stick);
  const arma::vec mean = arma::solve(
   arma::trimatu(chol_precision),
   arma::solve(arma::trimatl(chol_precision.t()), rhs));
  arma::vec noise(mean.n_elem);
  for (arma::uword i = 0; i < noise.n_elem; ++i) noise(i) = standard_normal(generator);
  const arma::vec draw = mean + arma::solve(arma::trimatu(chol_precision), noise);
  state.alpha.col(stick).zeros();
  state.alpha(0u, stick) = draw(0u);
  for (arma::uword column = 0; column < selected.n_elem; ++column) {
   state.alpha(selected(column) + 1u, stick) = draw(column + 1u);
  }
 }
}

inline void st_bayesrc_selection_update_hyperparameters(
 StBayesRCSelectionState& state,
 const StBayesRCSelectionHyperParameters& hyper,
 std::mt19937& generator
) {
 const arma::vec beta = st_bayesrc_selection_beta_parameters(
  state.delta, hyper.a_pi, hyper.b_pi);
 state.pi_a = st_bayesrc_selection_draw_beta(beta(0u), beta(1u), generator);
 const arma::mat inverse_gamma = st_bayesrc_selection_ig_parameters(
  state.alpha, state.delta, hyper.a_tau, hyper.b_tau);
 for (arma::uword stick = 0; stick < state.tau2.n_elem; ++stick) {
  state.tau2(stick) = st_bayesrc_selection_draw_ig(
   inverse_gamma(stick, 0u), inverse_gamma(stick, 1u), generator);
 }
}

#endif
