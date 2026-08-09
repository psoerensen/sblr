#ifndef SBLR_ST_LOGVAR_ANNOTATION_PRIOR_H
#define SBLR_ST_LOGVAR_ANNOTATION_PRIOR_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#include <armadillo>

namespace sblr {
namespace logvar {

constexpr double kThetaPriorSdV1 = 0.7;

struct EssUpdateDiagnostics {
  std::size_t theta_updates = 0;
  std::size_t likelihood_evaluations = 0;
  std::size_t bracket_contractions = 0;
  std::size_t max_likelihood_evaluations = 0;
  std::size_t max_bracket_contractions = 0;
  double min_log_q = std::numeric_limits<double>::infinity();
  double max_log_q = -std::numeric_limits<double>::infinity();

  void record_update(std::size_t evaluations, std::size_t contractions) {
    ++theta_updates;
    likelihood_evaluations += evaluations;
    bracket_contractions += contractions;
    max_likelihood_evaluations =
      std::max(max_likelihood_evaluations, evaluations);
    max_bracket_contractions =
      std::max(max_bracket_contractions, contractions);
  }

  void record_eta(const arma::vec& eta) {
    if (eta.n_elem == 0) return;
    min_log_q = std::min(min_log_q, eta.min());
    max_log_q = std::max(max_log_q, eta.max());
  }
};

inline void validate_annotation_inputs(
  const arma::mat& annotation,
  const arma::vec& theta,
  std::size_t marker_count
) {
  if (annotation.n_rows != marker_count) {
    throw std::invalid_argument(
      "log-variance annotation row count must match the marker count.");
  }
  if (annotation.n_cols == 0) {
    throw std::invalid_argument(
      "log-variance annotation design must contain at least one column.");
  }
  if (theta.n_elem != annotation.n_cols) {
    throw std::invalid_argument(
      "log-variance theta length must match the annotation column count.");
  }
  if (!annotation.is_finite() || !theta.is_finite()) {
    throw std::invalid_argument(
      "log-variance annotation design and theta must be finite.");
  }
}

inline arma::vec calculate_eta(
  const arma::mat& annotation,
  const arma::vec& theta
) {
  validate_annotation_inputs(annotation, theta, annotation.n_rows);
  arma::vec eta = annotation * theta;
  if (!eta.is_finite()) {
    throw std::runtime_error(
      "log-variance eta contains a non-finite value.");
  }
  return eta;
}

inline arma::vec calculate_prior_scale_from_eta(
  const arma::vec& eta,
  EssUpdateDiagnostics* diagnostics = nullptr
) {
  if (!eta.is_finite()) {
    throw std::runtime_error(
      "log-variance eta contains a non-finite value.");
  }
  // Match the frozen R oracle's explicit log-scale validity window. This is a
  // rejection/error guard, not a clamp.
  constexpr double kMaximumAbsoluteEta = 700.0;
  arma::vec scale(eta.n_elem);
  for (arma::uword j = 0; j < eta.n_elem; ++j) {
    const double value = eta(j);
    if (value > kMaximumAbsoluteEta || value < -kMaximumAbsoluteEta) {
      throw std::overflow_error(
        "log-variance exp(eta) is outside the positive finite double range.");
    }
    scale(j) = std::exp(value);
    if (!std::isfinite(scale(j)) || scale(j) <= 0.0) {
      throw std::overflow_error(
        "log-variance prior scale is not positive and finite.");
    }
  }
  if (diagnostics != nullptr) diagnostics->record_eta(eta);
  return scale;
}

inline arma::vec calculate_prior_scale(
  const arma::mat& annotation,
  const arma::vec& theta,
  EssUpdateDiagnostics* diagnostics = nullptr
) {
  return calculate_prior_scale_from_eta(
    calculate_eta(annotation, theta), diagnostics);
}

inline double scaled_square(
  double effect,
  double log_denominator,
  const char* model
) {
  if (!std::isfinite(effect) || !std::isfinite(log_denominator)) {
    throw std::runtime_error(
      std::string(model) + " theta likelihood received a non-finite state.");
  }
  if (effect == 0.0) return 0.0;
  const double log_value = 2.0 * std::log(std::abs(effect)) - log_denominator;
  const double log_max = std::log(std::numeric_limits<double>::max());
  if (!std::isfinite(log_value) || log_value > log_max) {
    throw std::overflow_error(
      std::string(model) + " theta likelihood quadratic term overflowed.");
  }
  const double value = std::exp(log_value);
  if (!std::isfinite(value)) {
    throw std::overflow_error(
      std::string(model) + " theta likelihood is non-finite.");
  }
  return value;
}

inline double theta_log_likelihood_bayesc(
  const arma::vec& theta,
  const arma::mat& annotation,
  const arma::vec& effect,
  const arma::ivec& state,
  double marker_variance
) {
  validate_annotation_inputs(annotation, theta, effect.n_elem);
  if (state.n_elem != effect.n_elem) {
    throw std::invalid_argument(
      "BayesC-LV state length must match the marker count.");
  }
  if (!effect.is_finite() || !std::isfinite(marker_variance) ||
      marker_variance <= 0.0) {
    throw std::invalid_argument(
      "BayesC-LV effects must be finite and marker variance positive finite.");
  }
  const arma::vec eta = calculate_eta(annotation, theta);
  const double log_vb = std::log(marker_variance);
  double result = 0.0;
  for (arma::uword j = 0; j < effect.n_elem; ++j) {
    if (state(j) <= 0) continue;
    if (eta(j) > 700.0 || eta(j) < -700.0) {
      throw std::overflow_error(
        "BayesC-LV log prior scale is outside the validated [-700, 700] range.");
    }
    result -= 0.5 * (eta(j) + scaled_square(
      effect(j), log_vb + eta(j), "BayesC-LV"));
    if (!std::isfinite(result)) {
      throw std::overflow_error(
        "BayesC-LV theta likelihood accumulated a non-finite value.");
    }
  }
  return result;
}

inline double theta_log_likelihood_bayesr(
  const arma::vec& theta,
  const arma::mat& annotation,
  const arma::vec& effect,
  const arma::ivec& component,
  double marker_variance,
  const arma::vec& gamma
) {
  validate_annotation_inputs(annotation, theta, effect.n_elem);
  if (component.n_elem != effect.n_elem || gamma.n_elem < 2) {
    throw std::invalid_argument(
      "BayesR-LV component or gamma dimensions are invalid.");
  }
  if (!effect.is_finite() || !gamma.is_finite() ||
      !std::isfinite(marker_variance) || marker_variance <= 0.0) {
    throw std::invalid_argument(
      "BayesR-LV effects/gamma must be finite and marker variance positive finite.");
  }
  const arma::vec eta = calculate_eta(annotation, theta);
  const double log_vb = std::log(marker_variance);
  double result = 0.0;
  for (arma::uword j = 0; j < effect.n_elem; ++j) {
    const int k = component(j);
    if (k == 0) continue;
    if (eta(j) > 700.0 || eta(j) < -700.0) {
      throw std::overflow_error(
        "BayesR-LV log prior scale is outside the validated [-700, 700] range.");
    }
    if (k < 0 || static_cast<arma::uword>(k) >= gamma.n_elem ||
        gamma(static_cast<arma::uword>(k)) <= 0.0) {
      throw std::invalid_argument(
        "BayesR-LV active component has no positive gamma scale.");
    }
    result -= 0.5 * (eta(j) + scaled_square(
      effect(j), log_vb + std::log(gamma(static_cast<arma::uword>(k))) +
        eta(j), "BayesR-LV"));
    if (!std::isfinite(result)) {
      throw std::overflow_error(
        "BayesR-LV theta likelihood accumulated a non-finite value.");
    }
  }
  return result;
}

template <class Generator>
inline arma::vec draw_theta_prior(
  std::size_t dimension,
  double prior_sd,
  Generator& generator
) {
  if (dimension == 0 || !std::isfinite(prior_sd) || prior_sd <= 0.0) {
    throw std::invalid_argument(
      "theta prior requires positive dimension and positive finite SD.");
  }
  std::normal_distribution<double> normal(0.0, prior_sd);
  arma::vec result(dimension);
  for (arma::uword j = 0; j < result.n_elem; ++j) result(j) = normal(generator);
  if (!result.is_finite()) {
    throw std::runtime_error("theta prior draw produced a non-finite value.");
  }
  return result;
}

template <class Generator, class LogLikelihood>
inline arma::vec elliptical_slice_update(
  const arma::vec& current,
  double prior_sd,
  bool active_set_empty,
  LogLikelihood log_likelihood,
  Generator& generator,
  EssUpdateDiagnostics& diagnostics,
  std::size_t maximum_evaluations = 10000
) {
  if (current.n_elem == 0 || !current.is_finite() ||
      !std::isfinite(prior_sd) || prior_sd <= 0.0) {
    throw std::invalid_argument(
      "elliptical slice inputs require finite theta and positive finite prior SD.");
  }
  if (active_set_empty) {
    arma::vec result = draw_theta_prior(current.n_elem, prior_sd, generator);
    diagnostics.record_update(0, 0);
    return result;
  }

  const double current_ll = log_likelihood(current);
  if (!std::isfinite(current_ll)) {
    throw std::runtime_error(
      "current theta has a non-finite conditional log likelihood.");
  }

  arma::vec direction = draw_theta_prior(current.n_elem, prior_sd, generator);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  const double u_threshold = std::max(
    uniform(generator), std::numeric_limits<double>::min());
  const double threshold = current_ll + std::log(u_threshold);
  const double two_pi = 2.0 * std::acos(-1.0);
  double angle = two_pi * uniform(generator);
  double lower = angle - two_pi;
  double upper = angle;
  std::size_t evaluations = 0;
  std::size_t contractions = 0;

  while (evaluations < maximum_evaluations) {
    arma::vec proposal = current * std::cos(angle) +
      direction * std::sin(angle);
    double proposal_ll = -std::numeric_limits<double>::infinity();
    try {
      proposal_ll = log_likelihood(proposal);
    } catch (const std::overflow_error&) {
      proposal_ll = -std::numeric_limits<double>::infinity();
    }
    ++evaluations;
    if (std::isfinite(proposal_ll) && proposal_ll > threshold) {
      diagnostics.record_update(evaluations, contractions);
      return proposal;
    }
    ++contractions;
    if (angle < 0.0) lower = angle;
    else upper = angle;
    angle = lower + (upper - lower) * uniform(generator);
  }
  throw std::runtime_error(
    "elliptical slice sampling exceeded the likelihood-evaluation guard.");
}

}  // namespace logvar
}  // namespace sblr

#endif
